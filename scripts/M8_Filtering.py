#!/usr/bin/env python3
# -*- coding: latin-1 -*-
from collections import defaultdict,Counter 
import argparse
import json
import sys 

def generate_all_annotation(M8_file) :
  def Parse_line(Line) :
    querry,subject,qstart,qend,qlen,sstart,send,slen,AL,Pid,evalue,bitscore=Line.rstrip().split()
    if int(qend)<int(qstart) :
        qstart,qend=qend,qstart
    Length=float(qend)-float(qstart)
    return querry,[querry,subject,float(evalue),int(qstart),int(qend),float(qend)-float(qstart),Line]
  # start here 
  Handle=open(M8_file) 
  querry,Annotation=Parse_line(next(Handle))
  current_querry=querry
  List_annotation=[Annotation]
  for line in Handle :
    querry,Annotation=Parse_line(line)
    if querry==current_querry :
      List_annotation.append(Annotation)
    else :
      yield current_querry,sorted(List_annotation,key=lambda x:x[2])
      current_querry=querry
      List_annotation=[Annotation]

def generate_filtered_Annotation(Annotation_generator) :
  for querry,List_annotation in Annotation_generator :
    # start by the smallest evalue use it to remove any highest overlaping evalue 
    # and go to the next entry not removed precedently, it will be the next best evalue 
    # and you can prune other leftover entries.
    List_Filtered_annotation=[]
    Index_del=set([])
    for index1,annotation1 in enumerate(List_annotation) :
      if index1 in Index_del :
        continue
      List_Filtered_annotation.append(annotation1[-1])
      for index2,annotation2 in enumerate(List_annotation) :
        if (index2 in Index_del)|(index1==index2) :
          continue
        _,_,evalue1,start1,end1,length1,_=annotation1
        _,_,evalue2,start2,end2,length2,_=annotation2
        if start1<=start2 :
          Overlap=(end1-start2)*2/(length1+length2)>0.1
          Inclusion=end2<=end1
        else :
          Overlap=(end2-start1)*2/(length1+length2)>0.1
          Inclusion=end1<=end2
        if Overlap|Inclusion :
          Index_del|={index2}
    yield sorted(List_Filtered_annotation,key=lambda x:int(x.split()[2]))

# min_Bitscore,max_Evalue,min_Pid,min_ref_pid,min_coverage,min_Query_coverage=(0,1e-5,0,0,0,0)
def main(m8_file,Database_file,min_Bitscore,max_Evalue,min_Pid,min_ref_pid,min_coverage,min_Query_coverage) :
	Header=[]
	if Database_file :
		# ----- get annotation information ---------------
		Handle=open(Database_file)
		Header=next(Handle).rstrip().split("\t")
		# ineficient if annotation file is big, but should still be overall faster.
		Dico_Ref_annotation={line.rstrip().split()[0]:line.rstrip().split()[1:] for line in open(Database_file)}
		Handle.close()
	# ---- output ----
	print("\t".join(["Query","Subject","Bitscore","PID","Subject_Pid","Coverage","Query_coverage"]+Header[2:]))
	for List_annotation in generate_filtered_Annotation(generate_all_annotation(m8_file)) :
		for annotation in List_annotation :
			#querry,subject,qstart,qend,qlen,sstart,send,slen,AL,Pid,evalue,bitscore
			Splitline=annotation.rstrip().split("\t")
			#--------Alignement info--------
			Len_AL=float(Splitline[8])
			PID=float(Splitline[9])/100 # must divide by 100. to get a real percentage
			Evalue=float(Splitline[10])
			Bitscore=float(Splitline[11])
			#--------querry info--------
			Query=Splitline[0]
			Len_query=float(Splitline[4])
			Query_coverage=Len_AL/Len_query
			#--------subject info--------
			Subject=Splitline[1]
			Len_subject=float(Splitline[7])
			Coverage=Len_AL/Len_subject
			#-------- usefull info --------
	 		## compute number of Match over length of reference, Subject_Pid
			Subject_Pid=Coverage*PID
			#-------- select only hit within criteria and output --------
			if (Bitscore>=min_Bitscore)&(max_Evalue>=Evalue)&(PID>=min_Pid)&(Subject_Pid>=min_ref_pid)&(Coverage>=min_coverage)&(Query_coverage>=min_Query_coverage) :
				if Database_file :
					try :
						Newname=Dico_Ref_annotation[Subject][0]
						Additional_annotation=Dico_Ref_annotation[Subject][1:]
					except :
						print("\t".join([Query,Subject,Bitscore,"{:.3f}".format(PID),"{:.3f}".format(Subject_Pid),"{:.3f}".format(Coverage),"{:.3f}".format(Query_coverage)]), file=sys.stderr)
						continue
				else :
					Newname=Subject
					Additional_annotation=[]
				Result= [Query,Newname,Bitscore,"{:.3f}".format(PID),"{:.3f}".format(Subject_Pid),"{:.3f}".format(Coverage),"{:.3f}".format(Query_coverage)]+Additional_annotation
				print("\t".join(map(str,Result)))
				
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("m8_file", help=" m8 output, columns should be 'qseqid sseqid qstart qend qlen sstart send slen length pident evalue bitscore'")  
	parser.add_argument("-D", help="Database, should be a .tsv file, should have one line header, column should be : sseqid,new name,additional annotation1,additional annotation2,additional annotation...",default=0)
	parser.add_argument("-B",help="minimum value cutoff for Bitscore",default=0)
	parser.add_argument("-E",help="maximum value for Evalue",default=1e-5)
	parser.add_argument("-P",help="cutoff for PID, between 0 and 1",default=0)
	parser.add_argument("-R",help="Subject pid : minimum (Subject_Coverage x percentage_identity), between 0 and 1",default=0)
	parser.add_argument("-C",help="Subject Coverage cutoff : minimum percentage of the subject the Alignment does cover, between 0 and 1 ",default=0)
	parser.add_argument("-Q",help="Querry Coverage cutoff : minimum percentage of the Querry the Alignment does cover, between 0 and 1 ",default=0)
	args = parser.parse_args()
	m8_file=args.m8_file
	Database_file=args.D
	min_Bitscore=float(args.B)
	max_Evalue=float(args.E)
	min_PID=float(args.P)
	min_ref_pid=float(args.R)
	min_coverage=float(args.C)
	min_Query_coverage=float(args.Q)
	try :
		handle=open(Database_file)
		handle.close()
	except :
		print("no valid database file given, it can't be opened for some reason",file=sys.stderr)
		Database_file=""
	main(m8_file,Database_file,min_Bitscore,max_Evalue,min_PID,min_ref_pid,min_coverage,min_Query_coverage)
