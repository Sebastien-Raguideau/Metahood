#!/usr/bin/env python
import os 
from Bio.SeqIO.FastaIO import *
from subprocess import Popen, PIPE
import argparse
import time 

# must change batch strategy to count number of characters not number of fasta entry. Otherwise it can be pretty stupid if part of entries are significantly bigger than others.

def Cant_launch_new_task(Number_task):
	# return false when new task can be launched
	process = Popen(['ps aux|grep '+os.getlogin()+'|grep " [p]rodigal" | wc -l '], stdout=PIPE, stderr=PIPE,shell=True)
	NB = int(process.communicate()[0].rstrip())
	return NB>=Number_task

def main(Fasta_file,Temp_location,Number_task,Splits,output_directory) :
	process = Popen(['grep -c ">" ' +Fasta_file ], stdout=PIPE, stderr=PIPE,shell=True)
	nb_line = int(process.communicate()[0].split()[0])
	os.system("mkdir -p " + Temp_location)
	Nb_line_batch=nb_line/float(Splits)
	Batch_num=0
	## Divide the initial fasta file in batchs and launch prodigal
	CWD=os.getcwd()
	Fasta_Directory="/".join(Fasta_file.split("/")[:-1])
	if Fasta_Directory :
		os.chdir(Fasta_Directory)
		AFD=os.getcwd()+"/"+Fasta_file.split("/")[-1]
		os.chdir(CWD)
	else :
		AFD=CWD+"/"+Fasta_file.split("/")[-1]
	os.chdir(Temp_location)
	file_name="Batch_0"
	Handle=open(file_name,"w")
	for (index,(header,seq)) in enumerate(SimpleFastaParser(open(AFD))) :
		if index > Nb_line_batch*(Batch_num+1) :
			if Batch_num==Splits-1 :
				Handle.write(">"+header+"\n"+seq+"\n")	
			else :
				Handle.write(">"+header+"\n"+seq+"\n")
				Handle.close()
				while Cant_launch_new_task(Number_task) :
					time.sleep(10)
				os.system("prodigal -i "+file_name+" -a "+file_name+".faa -d "+file_name+".fna -f gff -p meta -o "+file_name+".gff > "+file_name+"_prodigal.out 2>&1 &" )
				Batch_num+=1
				file_name="Batch_"+str(Batch_num)
				Handle=open("Batch_"+str(Batch_num),"w")
		else :
			Handle.write(">"+header+"\n"+seq+"\n")
	Handle.close()
	os.system("nohup prodigal -i "+file_name+" -a "+file_name+".faa -d "+file_name+".fna -f gff -p meta -o "+file_name+".gff > "+file_name+"_prodigal.out 2>&1 &" )
	process = Popen(['ps aux|grep '+os.getlogin()+'|grep " [p]rodigal" | wc -l '], stdout=PIPE, stderr=PIPE,shell=True)
	nb_prodigal_batch = int(process.communicate()[0].rstrip())
	while nb_prodigal_batch!=0 :
		time.sleep(60)
		process = Popen(['ps aux|grep '+os.getlogin()+'|grep "[p]rodigal" | grep -v Parallel_prodigal.py | wc -l'], stdout=PIPE, stderr=PIPE,shell=True)
		nb_prodigal_batch = int(process.communicate()[0].rstrip())
	Fasta_File_name=Fasta_file.split("/")[-1]
	Handle_faa=open(CWD+"/"+output_directory+"/"+Fasta_File_name.replace(".fa",".faa"),"w")
	Handle_fna=open(CWD+"/"+output_directory+"/"+Fasta_File_name.replace(".fa",".fna"),"w")
	Handle_gff=open(CWD+"/"+output_directory+"/"+Fasta_File_name.replace(".fa",".gff"),"w")
	Handle_gff.write(open("Batch_0.gff").readline())
	for nb_batch in range(Splits) :
		Handle_faa.write(open("Batch_"+str(nb_batch)+".faa").read())
		Handle_fna.write(open("Batch_"+str(nb_batch)+".fna").read())
		tmp=open("Batch_"+str(nb_batch)+".gff")
		_=tmp.readline()
		Handle_gff.write(tmp.read())
	Handle_faa.close()
	Handle_fna.close()
	Handle_gff.close()
	os.chdir(CWD)
	time.sleep(10)
	# os.system("rm -rf "+Temp_location[:-1]+"/*.gff" )
	# os.system("rm -rf "+Temp_location[:-1]+"/*.fna" )
	# os.system("rm -rf "+Temp_location[:-1]+"/*.out" )


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("n", help="Number of prodigal tasks")
	parser.add_argument("-s", help="Number of splits, default is same as number of tasks")
	parser.add_argument("-T", help="Temp file location ",default="./Split_Annotation")
	parser.add_argument("f", help="Fasta_file")
	parser.add_argument("-o", help="output Directory",default='./')
	args = parser.parse_args()
	Fasta_file=args.f
	Number_task=int(args.n)
	Temp_location=args.T
	Temp_location=Temp_location+(not Temp_location[-1]=="/")*"/"
	if args.s :
		Splits=int(args.s)
	else:
		Splits=Number_task
	output_directory=args.o
	main(Fasta_file,Temp_location,Number_task,Splits,output_directory)




