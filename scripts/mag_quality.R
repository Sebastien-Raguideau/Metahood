library(readr)
library(ggplot2)
library(RColorBrewer)


# all data summary
datSummary <- readr::read_tsv("/ei/.project-scratch/7/7a8b76af-019e-4731-8856-2ddea6f071a9/CD_TREATS/MAGs_CDTREATS/results/mag_to_dmags_summary.tsv", col_names=TRUE)  #MAGs with taxonomy and 
binwidth = c(2.5,.5)
p<-ggplot(datSummary, aes(completion,contamination))+ geom_bin2d(binwidth = binwidth) + stat_bin2d(geom="text",aes(label = ..count..,fontface = "bold"),size=3,binwidth = binwidth)+ scale_fill_gradientn(colours=rf(32), trans="log10") +ylim(-0.001,10.001) + xlim(49,101)
pdf("test_qual.pdf")
print(p)
dev.off()

