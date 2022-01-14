# Because "melanoma_2_cell_lines_af.txt" is a big file that cannot be uploaded 
# to github, we provided allele frequency matrixes of two cell lines with 
# "451Lu_af_var.txt" and "A357_af_var.txt". You will get "melanoma_2_cell_lines_af.txt", 
# which can be used as input of LINEAGE, by putting two matrix files and this 
# R script in the same directory and running command “Rscipt af.R”.

library(tidyr)

a=read.table('melanoma_451Lu_af_var.txt',header=T,row.names=1)
b=read.table('melanoma_A357_af_var.txt',header=T,row.names=1)

af=merge(a,b,by='var',all=T)
af[is.na(af)]=0
af=separate(af,var,c('pos','Allele'),sep='_r')
af=separate(af,Allele,c('refAllele','altAllele'),sep='_')
rownames(af)=paste(af$pos,af$altAllele,sep='_')
af=subset(af,select=-pos)

write.table(af,'melanoma_2_cell_lines_af.txt',col.names=T,row.names=T,sep='\t')

