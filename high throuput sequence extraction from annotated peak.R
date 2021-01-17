for (i in as.data.frame(annotatedPeak_ecr_6aph$"symbol")){
bed6aphstart<-c(as.data.frame(annotatedPeak_ecr_6aph)$"start"[as.data.frame(annotatedPeak_ecr_6aph$"symbol")==i])
bed6aphend<-c(as.data.frame(annotatedPeak_ecr_6aph)$"end"[as.data.frame(annotatedPeak_ecr_6aph$"symbol")==i])
bed6aphchr<-c(as.character(as.data.frame(annotatedPeak_ecr_6aph)$"seqnames"[as.data.frame(annotatedPeak_ecr_6aph$"symbol")==i]))
#print(i)
}


bed6aphstart=as.numeric(na.omit(bed6aphstart))
bed6aphend=as.numeric(na.omit(bed6aphend))
bed6aphchr=as.character(na.omit(bed6aphchr))

ecr6aph=GRanges(seqnames=as.factor(bed6aphchr),
	ranges=IRanges(start=as.numeric(bed6aphstart),
		end=as.numeric(bed6aphend)))


ecr_seq_fasta6aph=getSeq(Dmelanogaster, ecr6aph[1:length(ecr6aph)])
names(ecr_seq_fasta6aph)=as.character(ecr6aph[1:length(ecr6aph)])
#names(ecr_seq_fasta6aph)=as.factor(bed6aphchr)
ecr_seq_fasta6aph=unique(ecr_seq_fasta6aph)
#writeXStringSet(ecr_seq_fasta6aph[1:length(ecr6aph)], file="ecr_seq6aph.fasta")
#writeXStringSet(ecr_seq_fasta6aph[1:length(ecr_seq_fasta6aph)], file="ecr_seq6aph2.fasta")
writeXStringSet(ecr_seq_fasta6aph, file="ecr_seq6aph2.fasta")

motif<-c("AGTBG")

seqdna<-vmatchPattern(motif,ecr_seq_fasta6aph,fixed=F)
length(unique(unlist(seqdna)))
df<-as.data.frame(unlist(seqdna))
length(unique(df$"names"))

#writeXStringSet(seqdna, file="seqdnatest.fasta",append=T)
writeXStringSet(seqdna, file="seqdnatest.fasta")

start(seqdna[length(start(seqdna))>0])

fastaFile <- readDNAStringSet("my.fas")