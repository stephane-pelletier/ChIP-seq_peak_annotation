bed6aphstart<-c(as.data.frame(annotatedPeak_ecr_6aph)$"start"[as.data.frame(annotatedPeak_ecr_6aph$"symbol")=="Eip75B"])
bed6aphend<-c(as.data.frame(annotatedPeak_ecr_6aph)$"end"[as.data.frame(annotatedPeak_ecr_6aph$"symbol")=="Eip75B"])
#bed6aphchr<-c(as.data.frame(annotatedPeak_ecr_6aph)$"seqnames"[as.data.frame(annotatedPeak_ecr_6aph$"symbol")=="Eip75B"])

bed6aphstart=as.numeric(na.omit(bed6aphstart))
bed6aphend=as.numeric(na.omit(bed6aphend))

bed6aphchr<-c(rep("chr3L",length(bed6aphstart)))

ecr6aph=GRanges(seqnames=as.factor(bed6aphchr),
	ranges=IRanges(start=as.numeric(bed6aphstart),
		end=as.numeric(bed6aphend)))