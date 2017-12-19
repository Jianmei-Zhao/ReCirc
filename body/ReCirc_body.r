######read blast results of exon array probes with protein coding gene and circRNA body sequences separately#########

HuEx_pc<-read.table("HuEx2pc_19.txt",sep="\t",header=T,stringsAsFactors=FALSE)
HuEx_cir<-read.table("HuEx2circ.txt",sep="\t",header=T,stringsAsFactors=FALSE)

###########remove probes that match to both coding gene and circRNA#####################################################

     Overlap_cir_pc<-intersect(HuEx_cir[[1]],HuEx_pc[[1]]);
     PC_Del<-HuEx_pc[-which(HuEx_pc[[1]]%in%Overlap_cir_pc==TRUE),];
	 Cir_Del<-HuEx_cir[-which(HuEx_cir[[1]]%in%Overlap_cir_pc==TRUE),];
	 
#############################filter probes that match to more than one coding gene or circRNA##############################################
FilterProbeMultiTar<-function(ProbeMap){
	UniqueSet<-unlist(strsplit(unique(paste(ProbeMap[[1]],ProbeMap[[2]],sep=">")),split=">"));
	ProbeMap<-as.data.frame(cbind(UniqueSet[(1:(length(UniqueSet)/2))*2-1],UniqueSet[(1:(length(UniqueSet)/2))*2]),stringsAsFactors=FALSE);
	Probes<-unique(ProbeMap[[1]]);
	Probedelete<-c();
	for(i in 1:length(Probes)){ 
		if(length(which(ProbeMap[[1]]==Probes[i]))>1) Probedelete<-c(Probedelete,Probes[i]);
		if(i %% 5000 == 0) print(i);
	}
	locations<-which((ProbeMap[[1]]%in%Probedelete)==TRUE);
	probe2tar<-ProbeMap[-locations,];
	return(probe2tar);
}
probe2gene<-FilterProbeMultiTar(PC_Del);
probe2cir<-FilterProbeMultiTar(Cir_Del);

########################delete genes and circRNAs whose matched probes were less than 3 ########################################
Credible<-function(ProbeMap,n){
	tar<-unique(ProbeMap[[2]]);
	More3<-c();
	for(i in 1:length(tar)){
		locations<-which(ProbeMap[[2]]==tar[i]);
		if(length(locations)>=n) More3<-c(More3,locations);
		if(i %% 5000 == 0) print(i);
	}
	Probe2Tar<-ProbeMap[More3,];
	return(Probe2Tar);
}
Probe2Gene<-Credible(probe2gene,3);
write.table(Probe2Gene,"Probe2Gene.txt",sep="\t",row.names=F)
Probe2Cir<-Credible(probe2cir,3); 
write.table(Probe2Cir,"Probe2Cir.txt",sep="\t",row.names=F)	 

#################################filt probes matching circRNAs perfectly and uniquely with protein coding genes stored in UCSC and NCBI database########
HuEx_ucsc<-read.table("HuEx_ucsc.txt",sep="\t",stringsAsFactors=F)
HuEx_ncbi<-read.table("HuEx_NM_ncbi.txt",sep="\t",stringsAsFactors=F)
PR2Cir<-read.table("Probe2Cir.txt",sep="\t",stringsAsFactors=F)
a<-intersect(as.character(PR2Cir[,1]),as.character(HuEx_ncbi[,1]))
PR2Cir1<-PR2Cir[-which(PR2Cir[,1]%in%a),]
b<-intersect(as.character(PR2Cir1[,1]),as.character(HuEx_ucsc[,1]))
PR2Cir2<-PR2Cir1[-which(PR2Cir1[,1]%in%b),]
d<-as.data.frame(table(PR2Cir2[,2]))
e<-as.character(d[-which(d[,2]<3),1])
PR2Cir3<-PR2Cir2[which(PR2Cir2[,2]%in%e),]
write.table(PR2Cir3,"flit_Probe2Cir.txt",sep="\t",row.names=F,quote=F)
	 	 
###########################################################convert prboes to corrrsponding probe sets###########################################	 
probe2set<-read.table("probe2set.txt",header=T,sep="\t",stringsAsFactors=F)
Probe2Cir<-read.table("flit_Probe2Cir.txt",header=T,stringsAsFactors=F,sep="\t")
Cir_loc<-match(Probe2Cir[[1]],probe2set[[1]])
Probe_set2<-probe2set[Cir_loc,3]
Probe_set_Cir<-cbind(Probe2Cir[[1]],Probe_set2,Probe2Cir[[2]])
colnames(Probe_set_Cir)<-c("probe","probeset","circRNA")
Probeset_Cir<-unique(cbind(Probe_set_Cir[,2],Probe_set_Cir[,3]))
colnames(Probeset_Cir)<-c("probeset","circRNA")
write.table(Probeset_Cir,"filter_PRS2Cir.txt",sep="\t",row.names=F,quote=F)
