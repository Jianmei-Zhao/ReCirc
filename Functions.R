DelUncertainProbe <- function(Circ_data,Pc_data)
			 {
				 Overlap_cir_pc <- intersect(Circ_data[[1]],Pc_data[[1]]);
				 PC_Del <- Pc_data[-which(Pc_data[[1]]%in%Overlap_cir_pc==TRUE),];
				 Cir_Del <- Circ_data[-which(Circ_data[[1]]%in%Overlap_cir_pc==TRUE),];
				 Del_data <- list(PC_Del,Cir_Del);
				 names(Del_data)<-c('PC_Del','Cir_Del');
				 return(Del_data);
			 } # remove probes that match to both coding gene and circRNA

FilterProbeMultiTar <- function(ProbeMap)
			 {
				 UniqueSet <- unlist(strsplit(unique(paste(ProbeMap[[1]],ProbeMap[[2]],sep=">")),split=">"));
				 ProbeMap <- as.data.frame(cbind(UniqueSet[(1:(length(UniqueSet)/2))*2-1],UniqueSet[(1:(length(UniqueSet)/2))*2]),stringsAsFactors=FALSE);
				 Probes <- unique(ProbeMap[[1]]);
				 Probedelete <- c();
					 for(i in 1:length(Probes))
						 { 
							 if(length(which(ProbeMap[[1]]==Probes[i]))>1) Probedelete <- c(Probedelete,Probes[i]);
							 if(i %% 5000 == 0) print(i);
						 }
				 locations <- which((ProbeMap[[1]]%in%Probedelete)==TRUE);
				 probe2tar <- ProbeMap[-locations,];
				 return(probe2tar);
			 } # filter probes that match to more than one coding gene or circRNA

Credible <- function(ProbeMap,n)
			 {
				 tar <- unique(ProbeMap[[2]]);
				 More3 <- c()
				 for(i in 1:length(tar))
					 {
						 locations <- which(ProbeMap[[2]]==tar[i]);
						 if(length(locations)>=n) More3 <- c(More3,locations);
						 if(i %% 5000 == 0) 
						 (i);
					 }
				 Probe2Tar <- ProbeMap[More3,];
				 return(Probe2Tar);
			 } # delete genes and circRNAs whose matched probes were less than 3; convert probes in probe-gene pairs to corresponding probeset

Affy_Body_OtherFilter <- function(Probe2Cir,ncbi,ucsc,n=3)
			 {
				 inter1 <- intersect(as.character(Probe2Cir[,1]),as.character(ncbi[,1]))
				 Probe2Cir1 <- Probe2Cir[-which(Probe2Cir[,1]%in%inter1),]
				 inter2 <- intersect(as.character(Probe2Cir1[,1]),as.character(ucsc[,1]))
				 Probe2Cir2 <- Probe2Cir1[-which(Probe2Cir1[,1]%in%inter2),]
				 probe_number <- as.data.frame(table(Probe2Cir2[,2]))
				 probe_result <- as.character(probe_number[-which(probe_number[,2]<n),1])
				 Probe2Cir3 <- Probe2Cir2[which(Probe2Cir2[,2]%in%probe_result),]
				 return(Probe2Cir3);
			 } # filt probes matching circRNAs perfectly and uniquely with protein coding genes stored in UCSC and NCBI database
			 
Affy_Junc_OtherFilter <- function(Probe2Cir,ncbi,ucsc)
			 {
				 Probe2Circ_ucsc <- Probe2Cir[-match(unique(as.character(ucsc[,1])),as.character(Probe2Cir[,1])),]
				 Probe2Circ_ucsc_ncbi <- Probe2Circ_ucsc[-na.omit(match(unique(as.character(ncbi[,1])),as.character(Probe2Circ_ucsc[,1]))),]	
				 return(Probe2Circ_ucsc_ncbi);
			 } # filt probes matching circRNAs perfectly and uniquely with protein coding genes stored in UCSC and NCBI database
			 
Agilent_OtherFilter <- function(ncbi_data,ucsc_data,PR2Cir)
			 {
				 ncbi_data_PR2Cir <- merge(PR2Cir,ncbi_data,by.x=1,by.y=1)
				 PR2Cir <- PR2Cir[-na.omit(match(unique(ncbi_data_PR2Cir[,1]),PR2Cir[,1])),]
				 ucsc_data_PR2Cir <- merge(PR2Cir,ucsc_data,by.x=1,by.y=1)
				 PR2Cir <- PR2Cir[-na.omit(match(unique(ucsc_data_PR2Cir[,1]),PR2Cir[,1])),]
				 return(PR2Cir)
			 } # filt probes matching circRNAs perfectly and uniquely with protein coding genes stored in UCSC and NCBI database

ConvertProbe <- function(Probe2target,probe2set)
			 {
				 Target_loc <- match(Probe2target[[1]],probe2set[[1]])
				 Probe_set2 <- probe2set[Target_loc,3]
				 Probe_set_Target <- cbind(as.character(Probe2target[[1]]),Probe_set2,as.character(Probe2target[[2]]))
				 colnames(Probe_set_Target)<-c("probe","probeset","Target")
				 Probeset_Tar <- unique(cbind(Probe_set_Target[,2],Probe_set_Target[,3]))
				 colnames(Probeset_Tar)<-c("probeset","Target")
				 return(Probeset_Tar)
			 } # convert prboes in probe-circRNA pairs to corrrsponding probe sets