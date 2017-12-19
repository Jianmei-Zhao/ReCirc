creactecirdb<-function(pwd)
			 {
			     mysp<-function(x)
					   {
						   null<-c()
							for(i in 1:length(x))
							  {
								 null<-paste0(null,x[i])
							  }
							null
					   }
				mystr<-function(m)
					{
						 if(length(m)>=120)
							{
							  result<-paste0(mysp(m[(length(m)-23):length(m)]),mysp(m[1:24]))
							}
						  else
							{
							   result<-NA
							}
						   result
					 }

				mfile<-readLines(pwd)#��ȡ�ļ�
				location<-grepl(">",mfile)#�ҵ�̽������
				probe<-which(location==T)#�õ�̽���������ڵ�λ��
				gene<-mfile[location]##��ȡ̽������
				star<-as.numeric( probe)+1###��ȡ���е���ʼ��
				end<-c(as.numeric(probe[-1])-1,length(mfile))##
				null<-list()
				gc()
				for(i in 1:length(star))
				   {
					  null[[i]]<-unlist(strsplit(c(mfile[star[i]:end[i]]),""))
				   }
				mypg<-unlist(lapply(null,mystr))
				gc()
				myresult<-matrix("",2*length(gene),1)
				myresult[seq(1,2*length(gene),2),]<-gene
				myresult[seq(2,2*length(gene),2),]<-mypg
				myresult
			}

mfile<-creactecirdb("circRNA.fasta")
location<-which(mfile%in%NA)
newfile<-mfile[-c(location,location-1)]
mtrfile<-as.matrix(newfile)
write.table(mtrfile,"cirdb2.fa",row.names=F,col.names=F,quote=F)