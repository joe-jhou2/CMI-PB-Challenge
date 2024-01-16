setwd("C:/Users/Joseph/Desktop")
aa=readLines("gene.txt")
ftable=c()
for(i in 2:length(aa)){
	splited=strsplit(aa[i],split="\t")[[1]]
	splited02=splited[splited!=""]
	for(j in 6:length(splited02)){
		ftable=rbind(ftable,c(splited02[1],splited02[2],splited02[3],splited02[j]))
	}
}

dim(ftable)

ftable_df <- data.frame(ftable)

colnames(ftable_df) <- c("ID","ModuleTitle","Category","Gene")

length(unique(ftable_df$ID))  # yes, 346 modules

unique(ftable_df$ID)

length(which(ftable_df$ID =="M1.0"))

write.table(ftable_df,file="genetable.txt",sep="\t",quote=F,col.names=T,row.names=F)
