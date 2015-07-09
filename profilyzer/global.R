fc=read.delim("FunctionalComplexesTrimmed.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
fc=rbind(c("None","None","None",""),fc)
choiceList=as.list(1:length(fc$Notes))
names(choiceList)=fc$GroupName