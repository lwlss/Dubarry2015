library(qfa)
sgd=readSGD()

commonStrip=c("YDR173C","YER069W","YHR018C","YJL071W","YJL088W","YML099C","YMR042W","YMR062C","YOL058W","YOL140W","YBR248C","YCL030C","YFR025C","YER055C",
              "YIL020C","YIL116W","YCL018W","YGL009C","YHR002W","YLR451W","YNL104C","YOR108W","YBR115C","YDL131W","YDL182W","YDR034C","YDR234W","YGL154C",
              "YIL094C","YIR034C","YNR050C","YMR038C") # leu arg his lys
 
exptStrip=c("YBR089W","YBR190W","YDL041W","YDR355C","YGL214W","YLL044W","YLR076C","YOR302W",
              "YGL105W","YDL004W","YGL130W","YOR303W","YJR109C","YDL160C","YHR059W","YLR088W",
              "YGR191W","YIL125W","YHR194W","YDR461W","YOL115W","YDL090C","YJL204C","YMR235C",
              "YOL005C","YGL123W","YBL025W","YLR357W","YDL042C","YLR442C","YDR227W","YDL052C",
              "YLR362W","YDR410C","YFL026W","YHL007C","YOR212W","YDR103W","YDL159W","YMR146C",
              "YIL011W","YGL026C","YMR038C") # specific to the SGA experiments

keepGenes=read.delim("Stripping.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
keepPol3=keepGenes$POL3
keepPol2=keepGenes$POL2
keepPol1=keepGenes$POL1
stripLyp1=keepGenes$NOTQFA18

editFitnessFile=function(fname,keepORFs=c(),stripORFs=c(),scaleMedian=NULL,scaleMean=NULL,fpre=""){
	frep=file(fname,"r")
	frepLines=readLines(frep)
	close(frep)
	toSkip=grep("Gene",frepLines)-1
	header=paste(frepLines[1:toSkip],collapse="\n")
	dat=read.delim(fname,skip=toSkip,header=TRUE,stringsAsFactors=FALSE)
	if(length(keepORFs)>0) dat=dat[dat$ORF%in%keepORFs,]
	if(length(stripORFs)>0) dat=dat[!dat$ORF%in%stripORFs,]
	if(!is.null(scaleMedian)) dat$MedianFit=dat$MedianFit/scaleMedian
	if(!is.null(scaleMean)) dat$MeanFit=dat$MeanFit/scaleMean
	print(header)
	if(fpre!="") fname=paste(fpre,fname,sep="_")
	write(header,file=fname)
	write.table(dat,file="tmpdat.out",sep="\t",quote=FALSE,row.names=FALSE)
	file.append(fname,"tmpdat.out")
	file.remove("tmpdat.out")
	return(dat)
}

cdc13=editFitnessFile("QFA0015_FitnessReport_DJW_cdc13-1_DIL_28_SDM_rhlk_CTGNH_MDRMDP.txt",stripORFs=c(commonStrip,exptStrip,getNeighbours(c("cdc13","lyp1","can1"),20,sgd)$FName))

HU=editFitnessFile("QFA0044_FitnessReport_EJA_G4Quad_Drugs_30_100mM_HU_CSM _MDRMDP.txt",stripORFs=c(commonStrip,exptStrip))
CSM=editFitnessFile("QFA0044_FitnessReport_EJA_G4Quad_Drugs_30_CSM Lydall_MDRMDP.txt",stripORFs=c(commonStrip,exptStrip))

lyp1_30=editFitnessFile("QFA0018_FitnessReport_MGD_lyp1_HLN_30_SDM_rhlk_CTGNH_MDRMDP.txt",stripORFs=c(commonStrip,exptStrip,getNeighbours(c("lyp1","can1"),20,sgd)$FName))
lyp1_33=editFitnessFile("QFA0018_FitnessReport_MGD_lyp1_HLN_33_SDM_rhlk_CTGNH_MDRMDP.txt",stripORFs=c(commonStrip,exptStrip,getNeighbours(c("lyp1","can1"),20,sgd)$FName))
lyp1_36=editFitnessFile("QFA0018_FitnessReport_MGD_lyp1_HLN_36_SDM_rhlk_CTGNH_MDRMDP.txt",stripORFs=c(commonStrip,exptStrip,getNeighbours(c("lyp1","can1"),20,sgd)$FName))

# Just need to normalise HU data so that the control fitness (QFA0044_FitnessReport_EJA_G4Quad_Drugs_30_CSM Lydall_MDRMDP.txt) matches that of QFA0018_FitnessReport_MGD_lyp1_HLN_30_SDM_rhlk_CTGNH_MDRMDP.txt
mean_lyp1=mean(lyp1_30$MeanFit)
median_lyp1=median(lyp1_30$MedianFit)

mean_CSM=mean(CSM$MeanFit)
median_CSM=median(CSM$MedianFit)

Huadj=editFitnessFile("QFA0044_FitnessReport_EJA_G4Quad_Drugs_30_100mM_HU_CSM _MDRMDP.txt",scaleMedian=median_CSM/median_lyp1,scaleMean=mean_CSM/mean_lyp1,fpre="SCALED")
cdc13adj=editFitnessFile("QFA0015_FitnessReport_DJW_cdc13-1_DIL_28_SDM_rhlk_CTGNH_MDRMDP.txt",scaleMedian=2,scaleMean=2,fpre="SCALED")