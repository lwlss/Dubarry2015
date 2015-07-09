source("profilingFunctions.R")

# List of filenames to open
flist=list(
Control="FitnessReports/FitnessReport_MGD_CONT_CONC_30 SDM_rhlk_CTGNH.txt",
Cdc13="FitnessReports/FitnessReport_DJW_cdc13-1_DIL_28 SDM_rhlk_CTGNH.txt",
Pol_alpha="FitnessReports/FitnessReport_DAL_pol1-1_CONC_33 SDM_rhlk_CTGNH.txt",
Pol_epsilon="FitnessReports/FitnessReport_DAL_pol2-12_CONC_36 SDM_rhlk_CTGNH.txt",
Pol_delta="FitnessReports/FitnessReport_DAL_cdc2-2_CONC_30 SDM_rhlk_CTGNH.txt",
HU_30="FitnessReports/NormalisedHU.txt"
)

# Give columns pretty names
prettyNames=c("Control","cdc13-1","Pol \u03B1","Pol \u03B5","Pol \u03B4","HU 100mM","SD","ORF")

vals=list()
vals$prof=makeProfiles(flist,"MeanFit")
colnames(vals$prof)=prettyNames
vals$nExp=dim(vals$prof)[2]-2
vals$sims=similarities(vals$prof,"mahalanobis")
gdict=rownames(vals$prof)
names(gdict)=vals$prof$ORF

fc=read.delim("FunctionalComplexes.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)

# Discard rows where ORFs are not to be found in prof...
fcout=fc
for(i in 1:dim(fc)[1]){
	ORFs=strsplit(fc$GroupORFs[i]," ")[[1]]
	if(length(intersect(ORFs,names(gdict)))==0) reject=c(reject,i)
}
fc=fc[-1*reject,]
write.table(fc,"FunctionalComplexesTrimmed.txt",sep="\t",quote=FALSE,row.names=FALSE)


