library(shiny)
library(xtable)
source("profilingFunctions.R",local=TRUE)

# List of filenames to open
flist=list(
Control="FitnessReports/QFA0018_FitnessReport_30_MDRMDP.txt",
Cdc13="FitnessReports/QFA0015_FitnessReport_28_MDRMDP.txt",
Pol_alpha="FitnessReports/QFA0010_FitnessReport_33_MDRMDP.txt",
Pol_delta="FitnessReports/QFA0012_FitnessReport_30_MDRMDP.txt",
Pol_epsilon="FitnessReports/QFA0011_FitnessReport_36_MDRMDP.txt",
HU_30="FitnessReports/FIT_G4QUADDRUGS_100mM HU 30_SCALED.txt"
)

# Give columns pretty names
prettyNames=c("Control","cdc13-1","Pol \u03B1","Pol \u03B4","Pol \u03B5","HU 100mM","SD","ORF")

# Load fitness data consistent with default user input values
vals=list()
vals$prof=makeProfiles(flist,"MeanFit")
colnames(vals$prof)=prettyNames
vals$nExp=dim(vals$prof)[2]-2
vals$sims=similarities(vals$prof,"mahalanobis")
gdict=rownames(vals$prof)
names(gdict)=vals$prof$ORF

shinyServer(function(input, output) {

		updateScreens <- reactive({
			vals = list()
			vals$prof=makeProfiles(flist[as.numeric(input$checkGroup)],input$pType)
			vals$nExp=dim(vals$prof)[2]-2
			vals$sims=similarities(vals$prof,input$dType)
			colnames(vals$prof)=prettyNames[as.numeric(input$checkGroup)]
			vals$simtarg=c()
			return(vals)
		})
		
		updateTarget <- reactive({
			vals = updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			if((length(gene)==1)&(gene%in%rownames(vals$sims))) {
				simtarg = findSimilar(gene,vals$sims)
			}else{
				simtarg=c()
			}
			return(simtarg)
		})
		
		output$downloadPlot <- downloadHandler(
			filename = function() { "FitnessProfile.pdf" },
			content = function(file) {
			vals = updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			simtarg = updateTarget()
			near = head(simtarg,n=input$nearNum)
			far = tail(simtarg,n=input$farNum)
			cairo_pdf(file,width=14,height=10,onefile=TRUE)
				print(plotSimilarFit(gene,vals$prof,nearest=near,farthest=far))
			dev.off()
		})
		
		output$profiles <- renderPlot({
			vals = updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				print(fc$GroupORFs[input$ggroup])
				print(input$ggroup)
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			simtarg = updateTarget()
			near = head(simtarg,n=input$nearNum)
			far = tail(simtarg,n=input$farNum)
			plotSimilarFit(gene,vals$prof,nearest=near,farthest=far) 
		})
		
		output$ranks <- renderPlot({
			vals = updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			gene=intersect(gene,rownames(vals$sims))
			plotGenomewideSimilarity(gene,vals$sims)
		})
			
		output$nearest <- reactiveText(function(){
			vals = updateScreens()
			simtarg = updateTarget()
			if((length(simtarg)>0)&(input$nearNum>0)){
				near = head(simtarg,n=input$nearNum)
				dftab=data.frame('Gene'=paste('<a href="http://www.yeastgenome.org/locus/',tolower(names(near)),'/overview" target="_blank">',tolower(names(near)),'</a>',sep=""),Distance=near,row.names=NULL)
				colnames(dftab)=c("Gene Name","Distance")
				return(as.character(print(xtable(dftab,align=c("r","c","c")),type="html",sanitize.text.function = function(x){x})))
			}else{return("")}
		})
			
		output$farthest <- reactiveText(function(){
			vals = updateScreens()
			simtarg = updateTarget()
			if((length(simtarg)>0)&(input$farNum>0)){
				far = rev(tail(simtarg,n=input$farNum))
				dftab=data.frame(Gene=paste('<a href="http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=',tolower(names(far)),'" target="_blank">',tolower(names(far)),'</a>',sep=""),Distance=far,row.names=NULL)
				colnames(dftab)=c("Gene Name","Distance")
				return(as.character(print(xtable(dftab,align=c("r","c","c")),type="html",sanitize.text.function = function(x){x})))
			}else{return("")}
		})
})
