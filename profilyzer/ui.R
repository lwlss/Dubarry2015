library(shiny)

shinyUI(fluidPage(title="Profilyzer",
	titlePanel(h1("Profilyzer: Dubarry et al. (2015)")),
	HTML('<p>To explore evidence for genetic interaction in detail, use our fitness plot visualisation tool: <a href="http://bsu-srv.ncl.ac.uk/dixy-pol/">DIXY</a></p>'),
  
  fluidRow(
      
    column(2, 
      checkboxGroupInput("checkGroup", 
        label = h4("Screens to profile"), 
        choices = list("Control" = 1, "cdc13-1" = 2, "Pol alpha" = 3, "Pol delta" = 4, "Pol epsilon" = 5,  "HU 100mM" = 6),
        selected = c(1,2,3,4,5,6)),
		HTML('<h6>All screens presented here, except HU 100mM, include the <a href="http://dharmacon.gelifesciences.com/non-mammalian-cdna-and-orf/yeast-damp-collection/">DAmP</a> allele collection.</h6>')
		),
		
	column(2,
		radioButtons("pType",
			label=h4("Fitness summary"),
            choices=c("Mean" = "MeanFit","Median" = "MedianFit"),
			selected="MeanFit"
			   ),
		HTML('<h6>Method of summarising (<a href="https://en.wikipedia.org/wiki/Average">averaging</a>) replicate fitness observations for strains with individual genotypes.</h6>')
		),
	column(2,
		radioButtons("dType",
			label=h4("Distance Measure"),
            choices=c("Euclidean" = "distance","Correlation" = "correlation","Mahalanobis" = "mahalanobis"),
			selected="distance"
			   ),
		HTML('<h6><a href="https://en.wikipedia.org/wiki/Euclidean_distance">Euclidean distance</a> discovers profiles that are close to each other or overlap.  <a href="https://en.wikipedia.org/wiki/Correlation_and_dependence">Correlation</a> discovers profiles with the same pattern, but potentially offset from each other.  <a href="https://en.wikipedia.org/wiki/Mahalanobis_distance">Mahalanobis distance</a> can be considered as an intermediate between the two.</h6>')
			   
			   ),
		
	column(3,
		textInput("glist", label = h4("Target gene(s)"),value = "rad24"),
		HTML('<h6>Nearest or farthest profiles are only highlighted when just a single target gene name is specified in the box above.  Note that nearest or farthest profiles can only be identified when the fitness of strains carrying mutations in the target gene is measured in all selected screens.</h6>'),
		selectInput("ggroup", label = h4("Gene group"),choices = choiceList,selected = 1,selectize=TRUE),
		h6("Optionally select groups of functionally related target genes from the drop-down list above instead of specifying them manually.  Select 'None' to return to highlighting 'Target gene(s)'")
		), 
    
    column(2,
		numericInput("nearNum",label = h4("No. profiles nearest to target"),value = 12,min=0),
		numericInput("farNum",label = h4("No. profiles farthest from target"),value = 0,min=0),
		h6("Select number of nearest (or farthest) profiles to plot & tabulate alongside a single 'Target gene'.  Nearest genes may have a similar function to the target.  Farthest genes may have an opposite function, but such observations are probably more difficult to interpret.")
		) 
    ),
	  
   fluidRow(
	plotOutput("profiles", height="1000px", width = "100%")
	),
	downloadButton('downloadPlot', 'Download Plot'),
	fluidRow(
	column(2,
	h4("Nearest"),
	htmlOutput("nearest")
	),
	column(2,
	h4("Farthest"),
	htmlOutput("farthest")
	),
	column(8,
	plotOutput("ranks", height="500px", width = "100%")
	)	
   ),
   HTML('<p>Data, source code & documentation for this instance of profilyzer are hosted on <a href="https://github.com/lwlss/profilyzer">GitHub</a></p>')
))
