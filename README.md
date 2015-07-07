# profilyzer_Dubarry2015

Profilyzer with data from Dubarry et al. 2015

![A static version of a profilyzer plot](Demo.png?raw=true)

This is a shiny app for interactively exploring [QFA](http://research.ncl.ac.uk) fitness profiles.  The plot above is a static version of the plots produced by profilyzer.

Currently there is only a private, as-yet unpublished dataset describing a set of screens targetted towards understanding DNA replication available, but profilyzer has been built to allow other datasets to be loaded in a straightforward manner.

To see profilyzer in action, a live instance of this particular dataset can be found at [this page](http://research.ncl.ac.uk/Dubarry2015).  Alternatively, you can download the code & data from this repository and run a local session.

To run profilyzer in a local session: 
* Open R session and set the current working directory to the one just above the profilyzer_Dubarry2015 directory
* Ensure that you have shiny installed.  If you do not, then execute the following in the R terminal: `install.packages("shiny")`
* Load the shiny library by executing the following in the R terminal: `library("shiny")`
* Finally, use profilyzer to browse the included LydallLab dataset: `runApp("profilyzer_Dubarry2015")`
