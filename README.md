This package constitutes an interactive R problem set based on the RTutor package (https://github.com/skranz/RTutor). 

The problem set "Ecological Footprint of Poverty Alleviation" takes the user on a journey 
through the economic paper "The Ecological Footprint of Poverty Alleviation: 
Evidence from Mexico's Oportunidades Program" by J. Alix-Garcia et al. (2013). 
The economic findings along with the according analytic steps, as well as explanations of the economic theory
behind it and useful R commands in this context are covered in an interactive way.
The original paper can be found on http://www.mitpressjournals.org/doi/pdf/10.1162/REST_a_00349.

## 1. Installation

RTutor and this package is hosted on Github. To install everything, run the following code in your R console.
```s
if (!require(devtools))
  install.packages("devtools")
source_gist("gist.github.com/skranz/fad6062e5462c9d0efe4")
install.rtutor(update.github=TRUE)

devtools::install_github("kathkaufmann/RTutorEcologicalFootprintOfPovertyAlleviation", upgrade_dependencies=FALSE)
```

## 2. Show and work on the problem set
To start the problem set first create a working directory in which files like the data sets and your solution will be stored. Then adapt and run the following code.
```s
library(RTutorEcologicalFootprintOfPovertyAlleviation)

# Adapt your working directory to an existing folder
setwd("C:/problemsets/RTutorEcologicalFootprintOfPovertyAlleviation")
# Adapt your user name
run.ps(user.name="Jon Doe", package="RTutorEcologicalFootprintOfPovertyAlleviation",
       load.sav=TRUE, sample.solution=FALSE)
```
If everything works fine, a browser window should open, in which you can start exploring the problem set.
