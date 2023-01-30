#installing deps

# CountsQC

if(!require("optparse")) { 
    install.packages("optparse", repos = "http://cran.r-project.org") 
}

if(!require("NOISeq")) {
    BiocManager::install("NOISeq")
}

# Epigenetics

if(!require("XML")) { 
    install.packages("XML", repos = "http://cran.r-project.org")
}
   
if(!require("Repitools")) {
    BiocManager::install("Repitools")
}

if(!require("Rsamtools")) {
    BiocManager::install("Rsamtools")
}

if(!require("rtracklayer")) {
    BiocManager::install("rtracklayer")
}

