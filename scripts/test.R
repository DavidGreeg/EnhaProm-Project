
#Look for installed packages
installed.packages()

#Load library we gonna use
library(tidyverse)

#Set working directory
setwd("/home/davidfm/Epigenomica/UBMI-IFC")

#Load data
data <- read.csv("RefSeq_FuncElems_chr1.csv")

#View first data rows
head(data)
