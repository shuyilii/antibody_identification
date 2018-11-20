source("http://bioconductor.org/biocLite.R")

biocLite("msa", type = "source")

library(msa)
library("tools", lib.loc="/usr/lib/R/library")

setwd("~/Desktop/data/database")

mySequences <- readAAStringSet('alignment_0.fasta')

alignment <- msaMuscle(mySequences, gapOpening = -10000, gapExtension = -2000)

msaPrettyPrint(alignment, output="tex", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

texi2pdf("alignment.tex", clean=TRUE)

print(alignment)
