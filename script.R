library(readr)
name="r0_1_t0_1"
trace <- read_delim(paste(c("/home/fpinard/StatBay/BayesianPhylo/data/",name,".trace")), "\t", escape_double = FALSE, trim_ws = TRUE)
plot(trace$lnL,type="l")
acceptation <- read_delim(paste(c("/home/fpinard/StatBay/BayesianPhylo/data/",name,".acceptation")), "\t", escape_double = FALSE, trim_ws = TRUE)
plot(acceptation$acceptedRateMove/acceptation$NbCycle,type="l")
plot(acceptation$acceptedTimeMove/acceptation$NbCycle,type="l")
plot(acceptation$acceptedTopoMove/acceptation$NbCycle,type="l")
