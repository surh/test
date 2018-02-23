setwd("~/micropopgen/src/mytestpackage/")

devtools::document()

pi.small <- read.table("~/micropopgen/src/Rgroup/pi.txt",
                       sep = "\t",
                       header = TRUE)
head(pi.small)

devtools::use_data(pi.small)
