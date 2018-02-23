pi.small <- read.table("~/micropopgen/src/Rgroup/pi.txt", sep = "\t",
                       header = TRUE)
head(pi.small)

m1 <- linreg(Pi_content ~ Bacteria + EndP + Experiment, pi.small)
m1
summary(m1)

m2 <- lm(Pi_content ~ Bacteria + EndP + Experiment, pi.small)
summary(m2)


m1.perm <- permtest(m1, nperm = 999)
