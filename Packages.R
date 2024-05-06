library(VSURF)
data("toys")
RFtoys <- VSURF(toys$x, toys$y, mtree = 100)
RFtoys$varselect.interp
RFtoys$varselect.pred
plot(RFtoys)
RFtoys$mean.perf
test <- VSURF(toys$x, toys$y, mtree = 100, data = c(1,2,3,4,5))

library(Boruta)
RFtoys2 <- Boruta(toys$x, toys$y)
RFtoys2
which(RFtoys2$finalDecision=='Confirmed')
plot(RFtoys2)

library(varSelRF)
RFtoys3 <- varSelRF(toys$x, toys$y)
RFtoys3

library(vita)
# Janitza
PerVarImp1 <- CVPVI(toys$x, toys$y)
RFtoys4 <- NTA(PerVarImp1$cv_varim)
summary(RFtoys4,pless = 0.01)
which(RFtoys4$pvalue<0.05)

# Altmann
PerVarImp2 <- PIMP(toys$x, toys$y, randomForest(toys$x, toys$y))
RFtoys5 <- PimpTest(PerVarImp2)
which(RFtoys5$pvalue<0.05)

library(CoVVSURF)
RFtoys6 <- covsurf(toys$x, toys$y)
RFtoys6$vsurf_ptree$varselect.interp
