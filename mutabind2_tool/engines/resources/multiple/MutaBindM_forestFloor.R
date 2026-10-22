library(forestFloor)
library(randomForest)
load("minputfiles/mutabindm.RData")
test <- read.csv(file='test_outfeature',header = TRUE, sep = '\t')
test_info <- read.csv(file='test_sunddg',header = TRUE, sep = '\t') 
test_features <- test[,9:(ncol(test)-1)]
test_features$ACC_wt <- as.integer(test_features$ACC_wt)
test_features$ACCmp_wt <- as.integer(test_features$ACCmp_wt)

ff <- forestFloor(mutabindm.rf,data[-1],test_features,bootstrapFC = TRUE, calc_np = TRUE)
fc <- round(ff$FCmatrix,4)
fc <- as.data.frame(fc)
colnames(fc) <- c('DDE_vdw','DDG_solv','DDG_fold','SA_com_wt','SA_part_wt','CS','dE_vdw_wt')
fc <- fc[,c('DDE_vdw','DDG_solv','DDG_fold','CS','SA_com_wt','SA_part_wt','dE_vdw_wt')]
result <- cbind(test_info,fc[(nrow(fc)-nrow(test)+1):nrow(fc),])
write.table(result,file = 'test_outcontribution',sep = '\t',col.names = TRUE,row.names = FALSE,quote=FALSE)
