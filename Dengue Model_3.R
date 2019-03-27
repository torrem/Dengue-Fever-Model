library(randomForest)
library(corrplot)
library(imputeTS)
library(colorRamps)


setwd("C:/Users/Mike/Documents/Dengue Fever")

###------ input datasets-------------------------------------####
Sub_Form = read.csv('submission_format.csv')
TrainFeat = read.csv('dengue_features_train.csv')
TrainLab = read.csv('dengue_labels_train.csv')
TestFeat = read.csv('dengue_features_test.csv')


All = cbind(TrainLab, TrainFeat); All = All[,c(1:4, 8:length(All))]



###------ Seperate data #####_

#for San Juan
sj_Train = subset(All, city == 'sj') 
sj_Test = subset(TestFeat, city == 'sj'); sj_Test$total_cases = rep('TEST',nrow(sj_Test))


# for Iquitos
iq_Train = subset(All, city == 'iq') 
iq_Test = subset(TestFeat, city == 'iq');iq_Test$total_cases = rep('TEST',nrow(iq_Test))



##------ fill in data gaps and pick city to run #####_
# sapply(sj_Train, function(x)all(any(is.na(x)))) # check which variables have missing values
# sapply(iq_Train, function(x)all(any(is.na(x)))) # check which variables have missing values
# 
# #plot missing values
# plot(sj_Train$ndvi_ne, type='l') 
# plotNA.distribution(sj_Train$ndvi_ne)
# plotNA.distributionBar(sj_Train$ndvi_ne)
# 
# dd = na.interpolation(sj_Train$ndvi_ne, option = "spline")
# plot(dd, type='l')
# plotNA.distribution(dd)

for (i in 6:length(sj_Train)){
  sj_Train[,i] = na.interpolation(sj_Train[,i], option = "linear")
}

for (i in 6:length(iq_Train)){
  iq_Train[,i] = na.interpolation(iq_Train[,i], option = "linear")
}

for (i in 5:length(sj_Test)){
  sj_Test[,i] = na.interpolation(sj_Test[,i], option = "linear")
}

for (i in 5:length(iq_Test)){
  iq_Test[,i] = na.interpolation(iq_Test[,i], option = "linear")
}
 
# sapply(sj_Train, function(x)all(any(is.na(x)))) # check which variables have missing values
# sapply(iq_Train, function(x)all(any(is.na(x)))) # check which variables have missing values



sj_Test = sj_Test[,c(1:3,25,4:24)]
sj_both = rbind(sj_Train, sj_Test)

iq_Test = iq_Test[,c(1:3,25,4:24)]
iq_both = rbind(iq_Train, iq_Test)





###------calculate better variables for training data -------------------------------#####
## from literature search it appears that temp, accumulated rainfall, and sunshine are important
## whats the difference between accumulated rainfall & standing water?
## it seems to me that standing water can be a function of accumulated rainfall, humudity, and temp. 
## additionally, seepage through sediment would eb important (impervious surface)? I don't have a value for this so maybe we can 
## permute over a range of values for this. So i need to come up with an equation to calculate standing water.


###----------Add optimized lag variable function  ####_



AdOpLagVar = function(x, mthd,  maxLag = 25, varName, varInputNum, plot = FALSE){
  
  y = subset(x, total_cases == 'TEST')
  x = subset(x, !(total_cases %in% 'TEST')); x$total_cases = as.numeric(x$total_cases)
  OpLag = data.frame(lag = 3:maxLag, score = rep(NA, maxLag-2))
  x$lag = rep(NA, nrow(x))
  y$lag = rep(NA, nrow(y))
  n = length(x)

  
  for (k in 1:nrow(OpLag)){
    timeLag = OpLag[k,1]
    colnames(x)[n] = paste(varName,"_lag_", timeLag, sep="")
    for (i in 1:nrow(x)){
      if (i < timeLag & mthd == 'mean' ){x[i,n] = mean(x[1:timeLag,varInputNum])}
      if (i < timeLag & mthd == 'sum' ){x[i,n] = mean(x[1:timeLag,varInputNum])*(timeLag-1)}
      if (i >= timeLag & mthd =='mean' ){x[i,n] = mean(x[(i-(timeLag-1)):i,varInputNum])}
      if (i >= timeLag & mthd =='sum' ){x[i,n] = sum(x[(i-(timeLag-1)):i,varInputNum])}
      
    }
    OpLag[k,2] =  cor(x = x[n], y = x$total_cases, method="spearman")
  }
  
  if (plot==TRUE){plot(OpLag$score ~ OpLag$lag)}
  
  colnames(x)[n] = 'lag'
  x = rbind(x,y)
  
  
  timeLag = which(OpLag$score == max(OpLag$score))+2
  colnames(x)[n] = paste(varName,"_lag_", timeLag, sep="")
  for (i in 1:nrow(x)){
    if (i < timeLag & mthd == 'mean' ){x[i,n] = mean(x[1:timeLag,varInputNum])}
    if (i < timeLag & mthd == 'sum' ){x[i,n] = mean(x[1:timeLag,varInputNum])*(timeLag-1)}
    if (i >= timeLag & mthd =='mean' ){x[i,n] = mean(x[(i-(timeLag-1)):i,varInputNum])}
    if (i >= timeLag & mthd =='sum' ){x[i,n] = sum(x[(i-(timeLag-1)):i,varInputNum])}
  }
  
  
  return(x)
  
}



###----------put variables in TRAIN  ####_

iq_both = AdOpLagVar(iq_both, mthd = 'sum', varName = 'accum_rain', varInputNum = 16, plot=TRUE) 
iq_both = AdOpLagVar(iq_both, mthd = 'mean', varName = 'accum_humid', varInputNum = 19, plot=TRUE) 
iq_both = AdOpLagVar(iq_both, mthd = 'mean', varName = 'accum_temp', varInputNum = 15, plot=TRUE) 



sj_both = AdOpLagVar(sj_both, mthd = 'sum', varName = 'accum_rain', varInputNum = 16, plot=TRUE) 
sj_both = AdOpLagVar(sj_both, mthd = 'mean', varName = 'accum_humid', varInputNum = 19, plot=TRUE) 
sj_both = AdOpLagVar(sj_both, mthd = 'mean', varName = 'accum_temp', varInputNum = 15, plot=TRUE) 

## separate test and train datasets
iq_Test = subset(iq_both, total_cases == 'TEST')
iq_Train = subset(iq_both, !(total_cases %in% 'TEST')); iq_Train$total_cases = as.numeric(iq_Train$total_cases)

sj_Test = subset(sj_both, total_cases == 'TEST')
sj_Train = subset(sj_both, !(total_cases %in% 'TEST'));sj_Train$total_cases = as.numeric(sj_Train$total_cases)


###------ look at correlation between variables -------------------------------#####

sj_m<-cor(sj_Train[,c(4,6:length(sj_Train))], method = 'spearman')

iq_m<-cor(iq_Train[,c(4,6:length(iq_Train))], method = 'spearman')


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat_sj <- cor.mtest(sj_Train[,c(4,6:length(sj_Train))])
p.mat_iq <- cor.mtest(iq_Train[,c(4,6:length(iq_Train))])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(sj_m, method="color", col=col(200),
         type="upper", order="original",
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.5 , tl.col="black", tl.srt=15, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_sj, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=TRUE )

corrplot(iq_m, method="color", col=col(200),
         type="upper", order="original",
         addCoef.col = "black", # Add coefficient of correlation
         tl.cex = 0.5 , tl.col="black", tl.srt=15, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_sj, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=TRUE )








###------ Build Random Forest Regression Model ---------------------------   #####


## SJ model ##
rf_sj =randomForest(x = sj_Train[,c(13,15,16,19,26,27,28 )], y= sj_Train$total_cases,
                 mtry=2, ntree = 1000, data = sj_Train)  # all areas modeled together
rf_sj

sj_Train = cbind(sj_Train,predCases = predict(rf_sj, newdata=sj_Train, type='response'))
mean(abs(sj_Train$total_cases - sj_Train$predCases))

sj_Test = cbind(sj_Test,predCases = predict(rf_sj, newdata=sj_Test, type='response'))

## IQ model ##
rf_iq =randomForest(x = iq_Train[,c(13,15, 19, 24,26,27,28)], y= iq_Train$total_cases,
                      mtry=2, ntree = 1000, data = iq_Train)  # all areas modeled together
rf_iq

iq_Train = cbind(iq_Train,predCases = predict(rf_iq, newdata=iq_Train, type='response'))
mean(abs(iq_Train$total_cases - iq_Train$predCases))
iq_Test = cbind(iq_Test,predCases = predict(rf_iq, newdata=iq_Test, type='response'))


###------------------PLOT---------------------------------------####

par(mfrow=c(1,1))

## plot SJ ##
sj_Train$x = 1:nrow(sj_Train) 
sj_Test$x = (nrow(sj_Train) + 1):((nrow(sj_Train) + 1)+(nrow(sj_Test)-1))

plot(sj_Train$total_cases~sj_Train$x, xlim = c(0,(nrow(sj_Train) + nrow(sj_Test))),
     col = 'blue', type='l', 
     main = paste('San Juan', "        MSAE = ",mean(abs(sj_Train$total_cases - sj_Train$predCases)) ))
lines(sj_Train$predCases~sj_Train$x, col = 'red')
lines(sj_Test$predCases~sj_Test$x, col = 'green')
legend(x = 0.6*(nrow(sj_Train) + nrow(sj_Test)), y = 0.9*max(sj_Train$total_cases),
       legend=c("Total Cases", "Predicted Cases", "Test Cases"),
       col=c("blue", "red", "Green"), lty=1:1, cex=0.7)


## plot IQ ##
iq_Train$x = 1:nrow(iq_Train) 
iq_Test$x = (nrow(iq_Train) + 1):((nrow(iq_Train) + 1)+(nrow(iq_Test)-1))

plot(iq_Train$total_cases~iq_Train$x, xlim = c(0,(nrow(iq_Train) + nrow(iq_Test))),
     col = 'blue', type='l', 
     main = paste('Iquitos', "        MSAE = ",mean(abs(iq_Train$total_cases - iq_Train$predCases)) ))
lines(iq_Train$predCases~iq_Train$x, col = 'red')
lines(iq_Test$predCases~iq_Test$x, col = 'green')
legend(x = 0.6*(nrow(iq_Train) + nrow(iq_Test)), y = 0.9*max(iq_Train$total_cases),
       legend=c("Total Cases", "Predicted Cases", "Test Cases"),
       col=c("blue", "red", "Green"), lty=1:1, cex=0.8)






#-----------------------compile form to submit-------------------

Sub = rbind(sj_Test[,c(1:3,29)], iq_Test[,c(1:3,29)]); colnames(Sub)[4] = 'total_cases'

Sub$total_cases = round(Sub$total_cases)


write.csv(Sub, row.names = FALSE, 'sub3.csv') # may need to remove rows

#Test rounding














#-------------------------------END---------------------------------------------
# City and date indicators
# 
# city - City abbreviations: sj for San Juan and iq for Iquitos
# week_start_date - Date given in yyyy-mm-dd format
# 
# NOAA's GHCN daily climate data weather station measurements
# 
#     station_max_temp_c - Maximum temperature
#     station_min_temp_c - Minimum temperature
#     station_avg_temp_c - Average temperature
#     station_precip_mm - Total precipitation
#     station_diur_temp_rng_c - Diurnal temperature range
# 
# PERSIANN satellite precipitation measurements (0.25x0.25 degree scale)
# 
#     precipitation_amt_mm - Total precipitation
# 
# NOAA's NCEP Climate Forecast System Reanalysis measurements (0.5x0.5 degree scale)
# 
# reanalysis_sat_precip_amt_mm - Total precipitation
# reanalysis_dew_point_temp_k - Mean dew point temperature
# reanalysis_air_temp_k - Mean air temperature
# reanalysis_relative_humidity_percent - Mean relative humidity
# reanalysis_specific_humidity_g_per_kg - Mean specific humidity
# reanalysis_precip_amt_kg_per_m2 - Total precipitation
# reanalysis_max_air_temp_k - Maximum air temperature
# reanalysis_min_air_temp_k - Minimum air temperature
# reanalysis_avg_temp_k - Average air temperature
# reanalysis_tdtr_k - Diurnal temperature range
# 
# Satellite vegetation - Normalized difference vegetation index (NDVI) - NOAA's CDR Normalized Difference Vegetation Index (0.5x0.5 degree scale) measurements
# 
#     ndvi_se - Pixel southeast of city centroid
#     ndvi_sw - Pixel southwest of city centroid
#     ndvi_ne - Pixel northeast of city centroid
#     ndvi_nw - Pixel northwest of city centroid

## cool color ramp
## lette(c("black", "#202020", "#736AFF", "cyan", "yellow", "#F87431", "#FF7F00", "red", "#7E2217"))



## Use RF model to iteratively fill out recent # of cases with a time lag ##

#####need to build RF model and add on total cases iteratively

# ## first step is to combine Train and Test Data
# 
# sj_Train$TT = rep('Train', nrow(sj_Train))
# sj_Test$TT = rep('Test', nrow(sj_Test))
# 
# sj_Test$total_cases = rep('NA', nrow(sj_Test))
# sj_Test$predCases = rep('NA', nrow(sj_Test))
# sj_Test$predCases = rep('NA', nrow(sj_Test))
# sj_Test$lag = rep(NA, nrow(sj_Test));colnames(sj_Test)[30] = paste('Acum_cases_lag', timeLagCase, sep="_")
# 
# sj_Train = sj_Train[,c(1,30,2:29)]
# 
# sj_Test = sj_Test[,c(1,27,2,3,28, 4:26, 30, 29 )]
# 
# 
# 
# 
# sj_All = rbind(sj_Train, sj_Test)
# 
# 
# for (i in 937:nrow(sj_All)){
#   
#   sj_All[i,29] = sum(as.numeric(sj_All$total_cases[(i-(timeLagCase+20)):i-1]))
#      sj_All[i,5]          = predict(rf_sj, newdata=sj_All[i,], type='response')
#      sj_All[i,30] = sj_All[i,5]
# }
# 
# 
# plot(sj_All$total_cases, col = 'blue', type='l')
# lines(sj_Train$predCases, col = 'red')
# legend(x = 675,y = 400,legend=c("Total Cases", "Predicted Cases"),
#        col=c("blue", "red"), lty=1:1, cex=0.8)



## boosted regression tree###? which variables are important
## multicollinearity between variables


## The benchmark was a negative binomial relationship



