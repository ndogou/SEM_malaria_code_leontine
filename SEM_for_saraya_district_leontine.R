par(mfrow=c(1,1))
rm(list=ls())
library(GPArotation)
library(psych)
library(ggplot2)
library(lavaan)
library(lavaanPlot)
library(semPlot)
library(semTools)
library(caret)
setwd("C:/Users/USER/Desktop/refrence_travaux/refrence_equation_structurelle")
data_saraya<-read.csv("data_saraya_num.csv",header = TRUE,sep = ",")
data_sarayaA<-read.csv("data_saraya.csv",header = TRUE,sep = ",")


head(data_saraya)
str(data_saraya)
head(data_saraya)
names(data_saraya)
attach(data_saraya)
summary(is.na(data_saraya))
#### Convert variable integer to numerical#############
# Identify integer colones

# Identify integer colones
cols_integer <- sapply(data_saraya, is.integer)
cols_integer
# Convert this variable to numeric
data_saraya[ , cols_integer] <- lapply(data_saraya[ , cols_integer], as.numeric)

str(data_saraya)
data_saraya<-as.data.frame(data_saraya)

# Identifier les variables à faible variance
var_faible_variance_sar <- nearZeroVar(data_saraya, saveMetrics = TRUE)
print(var_faible_variance_sar)

# Afficher les variables identifiées
print(var_faible_variance_sar[var_faible_variance_sar$nzv, ])
#Remove thise variables
data_saraya <- data_saraya[, !var_faible_variance_sar$nzv]
names(data_saraya)
dim(data_saraya)

#selectionner des colones
summary(data_saraya)

data_saraya1 <- subset(data_saraya, select = -c(precipprob_min, precipprob_max,precipcover_min,feelslikemin_min,feelslikemax_min,feelslikemin_max,feelslikemax_max,feelslikemin_mean,tempmin_max,tempmax_max,feelslikemax_mean,tempmin_mean,tempmax_mean,tempmax_min,tempmin_min,
                                                sealevelpressure_mean,winddir_mean,winddir_max,winddir_min,mean_night_duration,temp_mean, humidity_mean,precipcover_mean,windspeed_mean,dew_mean,precipprob_mean,feelslike_mean,windgust_mean,precip_mean,cloudcover_mean,solarradiation_mean,solarenergy_mean,uvindex_mean))
names(data_saraya1)
dim(data_saraya1)
write.csv(data_saraya1 , "C:/Users/USER/Desktop/refrence_travaux/refrence_equation_structurelle/data_saraya1.csv", row.names = FALSE)

#
library(psych)
library(Hmisc)
KMO(cor(data_no_outliers1_reduced_sar))
#data_saraya1<- subset(data_saraya1, select = -c(moonphase_min,moonphase_mean,moonphase_max,VHI_mean))
#moonphase_min,moonphase_mean,moonphase_max sont supprimées 
# car il on un KMO inferieur à 0.5

# bartelette test
#  package "Hmisc"

# Correlation 
res_rcorr_sar=rcorr(as.matrix(data_no_outliers1_reduced_sar), type="pearson")
matrice_cor_sar = res_rcorr_sar$r
matrice_cor_sar
#test de sphéricité de Bartlett
cortest.bartlett(matrice_cor_sar,n=100)


####Multivariate normality check####
data_saraya1<-as.data.frame(data_saraya1)
ncol(data_saraya1) 
dim(data_saraya1)
library(MVN)###### For multivariate normality check##
Normality_sar<-mvn(data_no_outliers1_reduced_sar)
Normality_sar
sapply(data_saraya1, is.numeric)
str(data_saraya1)
data_saraya1 <- na.omit(data_saraya1)
# correlation  analysis
###For item level correlation matrix###
Correlation_sar<-round(cor(data_saraya1),2)
Correlation_sar[upper.tri(Correlation_sar)]<-''
Correlation_sar<-as.data.frame(Correlation_sar)
Correlation_sar

cor_matrix <- cor(data_saraya1, use = "pairwise.complete.obs")
cor_matrix

# 
lower_matrix <- lower.tri(cor_matrix, diag = TRUE)
cor_lower <- cor_matrix 
A<-cor_lower[!lower_matrix] <- NA
A




corrplot(cor_matrix, method = "color", type = "lower", 
         tl.col = "black", tl.srt =45,tl.cex=.5)

library(corrplot)
corrplot(cor(data_saraya1, use="complete.obs"), order = "hclust", tl.col='black', tl.cex=.41)

#avec  corrplot
corrplot(cor(data_saraya1, use="complete.obs"),order="hclust", tl.col='black', tl.cex=.41, addrect=6) # avec un clustering des variables

library(RcmdrMisc)

rcorr.adjust(data_saraya1) # This function is build into R Commander.


#data_saraya1 <- scale(data_saraya1)|> as.data.frame()

linear_model1 <- lm(cas_palu_confirme ~ ., data =data_saraya1)
summary(linear_model1)
library(dplyr)

# outliers with Mahalanobis distance

#  Mahalanobis distance
distances_sar <- mahalanobis(data_saraya1, colMeans(data_saraya1), cov(data_saraya1))

# Calculation of the threshold
p_sar <- ncol(data_saraya1)
seuil_sar <- qchisq(0.95, p_sar) 
#0.975

# outliers
outliers_sar <- which(distances_sar> seuil_sar)

# 
print(paste("Distances de Mahalanobis:", distances_sar))
print(paste("Seuil:", seuil_sar))
print(paste("Indices des outliers:", outliers_sar))

data_no_outliers_sar<- data_saraya1[-outliers_sar, ]
dim(data_no_outliers_sar)


linear_model1_sar <- lm(cas_palu_confirme ~ ., data =data_no_outliers_sar)
summary(linear_model1_sar)
library(dplyr)

vif(linear_model1_sar)
# Identify variable with VIF>10
vif_values_sar <- vif(linear_model1_sar)

# Keep only variables with VIF ≤ 10
variables_to_keep_sar <- names(vif_values_sar[vif_values_sar<= 10])
variables_to_keep_sar



data_no_outliers1_reduced_sar <- data_no_outliers_sar[, c("cas_palu_confirme", variables_to_keep_sar)]
dim(data_no_outliers1_reduced_sar)
names(data_no_outliers1_reduced_sar)
summary(data_no_outliers1_reduced_sar)
cor_matrix <- cor(data_no_outliers1_reduced_sar, use = "pairwise.complete.obs")
cor_matrix
#
cor_matrix_2dec <- round(cor_matrix, 2)
print(cor_matrix_2dec)

# 
cor_lower <- cor_matrix
cor_lower[upper.tri(cor_lower)] <- NA
cor_lower_2dec <- round(cor_lower, 2)
print(cor_lower_2dec)

# 
summary(data_no_outliers1_reduced_sar)
# 
variances <- sapply(data_no_outliers1_reduced_sar, var, na.rm = TRUE)
print(variances)
sd <- sapply(data_no_outliers1_reduced_sar, sd, na.rm = TRUE)
print(sd)
sd <- round(sd, 2)







#  Data splitting for EFA and CFA
set.seed(123)
n_total_sar <- nrow(data_no_outliers1_reduced_sar)
train_indices_sar <- sample(1:n_total_sar, size = floor(0.8 * n_total_sar))


#  AFE data 
data_sar_afe <- data_no_outliers1_reduced_sar[train_indices_sar, ]
subsetMal_sar_afe <- data_sar_afe[, apply(data_sar_afe, 2, function(x) var(x, na.rm = TRUE) > 0.008 & !all(is.na(x)))]
subsetMal_sar_afe <- scale(subsetMal_sar_afe) |> as.data.frame()

#  CFA data  
data_sar_cfa <- data_no_outliers1_reduced_sar[-train_indices_sar, ]
subsetMal_sar_cfa <- data_sar_cfa[, colnames(subsetMal_sar_afe)]
subsetMal_sar_cfa <- scale(subsetMal_sar_cfa) |> as.data.frame()
subsetMal_sar_cfa<- subsetMal_sar_cfa[, apply(subsetMal_sar_cfa, 2, function(x) var(x, na.rm = TRUE) > 0.008 & !all(is.na(x)))]







# Factorial exploratory analysis
set.seed(123)
data_no_outliers1_reduced_sar <- scale(data_no_outliers1_reduced_sar)|> as.data.frame()
frr_data_saraya2_sar<-char2numeric(data_no_outliers1_reduced_sar)

parallel_saraya <-fa.parallel(subsetMal_sar_afe,fm="ml",fa='both',n.iter=1000,
                              main="Parallel Analysis Scree Plots",SMC=T)

parallel_saraya <-fa.parallel(frr_data_saraya2_sar,fm="MLR",fa='both',n.iter=3000,
                              main="Parallel Analysis Scree Plots",SMC=T)

library(car)
parallel_saraya$nfact #number of factors
parallel_saraya$fa.values
#####Help to determine the exact number of factor###

scree(subsetMal_sar_afe,factors=TRUE,pc=TRUE,main='Scree plot',hline=NULL,add=FALSE)
ff_sar <- cor(subsetMal_sar_afe, use = "pairwise.complete")
VSS.scree(ff_sar)

# 
par(bg = "whitesmoke",                
    cex = 1,                   
    cex.main =1.5,               
    cex.lab = 1.3,              
    cex.axis = 1.5,               
    lwd = 3,                      
    font.main =2,                
    font.lab =1) 

#  scree plot
scree(subsetMal_sar_afe,
      factors = TRUE,
      pc = TRUE,
      main = 'saraya Makha',
      hline = NULL,
      add = FALSE)               


# 

ev_saraya <- eigen(cor(subsetMal_sar_afe)) # get eigenvalues
ev_saraya$values

######simple ###

library(plspm)
library(seminr)
library(Rcsdp)

set.seed(123)
subsetMal_numeric_saraya <- data_no_outliers1_reduced_sar[, apply(data_no_outliers1_reduced_sar, 2, function(x) var(x, na.rm = TRUE) > 0.008 & !all(is.na(x)))]
subsetMal_numeric_saraya <- scale(subsetMal_numeric_saraya)|> as.data.frame()

Twofactor_saraya_ml<- fa(subsetMal_sar_afe,nfactors =4,rotate ="oblimin",fm="ml",SMC =T)
Twofactor_saraya_ml 


#===============================================================

Twofactor_saraya_ml1<- fa(subsetMal_sar_afe,nfactors =5,rotate ="oblimin",fm="ml",SMC =T)
Twofactor_saraya_ml1 



# Modele  with 4 factors
chi2_4sar <- Twofactor_saraya_ml$chi
df_4sar <- Twofactor_saraya_ml$dof
p_4sar <- Twofactor_saraya_ml$PVAL
rmsea_4sar <- Twofactor_saraya_ml$RMSEA[1]
rmsea_ci_4sar <- Twofactor_saraya_ml$RMSEA[c(2,3)]  # Interval of confidence
tli_4sar <- Twofactor_saraya_ml$TLI
cfi_4sar <- Twofactor_saraya_ml$CFI
bic_4sar <- Twofactor_saraya_ml$BIC
aic_4sar <- Twofactor_saraya_ml$aic

# Model with 5 facteur  
chi2_5sar <- Twofactor_saraya_ml1$chi
df_5sar <- Twofactor_saraya_ml1$dof
p_5sar<- Twofactor_saraya_ml1$PVAL
rmsea_5sar <- Twofactor_saraya_ml1$RMSEA[1]
rmsea_ci_5sar <- Twofactor_saraya_ml$RMSEA[c(2,3)]  # Intervalle de confiance
tli_5sar <- Twofactor_saraya_ml1$TLI
cfi_5sar <- Twofactor_saraya_ml1$CFI
bic_5sar <- Twofactor_saraya_ml1$BIC
aic_5sar <- Twofactor_saraya_ml1$aic

# ===============================================================
#  MODELS Comparison
# ===============================================================


cat(sprintf("Chi-carré = %.3f (df = %d, p = %.10f)\n", chi2_4sar, df_4sar, p_4sar))
cat(sprintf("RMSEA = %.4f [%.4f; %.4f]\n", rmsea_4sar, rmsea_ci_4sar[1], rmsea_ci_4sar[2]))
cat(sprintf("TLI = %.4f\n", tli_4sar))
cat(sprintf("CFI = %.4f\n", cfi_4sar))
cat(sprintf("BIC = %.2f\n", bic_4sar))
cat(sprintf("AIC = %.2f\n", aic_4sar))

cat("\n MODeLE  with 5 FACTORS ---\n")
cat(sprintf("Chi-carré = %.3f (df = %d, p = %.10f)\n", chi2_5sar, df_5sar, p_5sar))
cat(sprintf("RMSEA = %.4f [%.4f; %.4f]\n", rmsea_5sar, rmsea_ci_5sar[1], rmsea_ci_5sar[2]))
cat(sprintf("TLI = %.4f\n", tli_5sar))
cat(sprintf("CFI = %.4f\n", cfi_5sar))
cat(sprintf("BIC = %.2f\n", bic_5sar))
cat(sprintf("AIC = %.2f\n", aic_5sar))


# ===============================================================
# 3.(LIKELIHOOD RATIO TEST)
# ===============================================================



#
chi2_diffsar <- chi2_4sar - chi2_5sar
df_diffsar <- df_4sar - df_5sar
p_diffsar <- 1 - pchisq(chi2_diffsar, df_diffsar)

cat("H0: The 4-factor model fits as well as the 5-factor model\n")
cat("H1: The 5-factor model fits significantly better\n\n")

cat(sprintf("Δχ² = %.3f\n", chi2_diffsar))
cat(sprintf("Δdf = %d\n", df_diffsar)) 
cat(sprintf("p-value = %.20f\n", p_diffsar))

#Interpretation of the test
if (p_diffsar < 0.001) {
  decision_lrsar<- "H0 is strongly rejected (p < 0.001)"
  conclusion_lrsar <- "The 5-factor model is much better."
} else if (p_diffsar < 0.01) {
  decision_lrsar<- "H0 is rejected at the 1% significance level"
  conclusion_lrsar <- "The 5-factor model provides a significantly better fit to the data"
} else if (p_diffsar < 0.05) {
  decision_lrsar <- "Rejection of H0 (p < 0.05)"
  conclusion_lrsar <- "The 5-factor model is better."
} else {
  decision_lrsar <- "NO Rejection of H0 (p ≥ 0.05)"
  conclusion_lrsar <- "No significant difference; the 4-factor model is preferred for parsimony."
}

cat(sprintf("\nDecision: %s\n", decision_lrsar))
cat(sprintf("Conclusion: %s\n", conclusion_lrsar))




















#===============================================


print(Twofactor_saraya_ml$communality)
print(Twofactor_saraya_ml$loadings)
communalities <- Twofactor_saraya_ml$communality
communalities


# Identifier les variables avec une communalité inférieure à 0,5
low_communalities <- which(communalities < 0.5)

names(low_communalities)

#===================================================================


summary(Twofactor_saraya_ml)
print(Twofactor_saraya_ml,cut = 0.3,digits =8)
fa.diagram(Twofactor_saraya_ml,sort=TRUE, cut=.3, simple=TRUE, errors=FALSE, digits=2, e.size=.01, rsize=0.11)
plot(Twofactor_saraya_ml)
g<-fa.plot(Twofactor_saraya_ml, labels=TRUE)
# Print factor loadings
print(Twofactor_saraya_ml)


png("C:/Users/USER/Desktop/refrence_equation_structurelle/resultat/saraya/fa_diagram.png", width = 2000, height = 1000, res = 100)
fa.diagram(Twofactor_saraya_ml,sort=TRUE, cut=.1, simple=TRUE, errors=FALSE, digits=2, e.size=.03, rsize=0.5,cex =5)
dev.off()





# essaye 1
cfaMod1_saraya1<-'f1=~windgust_max+windspeed_max
                   f2=~moonphase_max
                   f3=~ precip_max+precipcover_max+cloudcover_min
                   f4=~windgust_min+windspeed_min 
                   f5=~sealevelpressure_min+sealevelpressure_max
                   
                   
              '



data_saraya1_scaled <- scale(data_no_outliers1_reduced_sar)|> as.data.frame()
fit.cfaMod_saraya1<-cfa(cfaMod1_saraya1,data=subsetMal_sar_cfa,estimator = "ml",std.lv=T)
summary(fit.cfaMod_saraya1,fit.measures=TRUE,standardized=TRUE,rsq=TRUE)
summary(fit.cfaMod_saraya1, standardized = TRUE, fit.measures = TRUE)
fitMeasures(fit.cfaMod_saraya1, fit.measures = c("rmsea", "cfi", "tli", "srmr", "aic", "bic"))
cat(cfaMod1_saraya1)
modindices(fit.cfaMod_saraya1,sort = TRUE, maximum.number = 10)


#===================================================================
inspect(fit.cfaMod_saraya1, what = "std")
inspect(fit.cfaMod_saraya1, "converged")  
inspect(fit.cfaMod_saraya1, "theta") 

varTable(fit.cfaMod_saraya1)
lavResiduals(fit.cfaMod_saraya1)
fitMeasures(fit.cfaMod_saraya1)
#

# 
residus_standardises_saraya1 <- resid(fit.cfaMod_saraya1, type = "standardized")

#
residus_standardises_saraya1$cov

# 
inspect(fit.cfaMod_saraya1, what = "std")$lambda
# 
inspect(fit.cfaMod_saraya1, "std")$psi




#==========================
Mod_saraya<-'f1=~windgust_max+windspeed_max
f2=~moonphase_max
f3=~ precip_max+precipcover_max+cloudcover_min
f4=~windgust_min+windspeed_min 
f5=~sealevelpressure_min+sealevelpressure_max
cas_palu_confirme~ f1+f2+f4+f3'


mod.1_sar <- sem(Mod_saraya, data =data_saraya1_scaled,fixed.x=T,estimator = "MLR",std.lv=T)
fitMeasures(mod.1_sar, fit.measures = c("rmsea", "cfi", "tli", "srmr", "aic", "bic"))
mod.1_sar.fit=varTable(mod.1_sar)
library(lavaanPlot)
lavaanPlot(model=mod.1_sar, coefs=T, sig=0.01)
lavaanPlot(model=mod.1_sar, coefs=T,covs=T,stars='regress',stand=TRUE,sig=0.01, node_options = list(fontsize = 70),
           edge_options = list(fontsize = 80, penwidth = 6))

lavaanPlot(model = mod.1_sar, coefs = TRUE, covs = TRUE, stand = TRUE,sig=0.01)

semPaths(mod.1_sar,what="paths",whatLabels="par" )
semPaths(mod.1_sar,what="path",whatLabels="par",style = "ram",layout = "tree",
         rotation=2,sizeMan = 7,sizeLat = 7,color = "lightgray",edge.label.cex = 1.2,label.cex=1.3)
summary(mod.1_sar)
summary(mod.1_sar, stand=TRUE,fit.measures=TRUE,rsq=TRUE)





#===============================================================================



library(semPlot)
library(RColorBrewer)


node_labels <- c(  "windgust_min" = "Wind Gust Min",
                   "windspeed_min" = "Wind Speed Min",
                   "precip_max"="Precipitation Max",
                   "precipcover_max" = "Rain Cover Max",
                   "cloudcover_min" = "Cloud Cover Min",
                   "windgust_max" = "Wind Gust Max",
                   "windspeed_max" = "Wind Speed Max",
                   "sealevelpressure_min" = "Sea Pressure Min",
                   "sealevelpressure_max" = "Sea Pressure Max",
                   "moonphase_mean" = "Moon Phase Mean",
                   "cas_palu_confirme" = "Malaria cases \n confirmed by RDT",
                   "f4" = "F4",
                   "f3" = "F3", 
                   "f1" = "F1",
                   "f5" = "F5",
                   "f2" = "F2"
)

Mod_saraya<-'f4=~windgust_min+windspeed_min 
f3=~ precip_max+precipcover_max+cloudcover_min
f1=~windgust_max+windspeed_max
f5=~sealevelpressure_min+sealevelpressure_max
f2=~moonphase_max
cas_palu_confirme~ f1+f2+f4+f3'





#===============================================================
#===============================================================
par(bg = "white") ##
semPaths(mod.1_sar,
         what = "paths",
         whatLabels = "par",
         layout = "tree2",
         rotation = 2,
         
         # 
         color = list(man = "burlywood1", lat = "grey"),
         border.color = "#FF8800",
         edge.color = "black",
         
         # 
         sizeMan = 17,
         sizeLat = 10,
         sizeMan2 =6,
         sizeLat2 = 10,
         
         # 
         label.cex = 1.7,
         edge.label.cex = 1.4,
         
         # 
         nodeLabels = node_labels,
         nCharNodes = 0,
         nCharEdges = 0,
         nDigits =2,
         
         # 
         style = "ram",
         mar = c(4,3.5,4, 3.5),
         curve =2.3,
         asize =3,
         
         # 
         edge.label.bg = FALSE,
         optimizeLatRes = TRUE,
         residuals = TRUE,
         intercepts = FALSE,
         
         # 
         border.width = 3,
         minimum = 0,
         fade =TRUE,
         posCol = "black",
         negCol = "red"
)





