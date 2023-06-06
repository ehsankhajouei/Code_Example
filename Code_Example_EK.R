# Author "Ehsan Khajouei"
# Modified on 06June2023


# Cleaning the environment______________________________________________________
rm(list = ls()); cat("\f")


# Loading libraries_____________________________________________________________
# library(lme4)
library(lmerTest) # lmerTest is an expansion to lme4.
library(emmeans)
library(openxlsx)


# Reading in dataset with raw phenotypes and performing ANOVA___________________
physio_data <- read.csv("./raw_pheno.csv", header=TRUE)
physio_data[,c(1:4)] <- lapply(physio_data[,c(1:4)], factor) # Measurement_date is being treated as an integer.
sapply(physio_data, class) # Verifying the class of each variable. 
names(physio_data)

pheno_anova <- aov(Phenotype ~ Column + Row + Repetition + Measurement_Date + Taxa, data = physio_data) # ANOVA
anova(pheno_anova) # The effect of Measurement_Date is significant among all the factors. Its P-value is lower than 0.05. 


# BLUE__________________________________________________________________________
pheno_blue_1 <- lmer(Phenotype ~ Taxa + Measurement_Date + (1|Repetition:Column) + (1|Repetition:Row), data = physio_data)

anova(pheno_blue_1) # Measurement_Date was imported from the previous ANOVA step to the mixed model, but appeared to be insignificant after being tested
# under the new model, so it was removed from the mixed model. The corresponding P-value is 0.93 which is higher than 0.05.

summary(pheno_blue_1) # We can see in the 'Random effects' section, that 'Repetition:Row' contains zero variance, so we can drop it from the model. 
# But, before that, we can compare this model with the reduced model (not having 'Repetition:Row' as a random effect) to compare AIC.  


pheno_blue <- lmer(Phenotype ~ Taxa + Measurement_Date + (1|Repetition:Column), data = physio_data) # This is the reduced model.

anova(pheno_blue_1, pheno_blue) # Comparison between the models and having AIC criteria that are very similar, so adding another random part to the model 
# (in this case, 'Repetition:Row') did not improve the model and it is safe to drop it. The simpler model, which is 'pheno_blue', has a lower AIC which 
# is desirable. 
# Another drawback of the more complex model, i.e. 'pheno_blue_1', is that we face singulairty in the model which is unwanted. 
# Moreover, P-value of the likelihood ratio test is 1 which is high and indicates the addition of another random term did not improve the model.
anova(pheno_blue)
summary(pheno_blue)

pheno_blue_final <- lmer(Phenotype ~ Taxa + (1|Repetition:Column), data = physio_data) # The reduced model without the insignificant Measurement_Date term; 
# this model will be applied to predicting BLUEs.
summary(pheno_blue_final)


# Extracting BLUEs and forming its Data frame____________________________________
blue_final <- data.frame(emmeans(pheno_blue_final, 0~Taxa))

blue_f <- blue_final[,1:2]
colnames(blue_f)[2] <- "Adjusted_Pheno"
write.csv(blue_f, "./adjusted_pheno.csv", row.names=FALSE)


# Model and Normality Check______________________________________________________
physio_data <- na.omit(physio_data)
physio_data$Predicted_Pheno <- predict(pheno_blue_final, re.form=~0)
physio_data$Predicted_Pheno_Residual <- resid(pheno_blue_final, scaled=TRUE)

plot(physio_data$Predicted_Pheno_Residual ~ physio_data$Predicted_Pheno) # This plot checks for equal variance of residuals 
# (Level 1 residuals or the lowest level). 
# Having a systematic curvature or a funnel shape in this plot violates the assumption for residuals and the model. But, we see a cloud of spots which 
# shows constant variance of residuals, so the model assumption is valid. 

qqnorm(physio_data$Predicted_Pheno_Residual) # This Q-Q plot checks if the residuals are normally distributed. 
# Here, the normal Q-Q plot shows there is no outlier as residuals form a straight line which represents normally distributed residuals and the
# assumption of model is verified. 


# Performing GWAS________________________________________________________________

# Initializing Gapit by importing necessary libraries and then sourcing the main function either by using its data repository at:
# "http://zzlab.net/GAPIT/gapit_functions.txt" or saving the function locally and sourcing it. I used the version that I saved
# in July2022 which worked properly, therefore I source it from the local file. 

library(MASS)
library(multtest)
library(gplots)
library(compiler)
library(scatterplot3d) 
library(LDheatmap)  # To install LDheatmap package I had to first save it from https://cran.r-project.org/src/contrib/Archive/LDheatmap/
# and then I installed the version 0.99-4 using 'install.packages("./LDheatmap_0.99-4.tar.gz", repos=NULL, type="source")'. 

library(genetics)
library(ape)
library(EMMREML)
library(multtest)
library(snpStats)
library(bigmemory)

source("./gapit_functions_July2022.text")

# Reading in phenotypic and genotypic file. I used a hapmap fortmat of the SNP data. The data process to produce this hapmap file was performed by
# GATK and BCFtools using Linux environment and commmand line. This process involved:
# 1) Loci with missing values higher than 20% were removed
# 2) Loci with minor allele frequency less than 5% were removed
# 3) The filtered file was imputed using the Beagle package

myY = read.table("./adjusted_pheno.csv", sep = ",", header=TRUE)
myG = read.delim("./maf5_missed20_imputed_96.hmp.txt", head = FALSE)
#myKI <- read.table("./GAPIT.Kin.VanRaden.csv", head = FALSE, sep = ",") # GAPIT has the option to include a user defined kinship file into the main function.
# In this example, I asked GAPIT to calculate the kinship matrix based on the VanRaden method and commented out 'myKI' option in the main function. 

# Calling the main function of GAPIT. Here, I used 7 principal components (PCs) based on the previous results with different PC numbers, visualizing the 
# eigenvalue plot, and also checking the QQ plots of each GWAS based on a certain PC. The QQ plot shows false discoveries and it is desirable to select 
# a PC number that leads to less false discoveries. 
# Based on previous trials, the MLMM method shows to be an optimal method in controlling false discoveries. Therefore, here I am using MLMM and MLM which 
# is the most basic mixed-model in GAPIT that takes into account the population structure. 

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  kinship.algorithm = "VanRaden",
  #KI=myKI,
  model = c("MLM", "MLMM"),
  PCA.total = 7,
  Multiple_analysis = TRUE,
  file.output=TRUE
)

# Following GWAS and discovering marker-trait association (MTA), LD analysis and locating the associated loci, and gene ontology were performed.
