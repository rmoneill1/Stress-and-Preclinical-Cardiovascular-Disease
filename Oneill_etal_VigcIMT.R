setwd("C:/Users/Rmoneill/Box/SUMMER 2025")
require(psych)
require(tidyverse)
require(tidyr)
require(Hmisc)
require(car)
require(moments)
require(MASS)
require(DescTools)
library(mice)
library(pan)
library(ggplot2)
library(dplyr)
require(DescTools)
require(naniar)
require(scipub)
require(MASS)
library(lattice)
library(effects)
library(ggeffects)
library(emmeans)


#Read in data---- 
    f <- read.csv("Oneill_etal_VigcIMT.csv", na.strings = c("-999" , "", "NA"))
    
#Data prep (run every time) ----

  #Ensure necessary numeric/factor variables are stored accordingly  
      vars_to_convert <- c(
        "mEMASVQTV", "socsaf_mean", "MMFCCAV1T1", "MMBIFV1T1", "MMICAV1T1", "MMBIFICAV1T1",
        "MMFCCAV1T2", "MMBIFV1T2", "MMICAV1T2", "MMBIFICAV1T2", "CMMFCCAV1", "CMMBIFV1",
        "CMMICAV1", "CCBIFICAV1", "HIncome", "BMIT1", "SBPT1", "DBPT1", "Age"
      )
      
      f[vars_to_convert] <- lapply(f[vars_to_convert], as.numeric)
      
      f$Sex <- factor(f$Sex, levels=c('1','2'), labels=c('Male', 'Female'))
      f$Education<- factor(f$Education)
      f$Jobclass<- factor(f$Jobclass)
      f$CardiaxHxT1 <- factor(f$CardiaxHxT1)
      f$Marital <- factor(f$Marital)
      f$MedLipT1 <- factor(f$MedLipT1)
      f$MedDT1 <- factor(f$MedDT1)
      f$MedBT1 <- factor(f$MedBT1)
      f$Shift <- factor(f$Shift)
      f$ref <- factor(f$ref)
      f$ref <- relevel(f$ref, ref = "NHB")
      f$refNHW <- relevel(f$ref, ref = "NHW")
      f$HIncome_3level <- factor(f$HIncome_3level)
      f$HIncome_3level <- as.character(f$HIncome_3level)
      f$HIncome_3level[f$HIncome_3level == "Missing"] <- NA  
      f$HIncome_3level <- factor(f$HIncome_3level)  
      f$HIncome_3level <- relevel(f$HIncome_3level, ref = "<$30,000")
      
      #Label original Income variable 
          income_labels <- c(
            "Less than $10,000",
            "$10,001 - $15,000",
            "$15,001 – $20,000",
            "$20,001 - $25,000",
            "$25,001 - $29,999",
            "$30,000 - $40,000",
            "$40,001 - $50,000",
            "$50,001 - $75,000",
            "$75,001 - $100,000",
            "$100,001 - $150,000",
            "$150,001 - $200,000",
            "$200,001 - $250,000",
            "More than $250,000"
          )
        
          f$HIncome_f <- factor(f$HIncome, 
                              levels = 1:13, 
                              labels = income_labels)
      

      #Variable prep for mediation models 
      # Create dummy variables for ref
      dummies <- model.matrix(~ ref - 1, data = f)
      # Rename the columns for clarity
      colnames(dummies) <- c("ref_NHB", "ref_HL"  , "ref_NHO", "ref_NHW")
      # Add the dummy variables to the dataset
      f <- cbind(f, dummies)
      
      #create dummy variables for income(3level)
      # Convert NA values in HIncome_3level to "Missing"
      f$HIncome_3level_med <- factor(ifelse(is.na(f$HIncome_3level), "Missing", as.character(f$HIncome_3level)),
                                     levels = c("<$30,000", "[$30,000, $100,000]", ">$100,000", "Missing"))
      # Create dummy variables using model.matrix()
      dummies2 <- model.matrix(~ HIncome_3level_med - 1, data = f)
      
      # Rename the columns for clarity
      colnames(dummies2) <- c("Income_below_30K", "Income_30K_100K", "Income_above_100K", "Income_Missing")
      # Combine the dummy variables with the original dataset
      f <- cbind(f, dummies2)
      
      
      #Creating necessary interaction terms
      f$Sex_num <- ifelse(f$Sex == "Male", 1,
                          ifelse(f$Sex == "Female", 2, NA))
      f$Crime1z_Sex <- f$TotalCrime_1z * f$Sex_num  
      f$Crime5z_Sex <- f$TotalCrime_5z * f$Sex_num  
      f$Crime10z_Sex <- f$TotalCrime_10z * f$Sex_num  
      f$Crime1z_Age <- f$TotalCrime_1z * f$Age
      f$Crime5z_Age <- f$TotalCrime_5z * f$Age
      f$Crime10z_Age <- f$TotalCrime_10z * f$Age
      f$Age_Sex <- f$Age * f$Sex_num
      f$Crime1z_Sex_Age <- f$TotalCrime_1z * f$Sex_num * f$Age
      f$Crime5z_Sex_Age <- f$TotalCrime_5z * f$Sex_num  * f$Age
      f$Crime10z_Sex_Age <- f$TotalCrime_10z * f$Sex_num * f$Age
      
      f$Sex_num <- ifelse(f$Sex == "Male", 1,
                          ifelse(f$Sex == "Female", 2, NA))
      f$socsaf_Sex <- f$socsaf_mean * f$Sex_num  
      f$socsaf_Age <- f$socsaf_mean * f$Age 
      f$Age_Sex <- f$Age * f$Sex_num
      f$socsaf_Sex_Age <- f$socsaf_mean * f$Age * f$Sex_num
   
  #When analyzing crime variables, use subset g with participant with non-TX address removed: 
      g <- subset(f, !f$ID==5)
    

#Missingness ----
    ##Little's MCAR test <Tierney & Cook (2018); Little & Roderick (1988)> 
      install.packages("remotes")
      remotes::install_github("njtierney/naniar")
       require(naniar)
      #compares the estimated means of each variable between the different missing patterns. 

  #Initially assessed patterns of missingness among all main/substantive variables together. 
  #Data were not MCAR when examined all together. However, further examination demonstrated that cIMT were responsible for the systematic patterns of missingness. 
  #So, examined those variables separately for ease of understanding trends in missing data. 

#Missingness among main/substantive variables without cIMT variables 
    u <- subset(f, select = c(TotalCrime_1, TotalCrime_5, TotalCrime_10, mEMASVQTV, socsaf_mean))
    propmiss <- function(u) 
    {lapply(u,function(x) 
      data.frame(nmiss=sum(is.na(x)), n=length(x), propmiss=sum(is.na(x))/length(x)))} 
    require(mice)
    propmiss(u) 
    md.pattern(u)
    
    #Creating missingness indicators 
    f$mEMASVQTV_na <- ifelse(is.na(f$mEMASVQTV) == TRUE, 1, 0) #if case is missing = TRUE (1) ; if case is not missing = FALSE (0)
    table(f$mEMASVQTV_na)
    
    f$socsaf_mean_na <- ifelse(is.na(f$socsaf_mean) == TRUE, 1, 0)
    table(f$socsaf_mean_na)
    
  #Checking predictors of missingness for SVQ & social safety variables 
    #socsaf_mean_na
    m <- glm(socsaf_mean_na ~ Sex, data = f, family = binomial())
    summary(m)
    
    m <- glm(socsaf_mean_na ~ Age, data = f, family = binomial())
    summary(m)
    
    #mEMASVQTV_na
    m <- glm(mEMASVQTV_na ~ Sex, data = f, family = binomial())
    summary(m)
    
    m <- glm(mEMASVQTV_na ~ Age, data = f, family = binomial())
    summary(m)
      
    
#Missingness for cimt variables
  #Attrition was related with some of the missingness for Time 2 cIMT measurement & change scores 
      table(f$ParticipateT2) #61 people did not participate in Time 2 
      
  #Examining predictors of likelihood to return at time 2:
      #Sex
      m <- glm(ParticipateT2 ~ Sex, data = f, family = binomial())
      summary(m)   
      exp(m$coefficients) 
      exp(confint(m))
      
      #Age
      m <- glm(ParticipateT2 ~ Age, data = f, family = binomial())
      summary(m)
      exp(m$coefficients) 
      describeBy(f$Age, f$ParticipateT2) 
      exp(confint(m))
      
      #Household Income
      m <- glm(ParticipateT2 ~ HIncome_3level, data = f, family = binomial())
      summary(m)
      
      #Household Income (Releveled)
      m <- glm(ParticipateT2 ~ HIncome_3relevel, data = f, family = binomial())
      summary(m)

      #Household Income (controlling for BMI)
      m <- glm(ParticipateT2 ~ HIncome_3level +BMIT1, data = f, family = binomial())
      summary(m)
      exp(m$coefficients) 
      
      #Household Income (controlling for sex)
      m <- glm(ParticipateT2 ~ HIncome_3level + Sex, data = f, family = binomial())
      summary(m) 
      exp(m$coefficients)
      exp(confint(m))
      
      #Race/ethnicity 
      m <- glm(ParticipateT2 ~ refNHW, data = f, family = binomial())
      summary(m)
      
      #Caregiving status 
      m <- glm(ParticipateT2 ~ Caregiver, data = f, family = binomial())
      summary(m)
      
      #Educational attainment  
      m <- glm(ParticipateT2 ~ Education, data = f, family = binomial())
      summary(m)

      #Job Class
      m <- glm(ParticipateT2 ~ Jobclass, data = f, family = binomial())
      summary(m)
      
      #BMI
      m <- glm(ParticipateT2 ~ BMIT1, data = f, family = binomial())
      summary(m)
      
      #Cardiac History 
      m <- glm(ParticipateT2 ~ CardiaxHxT1, data = f, family = binomial())
      summary(m)
      
      #Neuroticism
      m <- glm(ParticipateT2 ~ NT1, data = f, family = binomial())
      summary(m) 
      
      #Currently smoking cigarettes 
      m <- glm(ParticipateT2 ~ Tobacco7, data = f, family = binomial())
      summary(m) 
      
      #Marital status
      m <- glm(ParticipateT2 ~ Marital, data = f, family = binomial())
      summary(m) 
      
      #Depressive symptoms 
      m <- glm(ParticipateT2 ~ cesdtott1, data = f, family = binomial())
      summary(m) 
      
      #Social vigilance 
      m <- glm(ParticipateT2 ~ SVQSTOTT1, data = f, family = binomial())
      summary(m) 
      
      #Shift work status 
      m <- glm(ParticipateT2 ~ Shift, data = f, family = binomial())
      summary(m)
      
#MCAR test on data without missing participants from non-returners 
      table(f$ParticipateT2) #61 people did not return 
      t2 <- subset(f, ParticipateT2=="1")
      ut2 <- subset(t2, select = c(TotalCrime_1, TotalCrime_5, TotalCrime_10, mEMASVQTV, socsaf_mean, CMMFCCAV1 , CMMBIFV1 , CMMICAV1 , CCBIFICAV1))
      
      mcar_test(ut2) #test significant - data are still not missing completely at random 
      
  #Assessing missingness patterns among cimt variables 
      c <- subset (f, select = c(MMFCCAV1T1, MMBIFV1T1, MMICAV1T1, MMBIFICAV1T1, MMFCCAV1T2, MMBIFV1T2, MMICAV1T2, MMBIFICAV1T2, CMMFCCAV1, CMMBIFV1, CMMICAV1, CCBIFICAV1))
      propmiss(c)
      md.pattern(c, rotate.names=TRUE) #Only 115 (38.33%) of all participants had complete cimt data 
      
      #MMFCCA variables 
      fcca <- subset (f, select = c(MMFCCAV1T1, MMFCCAV1T2, CMMFCCAV1))
      propmiss(fcca)
      md.pattern(fcca, rotate.names=TRUE)
   
      #MMBIF variables 
      bif <- subset (f, select = c(MMBIFV1T1, MMBIFV1T2, CMMBIFV1 ))
      propmiss(bif)
      md.pattern(bif, rotate.names=TRUE)

      #MMICA variables 
      ica <- subset (f, select = c(MMICAV1T1, MMICAV1T2, CMMICAV1))
      propmiss(ica)
      md.pattern(ica, rotate.names=TRUE)

      #MMBIFICA variables 
      comb <- subset (f, select = c(MMBIFICAV1T1, MMBIFICAV1T2, CCBIFICAV1))
      propmiss(comb)
      md.pattern(comb, rotate.names=TRUE)

  
 #cIMT data are not MCAR, therefore will compute variables to be indicators of missingness for assessing mechanisms of missingness: 

    f$CMMFCCAV1_na <- ifelse(is.na(f$CMMFCCAV1) == TRUE, 1, 0)
    table(f$CMMFCCAV1_na)
    
    f$CMMBIFV1_na <- ifelse(is.na(f$CMMBIFV1) == TRUE, 1, 0)
    table(f$CMMBIFV1_na)
    
    f$CMMICAV1_na <- ifelse(is.na(f$CMMICAV1) == TRUE, 1, 0)
    table(f$CMMICAV1_na)
    
    f$CCBIFICAV1_na <- ifelse(is.na(f$CCBIFICAV1) == TRUE, 1, 0)
    table(f$CCBIFICAV1_na)
    
  #Mechanisms of missing analysis - determining variables associated with likelihood of missing cIMT data 
    #CMMFCCAV1_na
        #Sex
        m <- glm(CMMFCCAV1_na ~ Sex, data = f, family = binomial())
        summary(m)
        exp(m$coefficients) 
        exp(coef(m))
        exp(confint(m))
        
        #Age
        m <- glm(CMMFCCAV1_na ~ Age, data = f, family = binomial())
        summary(m)
        exp(m$coefficients) 
        exp(coef(m))
        exp(confint(m))
        
        #Household income (3-level factor)
        m <- glm(CMMFCCAV1_na ~ HIncome_3level, data = f, family = binomial())
        summary(m) 
        exp(m$coefficients) 
        exp(coef(m))
        exp(confint(m))
        
        
        #Sex*Age interaction 
        m <- glm(CMMFCCAV1_na ~ Sex*Age, data = f, family = binomial())
        summary(m) 
        exp(m$coefficients) 
        exp(confint(m))
        
        #Household income*Sex interaction 
        m <- glm(CMMFCCAV1_na ~ HIncome_3level*Sex, data = f, family = binomial())
        summary(m)        
        exp(m$coefficients) 
        exp(confint(m))
        
        #Household income*Age 
        m <- glm(CMMFCCAV1_na ~ HIncome_3level*Age, data = f, family = binomial())
        summary(m)
        exp(m$coefficients) 
        exp(confint(m))
        
        #Household Income (controlling for BMI)
        m <- glm(CMMFCCAV1_na ~ HIncome_3level + BMIT1, data = f, family = binomial())
        summary(m) 
        exp(m$coefficients) 
        exp(confint(m))
        
        #Race/Ethnicity 
        m <- glm(CMMFCCAV1_na ~ refNHW, data = f, family = binomial())
        summary(m)
        exp(coef(m))
        exp(confint(m))
        
        
        #BMI
        m <- glm(CMMFCCAV1_na ~ BMIT1, data = f, family = binomial())
        summary(m)
        exp(coef(m))
        exp(confint(m))
          
        #Cardiac History 
        m <- glm(CMMFCCAV1_na ~ CardiaxHxT1, data = f, family = binomial())
        summary(m)
        exp(coef(m))
        exp(confint(m))  
        
  #CMMBIFV1_na
      
      #Sex 
      m <- glm(CMMBIFV1_na ~ Sex, data = f, family = binomial())
      summary(m)
      exp(m$coefficients) 
      exp(coef(m))
      exp(confint(m))  
      
      #Age
      m <- glm(CMMBIFV1_na ~ Age, data = f, family = binomial())
      summary(m)
      exp(m$coefficients) 
      exp(confint(m))  
      
      #BMI
      m <- glm(CMMBIFV1_na ~ BMIT1, data = f, family = binomial())
      summary(m)
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m)) 
      
      #Race/Ethnicity 
      m <- glm(CMMBIFV1_na ~ refNHW, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 
      
      #Household Income 
      m <- glm(CMMBIFV1_na ~ HIncome_3level, data = f, family = binomial())
      summary(m) 
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m)) 
      
      #Household income (controlling for BMI)
      m <- glm(CMMBIFV1_na ~ HIncome_3level + BMIT1, data = f, family = binomial())
      summary(m) 
      exp(coef(m))
      exp(confint(m)) 
      
      #Household income*Age interaction 
      m <- glm(CMMBIFV1_na ~ HIncome_3level*Age, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 
      
      #Household income*Sex interaction 
      m <- glm(CMMBIFV1_na ~ HIncome_3level*Sex, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 
      
      #Household income*BMI interaction 
      m <- glm(CMMBIFV1_na ~ HIncome_3level*BMIT1, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 
      
      #Sex*Age interaction 
      m <- glm(CMMBIFV1_na ~ Sex*Age, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 
      
      #Sex*BMI interaction 
      m <- glm(CMMBIFV1_na ~ Sex*BMIT1, data = f, family = binomial())
      summary(m) 
      exp(coef(m))
      exp(confint(m)) 
      
      #Age*BMI interaction 
      m <- glm(CMMBIFV1_na ~ Age*BMIT1, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 

      #Cardiac history 
      m <- glm(CMMBIFV1_na ~ CardiaxHxT1, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m)) 
      
    
  #CMMICAV1_na
      
      #BMI 
      m <- glm(CMMICAV1_na ~ BMIT1, data = f, family = binomial())
      summary(m)
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m)) 
      
      #Household income 
      m <- glm(CMMICAV1_na ~ HIncome_3level, data = f, family = binomial())
      summary(m)
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))
      
      #Household income (controlling for BMI)
      m <- glm(CMMICAV1_na ~ HIncome_3level + BMIT1, data = f, family = binomial())
      summary(m) 
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))

      #Sex 
      m <- glm(CMMICAV1_na ~ Sex, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Age
      m <- glm(CMMICAV1_na ~ Age, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Cardiac history 
      m <- glm(CMMICAV1_na ~ CardiaxHxT1, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Race/ethnicity 
      m <- glm(CMMICAV1_na ~ refNHW, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
  
  #CCBIFICAV1_na
      
      #BMI
      m <- glm(CCBIFICAV1_na ~ BMIT1, data = f, family = binomial())
      summary(m)
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))
      
      #Sex 
      m <- glm(CCBIFICAV1_na ~ Sex, data = f, family = binomial())
      summary(m) 
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))
      
      #Age 
      m <- glm(CCBIFICAV1_na ~ Age, data = f, family = binomial())
      summary(m)  
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))
      
      #Household income 
      m <- glm(CCBIFICAV1_na ~ HIncome_3level, data = f, family = binomial())
      summary(m) 
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))
      
      #Race/ethnicity 
      m <- glm(CCBIFICAV1_na ~ refNHW, data = f, family = binomial())
      summary(m)
      exp(m$coefficients)
      exp(coef(m))
      exp(confint(m))
      
      #Household income (controlling for BMI)
      m <- glm(CCBIFICAV1_na ~ HIncome_3level + BMIT1, data = f, family = binomial())
      summary(m) 
      exp(coef(m))
      exp(confint(m))
      
      #Household income*Age interaction 
      m <- glm(CCBIFICAV1_na ~ HIncome_3level*Age, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Household income*Sex 
      m <- glm(CCBIFICAV1_na ~ HIncome_3level*Sex, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Household income*BMI
      m <- glm(CCBIFICAV1_na ~ HIncome_3level*BMIT1, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #BMI*Sex
      m <- glm(CCBIFICAV1_na ~ BMIT1*Sex, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #BMI*Age
      m <- glm(CCBIFICAV1_na ~ BMIT1*Age, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Sex*Age
      m <- glm(CCBIFICAV1_na ~ Sex*Age, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      
      #Cardiac history
      m <- glm(CCBIFICAV1_na ~ CardiaxHxT1, data = f, family = binomial())
      summary(m)
      exp(coef(m))
      exp(confint(m))
      

#Outliers ----
#Assessed each substantive psychosocial variable to determine necessary transformations for addressing influential outliers
      
  #Used "find outlier" function with ggplot2 
      #(defined outlier as observation falls within the first quartile by 1.5 times the interquartile range (Q1) or exceeds the third quartile by 1.5 times the interquartile range (Q3) 
      
          findoutlier <- function(x) {
            return(x < quantile(x, .25, na.rm = TRUE) - 1.5*IQR(x, na.rm = TRUE) | x > quantile(x, .75, na.rm = TRUE) + 1.5*IQR(x, na.rm = TRUE))}
          
          #Social safety (6 outliers)
          describe(f$socsaf_mean)
          f$socsaf_mean_out <- ifelse(findoutlier(f$socsaf_mean), f$socsaf_mean, NA)
          describe(f$socsaf_mean_out)
          table(f$socsaf_mean_out)
          
          #Total Crime 1-mile (1 outlier)
          f$TotalCrime_1_out <- ifelse(findoutlier(f$TotalCrime_1), f$TotalCrime_1, NA)
          describe(f$TotalCrime_1_out)
          table(f$TotalCrime_1_out)

          #Total Crime 5-mile (6 outliers)
          f$TotalCrime_5_out <- ifelse(findoutlier(f$TotalCrime_5), f$TotalCrime_5, NA)
          describe(f$TotalCrime_5_out)
          table(f$TotalCrime_5_out)
          
          #Total Crime 10-mile (25 outliers)
          f$TotalCrime_10_out <- ifelse(findoutlier(f$TotalCrime_10), f$TotalCrime_10, NA)
          describe(f$TotalCrime_10_out)
          table(f$TotalCrime_10_out)
          
          #EMA-assessed social vigilance (23 outliers)
          f$mEMASVQTV_out <- ifelse(findoutlier(f$mEMASVQTV), f$mEMASVQTV, NA)
          describe(f$mEMASVQTV_out)
          table(f$mEMASVQTV_out)
          
          #FCCA (17 outliers)
          psych::describe(f$CMMFCCAV1)
          hist(f$CMMFCCAV1)  
          f$CMMFCCAV1_out <- ifelse(findoutlier(f$CMMFCCAV1), f$CMMFCCAV1, NA)
          describe(f$CMMFCCAV1_out)
          table(f$CMMFCCAV1_out) 
          plot(f$CMMFCCAV1_out)
          
          #BIF (7 outliers)
          plot(f$CMMBIFV1)
          f$CMMBIFV1_out <- ifelse(findoutlier(f$CMMBIFV1), f$CMMBIFV1, NA)
          describe(f$CMMBIFV1_out)
          table(f$CMMBIFV1_out) 
          plot(f$CMMBIFV1_out)
          
          #ICA (6 outliers)
          plot(f$CMMICAV1) 
          f$CMMICAV1_out <- ifelse(findoutlier(f$CMMICAV1), f$CMMICAV1, NA)
          describe(f$CMMICAV1_out)
          table(f$CMMICAV1_out) 
          plot(f$CMMICAV1_out)
          
          #BIF/ICA (13 outliers)
          plot(f$CCBIFICAV1) 
          f$CCBIFICAV1_out <- ifelse(findoutlier(f$CCBIFICAV1), f$CCBIFICAV1, NA)
          describe(f$CCBIFICAV1_out)
          table(f$CCBIFICAV1_out) #13 outliers for change in combined BIF/ICA V1 
          plot(f$CCBIFICAV1_out)
          
  
      
      #Creating necessary transformed variables 
        #Z-scored Total Crime 1-Mile 
          describe(f$TotalCrime_1)
          #f$TotalCrime_1z <- scale(f$TotalCrime_1, center = TRUE, scale = TRUE)
          #f$TotalCrime_1z <- as.numeric(f$TotalCrime_1z)
          #describe(f$TotalCrime_1z)
      
        #Z-scored Total Crime 5-Mile 
          describe(f$TotalCrime_5)
          #f$TotalCrime_5z <- scale(f$TotalCrime_5, center = TRUE, scale = TRUE)
          #f$TotalCrime_5z <- as.numeric(f$TotalCrime_5z)
          #describe(f$TotalCrime_5z)
        
        #Z-scored Total Crime 10-Mile 
          describe(f$TotalCrime_10)
          #f$TotalCrime_10z <- scale(f$TotalCrime_10, center = TRUE, scale = TRUE)
          #f$TotalCrime_10z <- as.numeric(f$TotalCrime_10z)
          #describe(f$TotalCrime_10z)
      
        #Box-Cox transformed mean EMA Total social vigilance 
          #boxcox(lm(f$mEMASVQTV ~ 1))
          
          #b <- boxcox(lm(f$mEMASVQTV ~ 1))
          #lambda <- b$x[which.max(b$y)]
          #lambda # -2 
          
          #f$mEMASVQTV_boxcox <- (f$mEMASVQTV ^ (-2) - 1) / -2
          hist(f$mEMASVQTV_boxcox)
          #psych::describe(f$mEMASVQTV_boxcox) #skew = 0.37 
          #boxplot(f$mEMASVQTV_boxcox, data = f, horizontal = TRUE)
          
        #Winsorized change in FCCA V1
          #f$CMMFCCAV1_w <- Winsorize(f$CMMFCCAV1, minval = NULL, maxval = NULL, probs = c(0.01, 0.999), na.rm = TRUE, type = 7)
        #Winsorized change in BIF V1
          #f$CMMBIFV1_w <- Winsorize(f$CMMBIFV1, minval = NULL, maxval = NULL, probs = c(0.01, 0.99), na.rm = TRUE, type = 7)
        #Winsorized change in ICA V1
          #f$CMMICAV1_w <- Winsorize(f$CMMICAV1, minval = NULL, maxval = NULL, probs = c(0.01, 0.999), na.rm = TRUE, type = 7)
        #Winsorized change in BIF/ICA V1
          #f$CCBIFICAV1_w <- Winsorize(f$CCBIFICAV1, minval = NULL, maxval = NULL, probs = c(0.01, 0.999), na.rm = TRUE, type = 7)

    
 
#Covariates (& Differences by Age and Sex)----

  #Hypothesis 1: Environmental safety positively related to social safety 
      #Candidate covariates: Age, Sex, Race, Ethnicity, Income 
      #Relevant covariates after analyses (see below): Age, Income 

    #Correlations between environmental safety, social safety, sv, & age  
      h1cov <- subset (f, select = c(TotalCrime_1z, TotalCrime_5z, TotalCrime_10z, socsaf_mean, mEMASVQTV_boxcox, Age))
      print((corr.test(h1cov)), short=FALSE)
      corPlot(h1cov, stars = TRUE)
      
    #Crime risk by age 
      m1 <- lm(TotalCrime_1z ~ Age, data = g)
      summary(m1)
      m2 <- lm(TotalCrime_5z ~ Age, data = g)
      summary(m2)
      m3 <- lm(TotalCrime_10z ~ Age, data = g)
      summary(m3)
      
    #Social safety by age 
      m4 <- lm(socsaf_mean ~ Age, data = f)
      summary(m4)

    #Safety differences by sex 
      t.test(TotalCrime_1z ~ Sex, data = g) 
      t.test(TotalCrime_5z ~ Sex, data = g) 
      t.test(TotalCrime_10z ~ Sex, data = g) 
      t.test(socsaf_mean ~ Sex, data = f) 
      t.test(mEMASVQTV_boxcox ~ Sex, data = f) 

  #Safety differences by race/ethnicity
      anova_result <- aov(f$socsaf_mean ~ f$ref)
      summary(anova_result) 
      TukeyHSD(anova_result) 
      
      anova_result <- aov(f$TotalCrime_1z ~ f$ref)
      summary(anova_result) 
      TukeyHSD(anova_result) 
      
      anova_result <- aov(f$TotalCrime_5z ~ f$ref)
      summary(anova_result) 
      TukeyHSD(anova_result) 
      
      anova_result <- aov(f$TotalCrime_10z ~ f$ref)
      summary(anova_result) 
      TukeyHSD(anova_result) 
      
      anova_result <- aov(f$mEMASVQTV_boxcox ~ f$ref)
      summary(anova_result) 
 
      
  #Safety differences by income 
      anova_result <- aov(f$TotalCrime_1z ~ f$HIncome_3level)
      summary(anova_result) 
      TukeyHSD(anova_result)
      pairwise.t.test(f$TotalCrime_1, f$HIncome_3level, p.adjust.method = "bonferroni")
      
      
      anova_result <- aov(f$TotalCrime_5z ~ f$HIncome_3level)
      summary(anova_result) 

      anova_result <- aov(f$TotalCrime_10z ~ f$HIncome_3level)
      summary(anova_result) 

      anova_result <- aov(f$socsaf_mean ~ f$HIncome_3level)
      summary(anova_result)  
      
      anova_result <- aov(f$mEMASVQTV_boxcox ~ f$HIncome_3level)
      summary(anova_result) 
      
  #Hypothesis 2: Environmental safety will be negatively associated with social vigilance behavior in daily life 
      #Candidate covariates: Age, Sex, Race, Ethnicity, Income 
      #Relevant covariates after analyses (see below): Age, Sex, HIncome 
      
      #Correlations between environmental safety, social vigilance, & age 
      h2cov <- subset (f, select = c(TotalCrime_1z, TotalCrime_5z, TotalCrime_10z, mEMASVQTV_boxcox, Age))
      print((corr.test(h2cov)), short=FALSE)
      corPlot(h2cov, stars = TRUE)
      
      #Social vigilance by age 
      m1 <- lm(mEMASVQTV_boxcox ~ Age, data = f)
      summary(m1)
     
      #Social vigilance differences by sex 
      t.test(mEMASVQTV_boxcox ~ Sex, data = f) 
      tapply(f$mEMASVQTV_boxcox, f$Sex, sd, na.rm = TRUE)
      result <- t.test(mEMASVQTV_boxcox ~ Sex, data = f) 
      result$estimate 
      diff(result$estimate)  

      #Social vigilance differences by race/ethnicity
      anova_result <- aov(f$mEMASVQTV_boxcox ~ f$ref)
      summary(anova_result)

      #Social vigilance differences by Income 
      anova_result <- aov(f$mEMASVQTV_boxcox ~ f$HIncome_3level)
      summary(anova_result)
           
   
  #Hypothesis 3: Social safety will be negatively associated with social vigilance behavior in daily life  
      #Relevant covariates after analyses (see above): Age, Sex 

  #H4: Social vigilance behavior at time 1 will be associated with greater increases in cIMT at 2-year follow-up
      #Candidate covariates: Age, Sex, Race, Ethnicity, Income , BMI, Resting SBP, Resting DBP, Time 1 cIMT outcome, menopausal status (MenopT1 - 0 = No, 1 = Yes), blood pressure meds (MedBT1), lipid meds (MedLipT1), other cardiac meds (MedCT1), diabetes meds (MedDT1)
      #Relevant covariates after analyses (see below): 
        #Social vigilance: Age, Sex 
        #Far common carotid artery:Age, SBP, DBP, T1 FCCA, Lipid meds  
        #Bifurcation: sex, time 1 bif
        #Internal Carotid Artery: time 1 ICA, lipid meds, diabetes meds
        #BIF/ICA: DBP, time 1 BIF/ICA

    #FCCA
      #Correlations between social vigilance, change in fcca, age, BMI, resting SBP, resting DBP, time 1 fcca 
      h4cov <- subset (f, select = c(mEMASVQTV_boxcox, CMMFCCAV1_w, Age, BMIT1, SBPT1, DBPT1, MMFCCAV1T1 ))
      print((corr.test(h4cov)), short=FALSE)
      corPlot(h4cov, stars = TRUE)
      
      #Change in FCCA differences by Age 
      m1 <- lm(CMMFCCAV1_w ~ Age, data = f)
      summary(m1)
  
      #Change in FCCA differences by Sex 
      t.test(CMMFCCAV1_w ~ Sex, data = f) 

      #Change in FCCA differences by Ethnicity 
      t.test(CMMFCCAV1_w ~ Ethnicity, data = f) 

      #Change in FCCA differences by menopausal status 
      f$Sex_num <- as.numeric(f$Sex)
      f$MenopT1_V2 <- ifelse(f$Sex_num == 1, 3, 
                             ifelse(f$MenopT1 == 0, 0, 
                                    ifelse(f$MenopT1 == 1, 1, NA))) #Creating new variable where 0 = no menopause (based on avg age of 0's), 1 = menopause, 3 = men, all others are NA 
      table(f$MenopT1_V2)
      f$MenopT1_V2 <- as.factor(f$MenopT1_V2)
      
      anova_result <- aov(f$CMMFCCAV1_w ~ f$MenopT1_V2)
      summary(anova_result) 

      #Change in FCCA differences by blood pressure meds 
      t.test(CMMFCCAV1_w ~ MedBT1, data = f)

      #Change in FCCA differences by lipid meds   
      t.test(CMMFCCAV1_w ~ MedLipT1, data = f)
      tapply(f$CMMFCCAV1_w, f$MedLipT1, sd, na.rm = TRUE)

      #Change in FCCA differences by other cardiac meds   
      t.test(CMMFCCAV1_w ~ MedCT1, data = f) 
      
      #Change in FCCA differences by diabetes meds  
      t.test(CMMFCCAV1_w ~ MedDT1, data = f)

      #Change in FCCA differences by Race/Ethnicity
      anova_result <- aov(f$CMMFCCAV1_w ~ f$ref)
      summary(anova_result) 

      #Change in FCCA differences by Income
      anova_result <- aov(f$CMMFCCAV1_w ~ f$HIncome_3level)
      summary(anova_result) 

    
    #BIF
      #Correlations between social vigilance, change in bif, age, BMI, resting SBP, resting DBP, time 1 bif
      h4covb <- subset (f, select = c(mEMASVQTV_boxcox, CMMBIFV1_w, Age, BMIT1, SBPT1, DBPT1, MMBIFV1T1 ))
      print((corr.test(h4covb)), short=FALSE)
      corPlot(h4covb, stars = TRUE)
      
      #Change in BIF differences by Age 
      m1 <- lm(CMMBIFV1_w ~ Age, data = f)
      summary(m1)

      #Change in BIF differences by Sex 
      t.test(CMMBIFV1_w ~ Sex, data = f)
      tapply(f$CMMBIFV1_w, f$Sex, sd, na.rm = TRUE)
      result <- t.test(CMMBIFV1_w ~ Sex, data = f)
      result$estimate 
      diff(result$estimate)  
      
      #Change in BIF differences by menopausal status 
      anova_result <- aov(f$CMMBIFV1_w ~ f$MenopT1_V2)
      summary(anova_result)
      TukeyHSD(anova_result) #Menopausal status not included in covariates because it is driven by sex difference (already included in covariates)
      
      #Change in BIF differences by blood pressure meds 
      t.test(CMMBIFV1_w ~ MedBT1, data = f)

      #Change in BIF differences by lipid meds   
      t.test(CMMBIFV1_w ~ MedLipT1, data = f) 

      #Change in BIF differences by other cardiac meds   
      t.test(CMMBIFV1_w ~ MedCT1, data = f) 
      
      #Change in BIF differences by diabetes meds   
      t.test(CMMBIFV1_w ~ MedDT1, data = f)

      #Change in BIF differences by Race/Ethnicity
      anova_result <- aov(f$CMMBIFV1_w ~ f$ref)
      summary(anova_result) 

      #Change in BIF differences by Income
      anova_result <- aov(f$CMMBIFV1_w ~ f$HIncome_3level)
      summary(anova_result) 

    
    #ICA
      #Correlations between social vigilance, change in ICA, age, BMI, resting SBP, resting DBP, time 1 ICA 
      h4covi <- subset (f, select = c(mEMASVQTV_boxcox, CMMICAV1_w, Age, BMIT1, SBPT1, DBPT1, MMICAV1T1 ))
      print((corr.test(h4covi)), short=FALSE)
      corPlot(h4covi, stars = TRUE)

      #Change in ICA differences by Sex 
      t.test(CMMICAV1_w ~ Sex, data = f) 

      #Change in ICA differences by menopausal status 
      anova_result <- aov(f$CMMICAV1_w ~ f$MenopT1_V2)
      summary(anova_result) 

      #Change in ICA differences by blood pressure meds 
      t.test(CMMICAV1_w ~ MedBT1, data = f) 

      #Change in ICA differences by lipid meds   
      t.test(CMMICAV1_w ~ MedLipT1, data = f) 
      tapply(f$CMMICAV1_w, f$MedLipT1, sd, na.rm = TRUE)

      #Change in ICA differences by other cardiac meds   
      t.test(CMMICAV1_w ~ MedCT1, data = f) 
      
      #Change in ICA differences by diabetes meds   
      t.test(CMMICAV1_w ~ MedDT1, data = f) 
      tapply(f$CMMICAV1_w, f$MedDT1, sd, na.rm = TRUE)

      #Change in ICA differences by Race/Ethnicity
      anova_result <- aov(f$CMMICAV1_w ~ f$ref)
      summary(anova_result)

      #Change in ICA differences by Income
      anova_result <- aov(f$CMMICAV1_w ~ f$HIncome_3level)
      summary(anova_result) #p = 0.294

    
    #BIF/ICA
      #Correlations between social vigilance, change in BIF/ICA, age, BMI, resting SBP, resting DBP, time 1 BIF/ICA 
      h4covc <- subset (f, select = c(mEMASVQTV_boxcox, CCBIFICAV1_w, Age, BMIT1, SBPT1, DBPT1, MMBIFICAV1T1))
      print((corr.test(h4covc)), short=FALSE)
      corPlot(h4covc, stars = TRUE)
      
      #Change in BIF/ICA differences by Sex 
      t.test(CCBIFICAV1_w ~ Sex, data = f) 

      #Change in BIF/ICA differences by menopausal status 
      anova_result <- aov(f$CCBIFICAV1_w ~ f$MenopT1_V2)
      summary(anova_result) 
   
      #Change in BIF/ICA differences by blood pressure meds 
      t.test(CCBIFICAV1_w ~ MedBT1, data = f)

      #Change in BIF/ICA differences by lipid meds   
      t.test(CCBIFICAV1_w ~ MedLipT1, data = f) 

      #Change in BIF/ICA differences by other cardiac meds   
      t.test(CCBIFICAV1_w ~ MedCT1, data = f) 

      #Change in BIF/ICA differences by diabetes meds
      t.test(CCBIFICAV1_w ~ MedDT1, data = f)

      #Change in BIF/ICA differences by Race/Ethnicity
      anova_result <- aov(f$CCBIFICAV1_w ~ f$ref)
      summary(anova_result) 

      #Change in BIF/ICA differences by Income
      anova_result <- aov(f$CCBIFICAV1_w ~ f$HIncome_3level)
      summary(anova_result) 

      

#Descriptives----
      f %>%
        count(Sex) %>%
        mutate(percentage = n / sum(n) * 100)
      
            f %>%
        count(ref) %>%
        mutate(percentage = n / sum(n) * 100)
      psych::describe(f$Age)
      f %>%
        count(HIncome_f) %>%
        mutate(percentage = n / sum(n) * 100)
      summary(f$HIncome)
      f %>%
        count(HIncome_3level) %>%
        mutate(percentage = n / sum(n) * 100)
      f %>%
        count(MedLipT1) %>%
        mutate(percentage = n / sum(n) * 100)
      f %>%
        count(MedDT1) %>%
        mutate(percentage = n / sum(n) * 100)
      psych::describe(f$TotalCrime_1)
      psych::describe(f$TotalCrime_5)
      psych::describe(f$TotalCrime_10)
      psych::describe(f$socsaf_mean)
      psych::describe(f$mEMASVQTV)
      psych::describe(f$CMMFCCAV1_w)
      psych::describe(f$CMMBIFV1_w)
      psych::describe(f$CMMICAV1_w)
      psych::describe(f$CCBIFICAV1_w)
      
      psych::describe(f$SBPT1)
      psych::describe(f$DBPT1)
      psych::describe(f$BMIT1)
      psych::describe(f$MMFCCAV1T1)
      psych::describe(f$MMBIFV1T1)
      psych::describe(f$MMICAV1T1)
      psych::describe(f$MMBIFICAV1T1)
      
#MLM note----
  #Initially, we intended to test all environmental safety hypotheses via a multilevel analytic approach due to possible hierarchical effects
    #of geographic proximity. We tested mixed models testing linear regression and mediation with random effects clustered by zip code using the “lme4” 
    #and “mediation” packages for R (Bates et al., 2015; R Core Team, 2022; Tingley et al., 2014). 
    
 #The mixed models indicated almost no variance was explained by random effects clustering (level-1 residual variances ≤ 0.01) &
    #removing participants with incomplete data for random intercepts models resulted in insufficient variance to accurately model the clusters. 
      
 #Therefore, all environmental safety hypotheses were tested with standard linear models (see below) and there were no differences in observed results. 

  #This code and associated dataset do not contain those analyses because the geographic data utilized to assess for hierarchical effects were too sensitive 
   #to be made widely publicly available. For any inquiry into those analyses and the affiliated data used to replicate our reported results,
   #please reach out to the paper authors Riley M. O'Neill (rmoneill@arizona.edu) and/or Dr. John M. Ruiz (johnruiz@arizona.edu).
      
 
#Aim 1----     
  #Partial correlations 
      h1cov <- subset (f, select = c(TotalCrime_1z, TotalCrime_5z, TotalCrime_10z, socsaf_mean, mEMASVQTV_boxcox, CMMFCCAV1_w,CMMBIFV1_w, CMMICAV1_w, CCBIFICAV1_w, Age, Sex, ref, HIncome_3level, SBPT1, DBPT1, MMFCCAV1T1, MMBIFV1T1, MMICAV1T1, MMBIFICAV1T1, MedLipT1, BMIT1, MedDT1 ))
      h1par <- partial_correltable(data = h1cov, vars = c("TotalCrime_1z", "TotalCrime_5z", "TotalCrime_10z", "socsaf_mean", "mEMASVQTV_boxcox", "CMMFCCAV1_w","CMMBIFV1_w", "CMMICAV1_w","CCBIFICAV1_w"), partialvars = c("Age","Sex","ref","HIncome_3level","SBPT1","DBPT1","MMFCCAV1T1","MMBIFV1T1","MMICAV1T1","MMBIFICAV1T1","MedLipT1","BMIT1","MedDT1"))
      h1par
      names(h1par)
   
      
#Hypothesis 1. Environmental safety will be positively associated with social safety. 
  #Hypothesis 1: Total Crime 1 Mile 
       #Unadjusted model 
       lm1.1 <- lm(socsaf_mean ~TotalCrime_1z, data=g)
       summary(lm1.1) 
       #plot(lm1.1)
       confint(lm1.1) 

       #Adjusted model 
       lm1.2 <- lm(socsaf_mean ~TotalCrime_1z + Age + HIncome_3level + ref, data=g)
       summary(lm1.2) 
       #plot(lm1.2)
       confint(lm1.2) 
       
      
       #Post-Hoc Analyses
       #Interaction 1: TotalCrime 1 mile * Sex (controlling for age)
       lm1.int1 <- lm(socsaf_mean ~TotalCrime_1z*Sex + Age + HIncome_3level + ref, data=g)
       summary(lm1.int1) 

       #Interaction 2: TotalCrime 1 mile * Age
       lm1.int2 <- lm(socsaf_mean ~TotalCrime_1z*Age + HIncome_3level + ref, data=g)
       summary(lm1.int2) 

       #Interaction 3: TotalCrime 1 mile with Sex*Age
       lm1.int3 <- lm(socsaf_mean ~TotalCrime_1z+ Age*Sex + HIncome_3level + ref, data=g)
       summary(lm1.int3) 

       #Interaction 4: TotalCrime 1 mile*Sex*Age
       lm1.int4 <- lm(socsaf_mean ~TotalCrime_1z*Age*Sex + HIncome_3level + ref, data=g)
       summary(lm1.int4) 

  #Hypothesis 1: Total Crime 5 Mile
       #Unadjusted 
       lm2.1 <- lm(socsaf_mean ~TotalCrime_5z, data=g)
       summary(lm2.1)
       confint(lm2.1) 

       #Adjusted 
       lm2.2 <- lm(socsaf_mean ~TotalCrime_5z + Age + HIncome_3level+ ref, data=g)
       summary(lm2.2)
       confint(lm2.2) 
       
       #Interaction 1: TotalCrime 5 mile * Sex (controlling for age)
       lm2.int1 <- lm(socsaf_mean ~TotalCrime_5z*Sex + Age + HIncome_3level + ref, data=g)
       summary(lm2.int1) 

       #Interaction 2: TotalCrime 5 mile * Age
       lm2.int2 <- lm(socsaf_mean ~TotalCrime_5z*Age + HIncome_3level + ref, data=g)
       summary(lm2.int2) 

       #Interaction 3: TotalCrime 5 mile with Sex*Age
       lm2.int3 <- lm(socsaf_mean ~TotalCrime_5z+ Age*Sex + HIncome_3level + ref, data=g)
       summary(lm2.int3) 

       #Interaction 4: TotalCrime 5 mile*Sex*Age
       lm2.int4 <- lm(socsaf_mean ~TotalCrime_5z*Age*Sex + HIncome_3level + ref, data=g)
       summary(lm2.int4) 

     
  #Hypothesis 1: Total Crime 10 Mile 
       #Unadjusted
       lm3.1 <- lm(socsaf_mean ~TotalCrime_10z, data=g)
       summary(lm3.1)
       confint(lm3.1)
   
       #Adjusted 
       lm3.2 <- lm(socsaf_mean ~TotalCrime_10z + Age + HIncome_3level + ref, data=g)
       summary(lm3.2)
       confint(lm3.2) 
  
       #Interaction 1: TotalCrime 10 mile * Sex (controlling for age)
       lm3.int1 <- lm(socsaf_mean ~TotalCrime_10z*Sex + Age + HIncome_3level + ref, data=g)
       summary(lm3.int1) 
  
       #Interaction 2: TotalCrime 10 mile * Age
       lm3.int2 <- lm(socsaf_mean ~TotalCrime_10z*Age + HIncome_3level + ref, data=g)
       summary(lm3.int2) 
  
       #Interaction 3: TotalCrime 10 mile with Sex*Age
       lm3.int3 <- lm(socsaf_mean ~TotalCrime_10z+ Age*Sex + HIncome_3level + ref, data=g)
       summary(lm3.int3) 
  
       #Interaction 4: TotalCrime 10 mile*Sex*Age
       lm3.int4 <- lm(socsaf_mean ~TotalCrime_10z*Age*Sex + HIncome_3level + ref, data=g)
       summary(lm3.int4) 

     
     
#Hypothesis 2: Environmental safety will be negatively associated with social vigilance behavior in daily life. 

  #Hypothesis 2: Total Crime 1 mile
     #Unadjusted 
     lm4.1 <- lm(mEMASVQTV_boxcox ~TotalCrime_1z, data=g)
     summary(lm4.1)
     confint(lm4.1) 
     #plot(lm4.1)
     hist(resid(lm4.1)) 
    
     #Adjusted 
     lm4.2 <- lm(mEMASVQTV_boxcox ~TotalCrime_1z + Age + HIncome_3level + Sex + ref, data=g)
     summary(lm4.2)
     confint(lm4.2) 
     #plot(lm4.2)  
     hist(resid(lm4.2)) 
    
     #Interaction 1: TotalCrime 1mile*Sex (controlling for Age) 
     lm4.int1 <- lm(mEMASVQTV_boxcox ~TotalCrime_1z*Sex + Age + HIncome_3level + ref, data=g)
     summary(lm4.int1)

     #Interaction 2: TotalCrime 1mile*Age 
     lm4.int2 <- lm(mEMASVQTV_boxcox ~TotalCrime_1z*Age + Sex + HIncome_3level + ref, data=g)
     summary(lm4.int2)
  
     #Interaction 3: TotalCrime 1mile + Sex*Age 
     lm4.int3 <- lm(mEMASVQTV_boxcox ~TotalCrime_1z +Sex*Age + HIncome_3level + ref, data=g)
     summary(lm4.int3)
  
     #Interaction 4: TotalCrime 1mile*Sex*Age 
     lm4.int4 <- lm(mEMASVQTV_boxcox ~TotalCrime_1z*Sex*Age + HIncome_3level + ref, data=g)
     summary(lm4.int4)

     
  #Hypothesis 2: Total Crime 5 mile
     #Unadjusted 
     lm5.1 <- lm(mEMASVQTV_boxcox ~TotalCrime_5z, data=g)
     summary(lm5.1)
     confint(lm5.1) 
     #plot(lm5.1) 
     hist(resid(lm5.1))
  
     
     #Adjusted 
     lm5.2 <- lm(mEMASVQTV_boxcox ~TotalCrime_5z + Age + HIncome_3level + Sex + ref, data=g)
     summary(lm5.2)
     confint(lm5.2) 
     #plot(lm5.2) 
     hist(resid(lm5.2))
   
     #Interaction 1: TotalCrime 5mile*Sex (controlling for Age) 
     lm5.int1 <- lm(mEMASVQTV_boxcox ~TotalCrime_5z*Sex + Age + HIncome_3level + ref, data=g)
     summary(lm5.int1)
     #significant negative effects of sex and age on SV 
     
     #Interaction 2: TotalCrime 5mile*Age 
     lm5.int2 <- lm(mEMASVQTV_boxcox ~TotalCrime_5z*Age + Sex + HIncome_3level + ref, data=g)
     summary(lm5.int2)
     #significant negative effects of age and sex on sV 
     
     #Interaction 3: TotalCrime 5mile + Sex*Age 
     lm5.int3 <- lm(mEMASVQTV_boxcox ~TotalCrime_5z +Sex*Age + HIncome_3level + ref, data=g)
     summary(lm5.int3)
     #no significant effects on SV 
     
     #Interaction 4: TotalCrime 5mile*Sex*Age 
     lm5.int4 <- lm(mEMASVQTV_boxcox ~TotalCrime_5z*Sex*Age + HIncome_3level + ref, data=g)
     summary(lm5.int4)
     #no significant effects on SV 
     #Total Crime 5 mile * Age approaches significance (p = .097)
     
     
     
  #Hypothesis 2: Total Crime 10 mile 
     #Unadjusted 
     lm6.1 <- lm(mEMASVQTV_boxcox ~TotalCrime_10lz, data=g)
     summary(lm6.1)
     confint(lm6.1) 
     #plot(lm6.1)
     hist(resid(lm6.1)) 
    
     #Adjusted 
     lm6.2 <- lm(mEMASVQTV_boxcox ~TotalCrime_10z + Age + HIncome_3level + Sex + ref, data=g)
     summary(lm6.2)
     confint(lm6.2) 
     #plot(lm6.2) 
     hist(resid(lm6.2)) 
     
     #Interaction 1: TotalCrime 10mile*Sex (controlling for Age) 
     lm6.int1 <- lm(mEMASVQTV_boxcox ~TotalCrime_10z*Sex + Age + HIncome_3level + ref, data=g)
     summary(lm6.int1)
     #Significant negative effect of sex on SV (<.001)
     #negative effect of Age on SV approached significance (p = .058)
     
     #Interaction 2: TotalCrime 10mile*Age 
     lm6.int2 <- lm(mEMASVQTV_boxcox ~TotalCrime_10z*Age + Sex + HIncome_3level + ref, data=g)
     summary(lm6.int2)
     #Significant negative effect of sex on SV (<.001)
     #total crime 10 mile * Age effect on SV approached significance (p = .060)
     
     #Interaction 3: TotalCrime 10mile + Sex*Age 
     lm6.int3 <- lm(mEMASVQTV_boxcox ~TotalCrime_10z +Sex*Age + HIncome_3level + ref, data=g)
     summary(lm6.int3)
     #no significant effects
     
     #Interaction 4: TotalCrime 10mile*Sex*Age 
     lm6.int4 <- lm(mEMASVQTV_boxcox ~TotalCrime_10z*Sex*Age + HIncome_3level + ref, data=g)
     summary(lm6.int4)
     #significant positive effect of total crime 10 mile*Age (B = 0.00, p = .042)
     
     
     

#Hypothesis 3: Social safety will be negatively associated with social vigilance behavior in daily life. 
      
     #Unadjusted
      lm7.1 <- lm(mEMASVQTV_boxcox ~ socsaf_mean, data=f)
      Anova(lm7.1, type="III") 
      summary(lm7.1) 
        
     #Adjusted
      lm7.2 <- lm(mEMASVQTV_boxcox ~ socsaf_mean + Age + Sex, data=f)
      Anova(lm7.2, type="III")
      summary(lm7.2) 
      confint(lm7.2)

      #APA-format adjusted effect plot, APA -style 
      pred <- ggpredict(lm7.2, terms = "socsaf_mean")
              
      ggplot(pred, aes(x = x, y = predicted)) +
              geom_line(linewidth = 1.2) +  # Use thicker line for clarity
              geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.05, fill = "gray30") +  # APA-style confidence interval
              labs(
                x = "Social Safety",  # Axis labels in sentence case
                y = "Social Vigilance (box-cox transformed)", 
                title = " "
              ) +
              theme_minimal() +  # Minimal theme is closest to APA style
              theme(
                # Customize text for APA style
                text = element_text(family = "serif", size = 12),
                plot.title = element_text(face = "bold", hjust = 0.5, size = 14),  # Center title, bolded, larger font
                axis.title = element_text(size = 12),  # Axis titles slightly larger
                axis.text = element_text(size = 10),  # Axis text smaller but readable
                panel.grid = element_blank(),  # Remove all grid lines
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border around plot
                plot.margin = margin(15, 15, 15, 15)  # Add some padding
              )
          

      #Social Safety x Social Vigilance (with sex*social safety interaction)
          mint1 <- lm(mEMASVQTV_boxcox ~ socsaf_mean*Sex, data=f)
          Anova(mint1, type="III")
          summary(mint1)
          
      #Social Safety x Social Vigilance (with sex*social safety interaction + controlling for Age)
          mint2 <- lm(mEMASVQTV_boxcox ~ socsaf_mean*Sex + Age, data=f)
          Anova(mint2, type="III") 
          summary(mint2) 
      
          
      #Social Safety x Social Vigilance (with age*social safety interaction)
          mint3 <- lm(mEMASVQTV_boxcox ~ socsaf_mean*Age, data=f)
          Anova(mint3, type="III") 
          summary(mint3)  

      #Social Safety x Social Vigilance (with age*social safety interaction + controlling for sex)
          mint4 <- lm(mEMASVQTV_boxcox ~ socsaf_mean*Age + Sex, data=f)
          Anova(mint4, type="III") 
          summary(mint4) 
          
      #Social Safety x Social Vigilance (with age*sex interaction)
          mint5 <- lm(mEMASVQTV_boxcox ~ socsaf_mean + Age*Sex, data=f)
          Anova(mint5, type="III") 
          summary(mint5) 
          
      #Social Safety x Social Vigilance (with socialsafety*age*sex interaction)
          mint6 <- lm(mEMASVQTV_boxcox ~ socsaf_mean*Age*Sex, data=f)
          Anova(mint6, type="III") 
          summary(mint6)    
          
#Hypothesis 4: Social vigilance behavior at time 1 will be associated with greater increases in cIMT at 2-year follow-up. 
    #Hypothesis 4: FCCA
          #Unadjusted
          lm9.1 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox, data=f)
          Anova(lm9.1, type="III")
          summary(lm9.1) 
         
          #Adjusted
          lm9.2 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox + Age + Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
          Anova(lm9.2, type="III") 
          summary(lm9.2) 
          confint(lm9.2)

       `` #Social Vigilance x FCCA V1 change (with SV*Sex interaction + other adjustments without Age)
                mint1 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox* Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
                Anova(mint1, type="III") 
                summary(mint1) 

          #Social Vigilance x FCCA V1 change (with SV*Sex interaction + other adjustments with Age)
              mint2 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox* Sex + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
              Anova(mint2, type="III") 
              summary(mint2) 
        
          #Social Vigilance x FCCA V1 change (with SV*Age + other adjustments without Sex)
              mint3 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox*Age + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
              Anova(mint3, type="III") 
              summary(mint3) 
  
          #Social Vigilance x FCCA V1 change (with SV*Age + other adjustments with Sex)
              mint4 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox*Age + Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
              Anova(mint4, type="III")  
              summary(mint4) 
    
          #Social Vigilance x FCCA V1 change (with SV + Age*Sex)
              mint5 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox + Age*Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
              Anova(mint5, type="III")   
              summary(mint5) 
   
          #Social Vigilance x FCCA V1 change (with SV*Age*Sex)
              mint6 <- lm(CMMFCCAV1_w ~ mEMASVQTV_boxcox*Age*Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + HIncome_3level + MedLipT1  , data=f)
              Anova(mint6, type="III")  
              summary(mint6)      
    
    #Hypothesis 4: BIF
          #Unadjusted
          lm10.1 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox, data=f)
          Anova(lm10.1, type="III") 
          summary(lm10.1) 
          
          #Adjusted
          lm10.2 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox + Age + Sex + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(lm10.2, type="III") 
          summary(lm10.2) 
          confint(lm10.2)
          
         #Social Vigilance x BIF V1 change (SV*sex + adjustments without age)
          mint1 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox*Sex + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(mint1, type="III") 
          summary(mint1) 
            
        
          #Social Vigilance x BIF V1 change (SV*sex + adjustments with age)
          mint2 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox*Sex + Age + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(mint2, type="III") 
          summary(mint2) 
          
          #Social Vigilance x BIF V1 change (SV*age + adjustments without sex)
          mint3 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox*Age + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(mint3, type="III") 
          summary(mint3) 
         
          #Social Vigilance x BIF V1 change (SV*age + adjustments with sex)
          mint4 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox*Age + Sex + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(mint4, type="III") 
          summary(mint4) 
          
          #Social Vigilance x BIF V1 change (SV + age*sex + adjustments)
          mint5 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox + Age*Sex + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(mint5, type="III")
          summary(mint5)  
          
          #Social Vigilance x BIF V1 change (SV*age*sex + adjustments)
          mint6 <- lm(CMMBIFV1_w ~ mEMASVQTV_boxcox*Age*Sex + BMIT1 + MMBIFV1T1 + HIncome_3level  , data=f)
          Anova(mint6, type="III") 
          summary(mint6)  
          confint(mint6)
          
          #Follow-up probing for significant 3-way interaction effect: 
              # Probing the slope of mEMASVQTV_boxcox at different ages for each sex
              emtrends(mint6, ~ Sex | Age, var = "mEMASVQTV_boxcox", at = list(Age = c(20, 40, 60)))
          
              emmip(mint6, Sex ~ mEMASVQTV_boxcox | Age, at = list(Age = c(20, 40, 60)), cov.reduce = FALSE)
              

    #Hypothesis 4: ICA
          #Unadjusted
          lm11.1 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox, data=f)
          Anova(lm11.1 , type="III") 
          summary(lm11.1 ) 
          
          
          #Adjusted 
          lm11.2 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox + Age + Sex + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(lm11.2, type="III") 
          summary(lm11.2) 
          confint(lm11.2)
        
          
          #Social Vigilance x ICA V1 change (SV*sex + adjustments without age)
          mint1 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox*Sex + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(mint1, type="III") 
          summary(mint1) 
          
  
          #Social Vigilance x ICA V1 change (SV*sex + adjustments with age)
          mint2 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox*Sex + Age + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(mint2, type="III") 
          summary(mint2) 
          
          #Social Vigilance x ICA V1 change (SV*age + adjustments without sex)
          mint3 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox*Age + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(mint3, type="III") 
          summary(mint3) 
          
          #Social Vigilance x ICA V1 change (SV*age + adjustments with sex)
          mint4 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox*Age + Sex + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(mint4, type="III") 
          summary(mint4) 
          
          #Social Vigilance x ICA V1 change (SV + age*sex + other adjustments)
          mint5 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox + Age*Sex + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(mint5, type="III") 
          summary(mint5) 
          
          #Social Vigilance x ICA V1 change (SV*age*sex + other adjustments)
          mint6 <- lm(CMMICAV1_w ~ mEMASVQTV_boxcox*Age*Sex + BMIT1 + MMICAV1T1 + HIncome_3level + MedLipT1 + MedDT1  , data=f)
          Anova(mint6, type="III") 
          summary(mint6) 
      
    #Hypothesis 4: BIF/ICA
          #Unadjusted
          lm12.1 <- lm(CCBIFICAV1_w ~ mEMASVQTV_boxcox, data=f)
          Anova(lm12.1, type="III") 
          summary(lm12.1) 
         
          #Adjusted
          lm12.2 <- lm(CCBIFICAV1_w ~ mEMASVQTV_boxcox + Age + Sex + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1, data=f)
          Anova(lm12.2, type="III") 
          summary(lm12.2) 
          confint(lm12.2)
          
          
         #Social Vigilance x BIF/ICA V1 change (SV*sex + adjustments without age)
          mint1 <- lm(CCBIFICAV1 ~ mEMASVQTV_boxcox*Sex + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1  , data=f)
          Anova(mint1, type="III") 
          summary(mint1) 
          
          #Social Vigilance x BIF/ICA V1 change (SV*sex + adjustments with age)
          mint2 <- lm(CCBIFICAV1 ~ mEMASVQTV_boxcox*Sex + Age + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1  , data=f)
          Anova(mint2, type="III") 
          summary(mint2) 
          
          
          #Social Vigilance x BIF/ICA V1 change (SV*age + adjustments without sex)
          mint3 <- lm(CCBIFICAV1 ~ mEMASVQTV_boxcox*Age + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1  , data=f)
          Anova(mint3, type="III") 
          summary(mint3) 
          
         
          #Social Vigilance x BIF/ICA V1 change (SV*age + adjustments with sex)
          mint4 <- lm(CCBIFICAV1 ~ mEMASVQTV_boxcox*Age + Sex + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1  , data=f)
          Anova(mint4, type="III") 
          summary(mint4) 
          
          #Social Vigilance x BIF/ICA V1 change (SV + age*sex + adjustments)
          mint5 <- lm(CCBIFICAV1 ~ mEMASVQTV_boxcox + Age*Sex + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1  , data=f)
          Anova(mint5, type="III") 
          summary(mint5) 
          
          #Social Vigilance x BIF/ICA V1 change (SV*age*sex + adjustments)
          mint6 <- lm(CCBIFICAV1 ~ mEMASVQTV_boxcox*Age*Sex + DBPT1 + MMBIFICAV1T1 + HIncome_3level + BMIT1  , data=f)
          Anova(mint6, type="III") 
          summary(mint6) 
          
   
#Aim 2 ----
 #Hypothesis 5. Social vigilance will mediate or partially mediate any observed effect of environmental safety on 2-year change in cIMT. 

    #Hypothesis 5: 1-mile risk & FCCA
      #Unadjusted 
    
      model <- ' #Direct effect between X and Y
                CMMFCCAV1_w ~ C*TotalCrime_1z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z
    
                #Path B
                CMMFCCAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B

                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
      #Adjusted linear regression for mediation model 
      model <- ' #Direct effect between X and Y
                CMMFCCAV1_w ~ C*TotalCrime_1z + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
    
                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)

      
    #Interaction 1: Crime*Sex 
        model_fcca_1mi_int1 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C1*TotalCrime_1z + C2*Sex_num + C3*Crime1z_Sex + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Sex_num + A3*Crime1z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A3*B
                  
                  #Total Effect 
                  TotalEffect := C3 + (A3*B)
                '
        
        mediationfit_fcca_1mi_int1 <- sem(model_fcca_1mi_int1, data = g) 
        #summary(mediationfit_1mi_int1)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_1mi_int1 <- sem(model_fcca_1mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
    
        summary(mediationfit_fcca_1mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
    
     #Interaction 2: Crime*Age
        model_fcca_1mi_int2 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Crime1z_Age + Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Crime1z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A3*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C3 + (A3*B)
                '
        
        mediationfit_fcca_1mi_int2 <- sem(model_fcca_1mi_int2, data = g) 
        #summary(mediationfit_fcca_1mi_int2)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_1mi_int2 <- sem(model_fcca_1mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
        
        summary(mediationfit_fcca_1mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
    
    #Interaction 3: Crime + Sex*Age
        model_fcca_1mi_int3 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C*TotalCrime_1z + Age + Sex_num + Age_Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex_num + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
        
        mediationfit_fcca_1mi_int3 <- sem(model_fcca_1mi_int3, data = g) 
        #summary(mediationfit_fcca_1mi_int3)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_1mi_int3 <- sem(model_fcca_1mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
        
        summary(mediationfit_fcca_1mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
    
    #Interaction 4: Crime*Sex*Age
        model_fcca_1mi_int4 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Sex_num + C4*Crime1z_Sex_Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Sex_num + A4*Crime1z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A4*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C4 + (A4*B)
                '
        
        mediationfit_fcca_1mi_int4 <- sem(model_fcca_1mi_int4, data = g) 
        #summary(mediationfit_1mi_int4)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_1mi_int4 <- sem(model_fcca_1mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
        
        summary(mediationfit_fcca_1mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
    
 
  
  
  #Hypothesis 5: 1-mile risk & BIF
            
      #Unadjusted
      model <- ' #Direct effect between X and Y
                CMMBIFV1_w ~ C*TotalCrime_1z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z
    
                #Path B
                CMMBIFV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
      #Adjusted
      model <- ' #Direct effect between X and Y
                CMMBIFV1_w ~ C*TotalCrime_1z + Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMBIFV1T1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
    #Interaction 1: Crime*Sex 
          model_bif_1mi_int1<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_1z + C2*Sex_num + C3*Crime1z_Sex + Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Sex_num + A3*Crime1z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A3*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C3 + (A3*B)
                                '
                        
          mediationfit_bif_1mi_int1 <- sem(model_bif_1mi_int1, data = g) 
          #summary(mediationfit_bif_1mi_int1)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_bif_1mi_int1 <- sem(model_bif_1mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_bif_1mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
      #Interaction 2: Crime*Age
          model_bif_1mi_int2<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Crime1z_Age + Sex_num + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Crime1z_Age + Sex_num + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A3*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C3 + (A3*B)
                                '
          
          mediationfit_bif_1mi_int2 <- sem(model_bif_1mi_int2, data = g) 
          #summary(mediationfit_bif_1mi_int2)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_bif_1mi_int2 <- sem(model_bif_1mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_bif_1mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
          
     
      #Interaction 3: Crime + Sex*Age
          model_bif_1mi_int3<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C*TotalCrime_1z + Age + Sex_num + + Age_Sex + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex_num + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C + (A*B)
                                '
          
          mediationfit_bif_1mi_int3 <- sem(model_bif_1mi_int3, data = g) 
          #summary(mediationfit_bif_1mi_int3)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_bif_1mi_int3 <- sem(model_bif_1mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_bif_1mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
      #Interaction 4: Crime*Sex*Age
          model_bif_1mi_int4<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Sex_num + C4*Crime1z_Sex_Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Sex_num + A4*Crime1z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A4*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C4 + (A4*B)
                                '
          
          mediationfit_bif_1mi_int4 <- sem(model_bif_1mi_int4, data = g) 
          #summary(mediationfit_bif_1mi_int4)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_bif_1mi_int4 <- sem(model_bif_1mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_bif_1mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
     
      
  #Hypothesis 5: 1-mile risk & ICA  
      #Unadjusted
      model <- ' #Direct effect between X and Y
                CMMICAV1_w ~ C*TotalCrime_1z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z
    
                #Path B
                CMMICAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
     
 
    #Adjusted
      model <- ' #Direct effect between X and Y
                CMMICAV1_w ~ C*TotalCrime_1z + Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
    #Interaction 1: Crime*Sex
          model_ica_1mi_int1 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_1z + C2*Sex_num + C3*Crime1z_Sex + Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Sex_num + A3*Crime1z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A3*B
    
                                      #Total Effect 
                                      TotalEffect := C3 + (A3*B)
                                    '
                            
          mediationfit_ica_1mi_int1 <- sem(model_ica_1mi_int1, data = g) 
          #summary(mediationfit_ica_1mi_int1)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_ica_1mi_int1 <- sem(model_ica_1mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_ica_1mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
          
      #Interaction 2: Crime*Age
          model_ica_1mi_int2 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Crime1z_Age + Sex + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Crime1z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A3*B
    
                                      #Total Effect 
                                      TotalEffect := C3 + (A3*B)
                                    '
          
          mediationfit_ica_1mi_int2 <- sem(model_ica_1mi_int2, data = g) 
          #summary(mediationfit_ica_1mi_int2)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_ica_1mi_int2 <- sem(model_ica_1mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_ica_1mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
     #Interaction 3: Crime + Age*Sex
          model_ica_1mi_int3 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C*TotalCrime_1z + Age + Sex + Age_Sex + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A*B
    
                                      #Total Effect 
                                      TotalEffect := C + (A*B)
                                    '
          
          mediationfit_ica_1mi_int3 <- sem(model_ica_1mi_int3, data = g) 
          #summary(mediationfit_ica_1mi_int3)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_ica_1mi_int3 <- sem(model_ica_1mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_ica_1mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
          
        
      #Interaction 4: Crime*Sex*Age
          model_ica_1mi_int4 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Sex_num + C4*Crime1z_Sex_Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Sex_num + A4*Crime1z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A4*B
    
                                      #Total Effect 
                                      TotalEffect := C4 + (A4*B)
                                    '
          
          mediationfit_ica_1mi_int4 <- sem(model_ica_1mi_int4, data = g) 
          #summary(mediationfit_ica_1mi_int4)
          
          # Modify the sem function to perform bootstrapping
          mediationfit_ica_1mi_int4 <- sem(model_ica_1mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
          
          # Display summary with confidence intervals
          summary(mediationfit_ica_1mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
  #Hypothesis 5: 1-mile risk & BIF/ICA  
    #Unadjusted 
        model <- ' #Direct effect between X and Y
                CCBIFICAV1_w ~ C*TotalCrime_1z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z
    
                #Path B
                CCBIFICAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
    #Adjusted
      model <- ' #Direct effect between X and Y
                CCBIFICAV1_w ~ C*TotalCrime_1z + Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMBIFICAV1T1 + DBPT1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
 
    #Interaction 1: Crime*Sex 
              model_bifica_1mi_int1 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_1z + C2*Sex_num + C3*Crime1z_Sex + Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Sex_num + A3*Crime1z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A3*B
                            #Creating new parameter testing for the indirect effects
                
                            #Total Effect 
                            TotalEffect := C3 + (A3*B)
                          '
                  
                mediationfit_bifica_1mi_int1 <- sem(model_bifica_1mi_int1, data = g) 
                #summary(mediationfit_bifica_1mi_int1)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_bifica_1mi_int1 <- sem(model_bifica_1mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_bifica_1mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
            
      #Interaction 2: Crime*Age 
                model_bifica_1mi_int2 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Crime1z_Age + Sex + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A1*TotalCrime_1z + A2*Age + A3*Crime1z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A3*B
                            #Creating new parameter testing for the indirect effects
                
                            #Total Effect 
                            TotalEffect := C3 + (A3*B)
                          '
                
                mediationfit_bifica_1mi_int2 <- sem(model_bifica_1mi_int2, data = g) 
                #summary(mediationfit_bifica_1mi_int2)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_bifica_1mi_int2 <- sem(model_bifica_1mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_bifica_1mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
     
          #Plotting simple slopes of Crime*Age interaction 
                # Coefficients from output
                b1 <- 0.086   # main effect of TtlCrm_1z
                b3 <- -0.002  # interaction effect: TtlCrm_1z * Age
                
                # Setting age categories
                ages <- c(20, 40, 60)
                
                # Calculating simple slopes
                simple_slopes <- b1 + b3 * ages
                
                df_slopes <- data.frame(
                  Age = ages,
                  SimpleSlope = simple_slopes
                )
                
                #Plotting simple slopes 
                ggplot(df_slopes, aes(x = Age, y = SimpleSlope)) +
                  geom_line(color = "#0072B2", size = 1.2) +
                  geom_point(size = 3, color = "#D55E00") +
                  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
                  labs(
                    title = "Simple Slopes of Total Crime on BIF/ICA Change by Age",
                    x = "Age",
                    y = "Effect of Total Crime (TtlCrm_1z)"
                  ) +
                  theme_minimal(base_size = 14)
                
                
      #Interaction 3: Crime + Sex*Age 
                model_bifica_1mi_int3 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C*TotalCrime_1z + Age + Sex + Age_Sex + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A*TotalCrime_1z + Age + Sex + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A*B
                            
                            #Total Effect 
                            TotalEffect := C + (A*B)
                          '
                
                mediationfit_bifica_1mi_int3 <- sem(model_bifica_1mi_int3, data = g) 
                #summary(mediationfit_bifica_1mi_int3)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_bifica_1mi_int3 <- sem(model_bifica_1mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_bifica_1mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
    #Interaction 4: Crime*Sex*Age 
                model_bifica_1mi_int4 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_1z + C2*Age + C3*Sex_num + C4*Crime1z_Sex_Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~  A1*TotalCrime_1z + A2*Age + A3*Sex_num + A4*Crime1z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A4*B
                            
                            #Total Effect 
                            TotalEffect := C4 + (A4*B)
                          '
                
                mediationfit_bifica_1mi_int4 <- sem(model_bifica_1mi_int4, data = g) 
                #summary(mediationfit_bifica_1mi_int4)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_bifica_1mi_int4 <- sem(model_bifica_1mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_bifica_1mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
        
  #Hypothesis 5: 5-mile risk & FCCA 
      #Unadjusted  
        model <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C*TotalCrime_5z

                  #Path A 
                  mEMASVQTV_boxcox ~ A*TotalCrime_5z
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
        mediationfit <- sem(model, data = g) 
        #summary(mediationfit)
        
        # Modify the sem function to perform bootstrapping
        mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
        
        # Display summary with confidence intervals
        summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
    
      #Adjusted 
        model <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C*TotalCrime_5z + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 

                  #Path A 
                  mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
        
        mediationfit <- sem(model, data = g) 
        #summary(mediationfit)
        
        # Modify the sem function to perform bootstrapping
        mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
        
        # Display summary with confidence intervals
        summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
   
    #Interaction 1: Crime*Sex 
        model_fcca_5mi_int1 <- ' #Direct effect between X and Y
                    CMMFCCAV1_w ~ C1*TotalCrime_5z + C2*Sex_num + C3*Crime5z_Sex + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                    
                    #Path A 
                    mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Sex_num + A3*Crime5z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
        
                    #Path B
                    CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
              
                    #Indirect Effect
                    IndirectEffect := A3*B
                    #Creating new parameter testing for the indirect effects
        
                    #Total Effect 
                    TotalEffect := C3 + (A3*B)
                  '
        
        mediationfit_fcca_5mi_int1 <- sem(model_fcca_5mi_int1, data = g) 
        #summary(mediationfit_5mi_int1)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_5mi_int1 <- sem(model_fcca_5mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
        
        summary(mediationfit_fcca_5mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
      #Interaction 2: Crime*Age
        model_fcca_5mi_int2 <- ' #Direct effect between X and Y
                    CMMFCCAV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Crime5z_Age + Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                    
                    #Path A 
                    mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Crime5z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
        
                    #Path B
                    CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
              
                    #Indirect Effect
                    IndirectEffect := A3*B
                    #Creating new parameter testing for the indirect effects
        
                    #Total Effect 
                    TotalEffect := C3 + (A3*B)
                  '
        
        mediationfit_fcca_5mi_int2 <- sem(model_fcca_5mi_int2, data = g) 
        #summary(mediationfit_fcca_5mi_int2)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_5mi_int2 <- sem(model_fcca_5mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
        
        summary(mediationfit_fcca_5mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
       
    #Interaction 3: Crime + Sex*Age
        model_fcca_5mi_int3 <- ' #Direct effect between X and Y
                    CMMFCCAV1_w ~ C*TotalCrime_5z + Age + Sex_num + Age_Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                    
                    #Path A 
                    mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex_num + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
        
                    #Path B
                    CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
              
                    #Indirect Effect
                    IndirectEffect := A*B
                    #Creating new parameter testing for the indirect effects
        
                    #Total Effect 
                    TotalEffect := C + (A*B)
                  '
        
        mediationfit_fcca_5mi_int3 <- sem(model_fcca_5mi_int3, data = g) 
        #summary(mediationfit_fcca_5mi_int3)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_5mi_int3 <- sem(model_fcca_5mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
        
        summary(mediationfit_fcca_5mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
        
   
    #Interaction 4: Crime*Sex*Age
        model_fcca_5mi_int4 <- ' #Direct effect between X and Y
                    CMMFCCAV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Sex_num + C4*Crime5z_Sex_Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                    
                    #Path A 
                    mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Sex_num + A4*Crime5z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
        
                    #Path B
                    CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
              
                    #Indirect Effect
                    IndirectEffect := A4*B
                    #Creating new parameter testing for the indirect effects
        
                    #Total Effect 
                    TotalEffect := C4 + (A4*B)
                  '
        
        mediationfit_fcca_5mi_int4 <- sem(model_fcca_5mi_int4, data = g) 
        #summary(mediationfit_5mi_int4)
        
        # Modify the sem function to perform bootstrapping
        mediationfit_fcca_5mi_int4 <- sem(model_fcca_5mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
        
     
  #Hypothesis 5: 5-mile risk & BIF
      #Unadjusted
      model <- ' #Direct effect between X and Y
                CMMBIFV1_w ~ C*TotalCrime_5z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_5z
    
                #Path B
                CMMBIFV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
   #Adjusted
      model <- ' #Direct effect between X and Y
                CMMBIFV1_w ~ C*TotalCrime_5z + Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMBIFV1T1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
 
    #Interaction 1: Crime*Sex 
      model_bif_5mi_int1<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_5z + C2*Sex_num + C3*Crime5z_Sex + Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Sex_num + A3*Crime5z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A3*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C3 + (A3*B)
                                '
      
      mediationfit_bif_5mi_int1 <- sem(model_bif_5mi_int1, data = g) 
      #summary(mediationfit_bif_5mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_5mi_int1 <- sem(model_bif_5mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_5mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE) 
      
  
   #Interaction 2: Crime*Age
      model_bif_5mi_int2<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Crime5z_Age + Sex_num + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Crime5z_Age + Sex_num + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A3*B
                                  
                      
                                  #Total Effect 
                                  TotalEffect := C3 + (A3*B)
                                '
      
      mediationfit_bif_5mi_int2 <- sem(model_bif_5mi_int2, data = g) 
      #summary(mediationfit_bif_5mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_5mi_int2 <- sem(model_bif_5mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_5mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)   
      
    #Interaction 3: Crime + Sex*Age
      model_bif_5mi_int3<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C*TotalCrime_5z + Age + Sex_num + + Age_Sex + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex_num + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C + (A*B)
                                '
      
      mediationfit_bif_5mi_int3 <- sem(model_bif_5mi_int3, data = g) 
      #summary(mediationfit_bif_5mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_5mi_int3 <- sem(model_bif_5mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_5mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
   #Interaction 4: Crime*Sex*Age
      model_bif_5mi_int4 <- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Sex_num + C4*Crime5z_Sex_Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Sex_num + A4*Crime5z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A4*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C4 + (A4*B)
                                '
      
      mediationfit_bif_5mi_int4 <- sem(model_bif_5mi_int4, data = g) 
      #summary(mediationfit_bif_5mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_5mi_int4 <- sem(model_bif_5mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_5mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
  
  #Hypothesis 5: 5-mile risk & ICA  
      #Unadjusted
      model <- ' #Direct effect between X and Y
                CMMICAV1_w ~ C*TotalCrime_5z
                #C* means I am giving the parameter the name C
    
                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_5z
    
                #Path B
                CMMICAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)

      
      
      #Adjusted
      model <- ' #Direct effect between X and Y
                CMMICAV1_w ~ C*TotalCrime_5z + Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)

      
    #Interaction 1: Crime*Sex
      model_ica_5mi_int1 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_5z + C2*Sex_num + C3*Crime5z_Sex + Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Sex_num + A3*Crime5z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A3*B
    
                                      #Total Effect 
                                      TotalEffect := C3 + (A3*B)
                                    '
      
      mediationfit_ica_5mi_int1 <- sem(model_ica_5mi_int1, data = g) 
      #summary(mediationfit_ica_5mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_5mi_int1 <- sem(model_ica_5mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_5mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
  #Interaction 2: Crime*Age
      model_ica_5mi_int2 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Crime5z_Age + Sex + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Crime5z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A3*B
    
                                      #Total Effect 
                                      TotalEffect := C3 + (A3*B)
                                    '
      
      mediationfit_ica_5mi_int2 <- sem(model_ica_5mi_int2, data = g) 
      #summary(mediationfit_ica_5mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_5mi_int2 <- sem(model_ica_5mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_5mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
   #Interaction 3: Crime + Age*Sex
      model_ica_5mi_int3 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C*TotalCrime_5z + Age + Sex + Age_Sex + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A*B
    
                                      #Total Effect 
                                      TotalEffect := C + (A*B)
                                    '
      
      mediationfit_ica_5mi_int3 <- sem(model_ica_5mi_int3, data = g) 
      #summary(mediationfit_ica_5mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_5mi_int3 <- sem(model_ica_5mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_5mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
      
    #Interaction 4: Crime*Sex*Age
      model_ica_5mi_int4 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Sex_num + C4*Crime5z_Sex_Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Sex_num + A4*Crime5z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A4*B
    
                                      #Total Effect 
                                      TotalEffect := C4 + (A4*B)
                                    '
      
      mediationfit_ica_5mi_int4 <- sem(model_ica_5mi_int4, data = g) 
      #summary(mediationfit_ica_5mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_5mi_int4 <- sem(model_ica_5mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_5mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
  #Hypothesis 5: 5-mile risk & BIF/ICA 
      #Unadjusted 
      model <- ' #Direct effect between X and Y
                CCBIFICAV1_w ~ C*TotalCrime_5z
                #C* means I am giving the parameter the name C
    
                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_5z
    
                #Path B
                CCBIFICAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
     
      #Adjusted
      model <- ' #Direct effect between X and Y
                CCBIFICAV1_w ~ C*TotalCrime_5z + Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMBIFICAV1T1 + DBPT1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
      #Interaction 1: Crime*Sex 
      model_bifica_5mi_int1 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_5z + C2*Sex_num + C3*Crime5z_Sex + Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Sex_num + A3*Crime5z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A3*B
                            #Creating new parameter testing for the indirect effects
                
                            #Total Effect 
                            TotalEffect := C3 + (A3*B)
                          '
      
      mediationfit_bifica_5mi_int1 <- sem(model_bifica_5mi_int1, data = g) 
      #summary(mediationfit_bifica_5mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_5mi_int1 <- sem(model_bifica_5mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_5mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE) 
   
      
    #Interaction 2: Crime*Age 
      model_bifica_5mi_int2 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Crime5z_Age + Sex + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A1*TotalCrime_5z + A2*Age + A3*Crime5z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A3*B
                            #Creating new parameter testing for the indirect effects
                
                            #Total Effect 
                            TotalEffect := C3 + (A3*B)
                          '
      
      mediationfit_bifica_5mi_int2 <- sem(model_bifica_5mi_int2, data = g) 
      #summary(mediationfit_bifica_5mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_5mi_int2 <- sem(model_bifica_5mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_5mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
      
    #Interaction 3: Crime + Sex*Age 
      model_bifica_5mi_int3 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C*TotalCrime_5z + Age + Sex + Age_Sex + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A*TotalCrime_5z + Age + Sex + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A*B
                            
                            #Total Effect 
                            TotalEffect := C + (A*B)
                          '
      
      mediationfit_bifica_5mi_int3 <- sem(model_bifica_5mi_int3, data = g) 
      #summary(mediationfit_bifica_5mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_5mi_int3 <- sem(model_bifica_5mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_5mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
  #Interaction 4: Crime*Sex*Age 
      model_bifica_5mi_int4 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_5z + C2*Age + C3*Sex_num + C4*Crime5z_Sex_Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~  A1*TotalCrime_5z + A2*Age + A3*Sex_num + A4*Crime5z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A4*B
                            
                            #Total Effect 
                            TotalEffect := C4 + (A4*B)
                          '
      
      mediationfit_bifica_5mi_int4 <- sem(model_bifica_5mi_int4, data = g) 
      #summary(mediationfit_bifica_5mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_5mi_int4 <- sem(model_bifica_5mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_5mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
    
  #Hypothesis 5: 10-mile risk & FCCA  
      #Unadjusted 
      model <- ' #Direct effect between X and Y
                CMMFCCAV1_w ~ C*TotalCrime_10z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z
    
                #Path B
                CMMFCCAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
 
  #Adjusted linear regression for mediation model 
      model <- ' #Direct effect between X and Y
                CMMFCCAV1_w ~ C*TotalCrime_10z + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
    
      
    #Interaction 1: Crime*Sex 
      model_fcca_10mi_int1 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C1*TotalCrime_10z + C2*Sex_num + C3*Crime10z_Sex + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Sex_num + A3*Crime10z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A3*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C3 + (A3*B)
                '
      
      mediationfit_fcca_10mi_int1 <- sem(model_fcca_10mi_int1, data = g) 
      #summary(mediationfit_10mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_fcca_10mi_int1 <- sem(model_fcca_10mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      summary(mediationfit_fcca_10mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
   
  #Interaction 2: Crime*Age
      model_fcca_10mi_int2 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Crime10z_Age + Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Crime10z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A3*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C3 + (A3*B)
                '
      
      mediationfit_fcca_10mi_int2 <- sem(model_fcca_10mi_int2, data = g) 
      #summary(mediationfit_fcca_10mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_fcca_10mi_int2 <- sem(model_fcca_10mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      summary(mediationfit_fcca_10mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
    #Interaction 3: Crime + Sex*Age
      model_fcca_10mi_int3 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C*TotalCrime_10z + Age + Sex_num + Age_Sex + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex_num + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
      
      mediationfit_fcca_10mi_int3 <- sem(model_fcca_10mi_int3, data = g) 
      #summary(mediationfit_fcca_10mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_fcca_10mi_int3 <- sem(model_fcca_10mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
  
    #Interaction 4: Crime*Sex*Age
      model_fcca_10mi_int4 <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Sex_num + C4*Crime10z_Sex_Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K
                  
                  #Path A 
                  mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Sex_num + A4*Crime10z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A4*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C4 + (A4*B)
                '
      
      mediationfit_fcca_10mi_int4 <- sem(model_fcca_10mi_int4, data = g) 
      #summary(mediationfit_10mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_fcca_10mi_int4 <- sem(model_fcca_10mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      summary(mediationfit_fcca_10mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
  #Hypothesis 5: 10-mile risk & BIF  
      #Unadjusted
      model <- ' #Direct effect between X and Y
                CMMBIFV1_w ~ C*TotalCrime_10z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z
    
                #Path B
                CMMBIFV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
    #Adjusted
      model <- ' #Direct effect between X and Y
                CMMBIFV1_w ~ C*TotalCrime_10z + Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMBIFV1T1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
      
    #Interaction 1: Crime*Sex 
      model_bif_10mi_int1<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_10z + C2*Sex_num + C3*Crime10z_Sex + Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Sex_num + A3*Crime10z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A3*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C3 + (A3*B)
                                '
      
      mediationfit_bif_10mi_int1 <- sem(model_bif_10mi_int1, data = g) 
      #summary(mediationfit_bif_10mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_10mi_int1 <- sem(model_bif_10mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_10mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE) 
      
   
    #Interaction 2: Crime*Age
      model_bif_10mi_int2<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Crime10z_Age + Sex_num + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Crime10z_Age + Sex_num + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A3*B
                                  
                      
                                  #Total Effect 
                                  TotalEffect := C3 + (A3*B)
                                '
      
      mediationfit_bif_10mi_int2 <- sem(model_bif_10mi_int2, data = g) 
      #summary(mediationfit_bif_10mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_10mi_int2 <- sem(model_bif_10mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_10mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)   
     
      
    #Interaction 3: Crime + Sex*Age
      model_bif_10mi_int3<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C*TotalCrime_10z + Age + Sex_num + + Age_Sex + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex_num + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C + (A*B)
                                '
      
      mediationfit_bif_10mi_int3 <- sem(model_bif_10mi_int3, data = g) 
      #summary(mediationfit_bif_10mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_10mi_int3 <- sem(model_bif_10mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_10mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
  
      
    #Interaction 4: Crime*Sex*Age
      model_bif_10mi_int4<- ' #Direct effect between X and Y
                                  CMMBIFV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Sex_num + C4*Crime10z_Sex_Age + MMBIFV1T1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                  
                                  #Path A 
                                  mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Sex_num + A4*Crime10z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                      
                                  #Path B
                                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFV1T1 + BMIT1
                            
                                  #Indirect Effect
                                  IndirectEffect := A4*B
                                  #Creating new parameter testing for the indirect effects
                      
                                  #Total Effect 
                                  TotalEffect := C4 + (A4*B)
                                '
      
      mediationfit_bif_10mi_int4 <- sem(model_bif_10mi_int4, data = g) 
      #summary(mediationfit_bif_10mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bif_10mi_int4 <- sem(model_bif_10mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bif_10mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
    
  #Hypothesis 5: 10-mile risk & ICA 
      #Unadjusted
      model <- ' #Direct effect between X and Y
                CMMICAV1_w ~ C*TotalCrime_10z
                #C* means I am giving the parameter the name C
    
                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z
    
                #Path B
                CMMICAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
      
    #Adjusted
      model <- ' #Direct effect between X and Y
                CMMICAV1_w ~ C*TotalCrime_10z + Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                #C* means I am giving the parameter the name C
    
                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)

    #Interaction 1: Crime*Sex
      #library(lavaan)
      model_ica_10mi_int1 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_10z + C2*Sex_num + C3*Crime10z_Sex + Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Sex_num + A3*Crime10z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A3*B
    
                                      #Total Effect 
                                      TotalEffect := C3 + (A3*B)
                                    '
      
      mediationfit_ica_10mi_int1 <- sem(model_ica_10mi_int1, data = g) 
      #summary(mediationfit_ica_10mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_10mi_int1 <- sem(model_ica_10mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_10mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
      
    #Interaction 2: Crime*Age
      model_ica_10mi_int2 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Crime10z_Age + Sex + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Crime10z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A3*B
    
                                      #Total Effect 
                                      TotalEffect := C3 + (A3*B)
                                    '
      
      mediationfit_ica_10mi_int2 <- sem(model_ica_10mi_int2, data = g) 
      #summary(mediationfit_ica_10mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_10mi_int2 <- sem(model_ica_10mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_10mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
     
    #Interaction 3: Crime + Age*Sex
      model_ica_10mi_int3 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C*TotalCrime_10z + Age + Sex + Age_Sex + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A*B
    
                                      #Total Effect 
                                      TotalEffect := C + (A*B)
                                    '
      
      mediationfit_ica_10mi_int3 <- sem(model_ica_10mi_int3, data = g) 
      #summary(mediationfit_ica_10mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_10mi_int3 <- sem(model_ica_10mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_10mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
    
    #Interaction 4: Crime*Sex*Age
      model_ica_10mi_int4 <- ' #Direct effect between X and Y
                                      CMMICAV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Sex_num + C4*Crime10z_Sex_Age + MMICAV1T1 + MedLipT1 + MedDT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                                      
                                      #Path A 
                                      mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Sex_num + A4*Crime10z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                          
                                      #Path B
                                      CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMICAV1T1 + MedLipT1 + MedDT1 + BMIT1
                                
                                      #Indirect Effect
                                      IndirectEffect := A4*B
    
                                      #Total Effect 
                                      TotalEffect := C4 + (A4*B)
                                    '
      
      mediationfit_ica_10mi_int4 <- sem(model_ica_10mi_int4, data = g) 
      #summary(mediationfit_ica_10mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_ica_10mi_int4 <- sem(model_ica_10mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_ica_10mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
    
#Hypothesis 5: 10-mile risk & BIF/ICA 
      #Unadjusted 
      model <- ' #Direct effect between X and Y
                CCBIFICAV1_w ~ C*TotalCrime_10z

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z
    
                #Path B
                CCBIFICAV1_w ~ B*mEMASVQTV_boxcox
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
    #Adjusted
      model <- ' #Direct effect between X and Y
                CCBIFICAV1_w ~ C*TotalCrime_10z + Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1

                #Path A 
                mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
    
                #Path B
                CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + MMBIFICAV1T1 + DBPT1 + BMIT1
          
                #Indirect Effect
                IndirectEffect := A*B
                #Creating new parameter testing for the indirect effects
    
                #Total Effect 
                TotalEffect := C + (A*B)
              '
      
      mediationfit <- sem(model, data = g) 
      #summary(mediationfit)
      
      # Modify the sem function to perform bootstrapping
      mediationfit <- sem(model, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
      
  #Interaction 1: Crime*Sex 
      model_bifica_10mi_int1 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_10z + C2*Sex_num + C3*Crime10z_Sex + Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Sex_num + A3*Crime10z_Sex + Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A3*B
                            #Creating new parameter testing for the indirect effects
                
                            #Total Effect 
                            TotalEffect := C3 + (A3*B)
                          '
      
      mediationfit_bifica_10mi_int1 <- sem(model_bifica_10mi_int1, data = g) 
      #summary(mediationfit_bifica_10mi_int1)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_10mi_int1 <- sem(model_bifica_10mi_int1, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_10mi_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
      
  
    #Interaction 2: Crime*Age 
      model_bifica_10mi_int2 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Crime10z_Age + Sex + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A1*TotalCrime_10z + A2*Age + A3*Crime10z_Age + Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A3*B
                            #Creating new parameter testing for the indirect effects
                
                            #Total Effect 
                            TotalEffect := C3 + (A3*B)
                          '
      
      mediationfit_bifica_10mi_int2 <- sem(model_bifica_10mi_int2, data = g) 
      #summary(mediationfit_bifica_10mi_int2)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_10mi_int2 <- sem(model_bifica_10mi_int2, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_10mi_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
      
    #Interaction 3: Crime + Sex*Age 
      model_bifica_10mi_int3 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C*TotalCrime_10z + Age + Sex + Age_Sex + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~ A*TotalCrime_10z + Age + Sex + Age_Sex + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A*B
                            
                            #Total Effect 
                            TotalEffect := C + (A*B)
                          '
      
      mediationfit_bifica_10mi_int3 <- sem(model_bifica_10mi_int3, data = g) 
      #summary(mediationfit_bifica_10mi_int3)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_10mi_int3 <- sem(model_bifica_10mi_int3, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_10mi_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
      
    #Interaction 4: Crime*Sex*Age 
      model_bifica_10mi_int4 <- ' #Direct effect between X and Y
                            CCBIFICAV1_w ~ C1*TotalCrime_10z + C2*Age + C3*Sex_num + C4*Crime10z_Sex_Age + MMBIFICAV1T1 + DBPT1 + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K + BMIT1
                           
                            #Path A 
                            mEMASVQTV_boxcox ~  A1*TotalCrime_10z + A2*Age + A3*Sex_num + A4*Crime10z_Sex_Age + ref_NHB + ref_NHO + ref_HL + Income_below_30K + Income_30K_100K 
                
                            #Path B
                            CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + MMBIFICAV1T1 + DBPT1 + BMIT1
                      
                            #Indirect Effect
                            IndirectEffect := A4*B
                            
                            #Total Effect 
                            TotalEffect := C4 + (A4*B)
                          '
      
      mediationfit_bifica_10mi_int4 <- sem(model_bifica_10mi_int4, data = g) 
      #summary(mediationfit_bifica_10mi_int4)
      
      # Modify the sem function to perform bootstrapping
      mediationfit_bifica_10mi_int4 <- sem(model_bifica_10mi_int4, data = g, se = "bootstrap", bootstrap = 5000)
      
      # Display summary with confidence intervals
      summary(mediationfit_bifica_10mi_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
            
#Hypothesis 6: Social vigilance will mediate or partially mediate any observed effect of social safety on 2-year change in cIMT. 

      
    #Hypothesis 6: FCCA 
        #Unadjusted 
                  model <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C*socsaf_mean

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                #summary(mediationfit)
      
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE) 
     
          
          #Adjusted 
                model <- ' #Direct effect between X and Y
                  CMMFCCAV1_w ~ C*socsaf_mean + Age + Sex + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean + Age + Sex 
      
                  #Path B
                  CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex+ SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
            
      
        #Interactions: FCCA mediation model with SS*Sex and adjustment without age l 
              model_int1 <- '
              # Direct effect
              CMMFCCAV1_w ~ C1*socsaf_mean + C2*Sex + C3*socsaf_Sex + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex + A3*socsaf_Sex 
            
              # Path B
              CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int1 <- sem(model_int1, data = f) 
                summary(mediationfit_int1)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int1 <- sem(model_int1, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
        #Interactions: FCCA mediation model with SS*Sex and adjustment with age
                model_int2 <- '
              # Direct effect
              CMMFCCAV1_w ~ C1*socsaf_mean + C2*Sex_num + C3*socsaf_Sex + Age + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex_num + A3*socsaf_Sex + Age
            
              # Path B
              CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Age + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int2 <- sem(model_int2, data = f) 
                summary(mediationfit_int2)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int2 <- sem(model_int2, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
        #Interactions: FCCA mediation model with SS*Age and adjustment without sex
              model_int3 <- '
              # Direct effect
              CMMFCCAV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age
            
              # Path B
              CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int3 <- sem(model_int3, data = f) 
                summary(mediationfit_int3)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int3 <- sem(model_int3, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
                
        #Interactions: FCCA mediation model with SS*Age and adjustment with sex
              model_int4 <- '
              # Direct effect
              CMMFCCAV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + Sex + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age + Sex
            
              # Path B
              CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + Sex + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int4 <- sem(model_int4, data = f) 
                summary(mediationfit_int4)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int4 <- sem(model_int4, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
  
        #Interactions: FCCA mediation model with SS and age*sex and adjustments
              model_int5 <- '
              # Direct effect
              CMMFCCAV1_w ~ C*socsaf_mean + Age_Sex + Age + Sex + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Path A
              mEMASVQTV_boxcox ~ A*socsaf_mean + Age_Sex + Age + Sex 
            
              # Path B
              CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Indirect effect of interaction
              IndirectEffect := A*B
            
              # Total effect of interaction
              TotalEffect := C + (A*B)
            '
                mediationfit_int5 <- sem(model_int5, data = f) 
                summary(mediationfit_int5)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int5 <- sem(model_int5, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int5, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
      #Interactions: FCCA mediation model with SS*age*sex and adjustments
            model_int6 <- '
              # Direct effect
              CMMFCCAV1_w ~ C1*socsaf_mean + C2*Age + C3*Sex_num + C4*socsaf_Sex_Age + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*Sex_num + A4*socsaf_Sex_Age 
            
              # Path B
              CMMFCCAV1_w ~ B*mEMASVQTV_boxcox + SBPT1 + DBPT1 + HIncome + MMFCCAV1T1 + MedLipT1
            
              # Indirect effect of interaction
              IndirectEffect := A4*B
            
              # Total effect of interaction
              TotalEffect := C4 + (A4*B)
            '
                mediationfit_int6 <- sem(model_int6, data = f) 
                summary(mediationfit_int6)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int6 <- sem(model_int6, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int6, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
  #Hypothesis 6: BIF 
        #Unadjusted  
                  model <- ' #Direct effect between X and Y
                  CMMBIFV1_w ~ C*socsaf_mean

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean
      
                  #Path B
                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
  
        #Adjusted 
                model <- ' #Direct effect between X and Y
                  CMMBIFV1_w ~ C*socsaf_mean + Age + Sex + BMIT1 + HIncome + MMBIFV1T1 

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean + Age + Sex 
      
                  #Path B
                  CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + BMIT1 + HIncome + MMBIFV1T1 
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                

      #Interactions: BIF mediation model with SS*Sex and adjustment without age 
              model_int1 <- '
              # Direct effect
              CMMBIFV1_w ~ C1*socsaf_mean + C2*Sex + C3*socsaf_Sex + BMIT1 + HIncome + MMBIFV1T1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex + A3*socsaf_Sex 
            
              # Path B
              CMMBIFV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFV1T1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int1 <- sem(model_int1, data = f) 
                summary(mediationfit_int1)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int1 <- sem(model_int1, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                

      #Interactions: BIF mediation model with SS*Sex and adjustment with age
              model_int2 <- '
              # Direct effect
              CMMBIFV1_w ~ C1*socsaf_mean + C2*Sex_num + C3*socsaf_Sex + Age + BMIT1 + HIncome + MMBIFV1T1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex_num + A3*socsaf_Sex + Age
            
              # Path B
              CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Age + BMIT1 + HIncome + MMBIFV1T1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int2 <- sem(model_int2, data = f) 
                summary(mediationfit_int2)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int2 <- sem(model_int2, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
              
      
      #Interactions: BIF mediation model with SS*Age and adjustment without sex
              model_int3 <- '
              # Direct effect
              CMMBIFV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + BMIT1 + HIncome + MMBIFV1T1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age
            
              # Path B
              CMMBIFV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFV1T1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int3 <- sem(model_int3, data = f) 
                summary(mediationfit_int3)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int3 <- sem(model_int3, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
   
      #Interactions: BIF mediation model with SS*Age and adjustment with sex
              model_int4 <- '
              # Direct effect
              CMMBIFV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + Sex + BMIT1 + HIncome + MMBIFV1T1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age + Sex 
            
              # Path B
              CMMBIFV1_w ~ B*mEMASVQTV_boxcox + Sex + BMIT1 + HIncome + MMBIFV1T1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int4 <- sem(model_int4, data = f) 
                summary(mediationfit_int4)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int4 <- sem(model_int4, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
      #Interactions: BIF mediation model with SS and Age*Sex and adjustments
              model_int5 <- '
              # Direct effect
              CMMBIFV1_w ~ C*socsaf_mean + Age_Sex + Age + Sex + BMIT1 + HIncome + MMBIFV1T1 
            
              # Path A
              mEMASVQTV_boxcox ~ A*socsaf_mean + Age_Sex + Age + Sex
            
              # Path B
              CMMBIFV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFV1T1 
            
              # Indirect effect of interaction
              IndirectEffect := A*B
            
              # Total effect of interaction
              TotalEffect := C + (A*B)
            '
                mediationfit_int5 <- sem(model_int5, data = f) 
                summary(mediationfit_int5)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int5 <- sem(model_int5, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int5, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                

      #Interactions: BIF mediation model with SS*age*sex and adjustments
              model_int6 <- '
              # Direct effect
              CMMBIFV1_w ~ C1*socsaf_mean + C2*Age + C3*Sex_num + C4*socsaf_Sex_Age + BMIT1 + HIncome + MMBIFV1T1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*Sex_num + A4*socsaf_Sex_Age 
            
              # Path B
              CMMBIFV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFV1T1 
            
              # Indirect effect of interaction
              IndirectEffect := A4*B
            
              # Total effect of interaction
              TotalEffect := C4 + (A4*B)
            '
                mediationfit_int6 <- sem(model_int6, data = f) 
                summary(mediationfit_int6)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int6 <- sem(model_int6, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int6, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
 
  #Hypothesis 6: ICA 
      #Unadjusted 
                model <- ' #Direct effect between X and Y
                  CMMICAV1_w ~ C*socsaf_mean

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean
      
                  #Path B
                  CMMICAV1_w ~ B*mEMASVQTV_boxcox
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
              
          
        #Adjusted 
                model <- ' #Direct effect between X and Y
                  CMMICAV1_w ~ C*socsaf_mean + Age + Sex + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean + Age + Sex 
      
                  #Path B
                  CMMICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
      
      #Interactions: ICA mediation model with SS*Sex and adjustment without age 
              model_int1 <- '
              # Direct effect
              CMMICAV1_w ~ C1*socsaf_mean + C2*Sex_num + C3*socsaf_Sex + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex_num + A3*socsaf_Sex 
            
              # Path B
              CMMICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int1 <- sem(model_int1, data = f) 
                summary(mediationfit_int1)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int1 <- sem(model_int1, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
       
      #Interactions: ICA mediation model with SS*Sex and adjustment with age 
              model_int2 <- '
              # Direct effect
              CMMICAV1_w ~ C1*socsaf_mean + C2*Sex_num + C3*socsaf_Sex + Age + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex_num + A3*socsaf_Sex 
            
              # Path B
              CMMICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int2 <- sem(model_int2, data = f) 
                summary(mediationfit_int2)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int2 <- sem(model_int2, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                

      #Interactions: ICA mediation model with SS*Age and adjustment without sex  
                model_int3 <- '
              # Direct effect
              CMMICAV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age 
            
              # Path B
              CMMICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int3 <- sem(model_int3, data = f) 
                summary(mediationfit_int3)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int3 <- sem(model_int3, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
  
      #Interactions: ICA mediation model with SS*Age and adjustment with sex  
                model_int4 <- '
              # Direct effect
              CMMICAV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + Sex + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age + Sex
            
              # Path B
              CMMICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int4 <- sem(model_int4, data = f) 
                summary(mediationfit_int4)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int4 <- sem(model_int4, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
              
      
      #Interactions: ICA mediation model with SS and age*sex and adjustments 
                model_int5 <- '
              # Direct effect
              CMMICAV1_w ~ C*socsaf_mean + Age_Sex + Age + Sex + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Path A
              mEMASVQTV_boxcox ~ A*socsaf_mean + Age_Sex + Age + Sex
            
              # Path B
              CMMICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Indirect effect of interaction
              IndirectEffect := A*B
            
              # Total effect of interaction
              TotalEffect := C + (A*B)
            '
                mediationfit_int5 <- sem(model_int5, data = f) 
                summary(mediationfit_int5)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int5 <- sem(model_int5, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int5, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
      #Interactions: ICA mediation model with SS*age*sex and adjustments 
                model_int6 <- '
              # Direct effect
              CMMICAV1_w ~ C1*socsaf_mean + C2*Age + C3*Sex_num + C4*socsaf_Sex_Age + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*Sex_num + A4*socsaf_Sex_Age
            
              # Path B
              CMMICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMICAV1T1 + MedLipT1 + MedDT1 
            
              # Indirect effect of interaction
              IndirectEffect := A4*B
            
              # Total effect of interaction
              TotalEffect := C4 + (A4*B)
            '
                mediationfit_int6 <- sem(model_int6, data = f) 
                summary(mediationfit_int6)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int6 <- sem(model_int6, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int6, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
  #Hypothesis 6: BIF/ICA 
         #Unadjusted

                model <- ' #Direct effect between X and Y
                  CCBIFICAV1_w ~ C*socsaf_mean

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean
      
                  #Path B
                  CCBIFICAV1_w ~ B*mEMASVQTV_boxcox
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
              
                
        #Adjusted
                
                model <- ' #Direct effect between X and Y
                  CCBIFICAV1_w ~ C*socsaf_mean + Age + Sex + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1 

                  #Path A 
                  mEMASVQTV_boxcox ~ A*socsaf_mean + Age + Sex 
      
                  #Path B
                  CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + Sex + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1 
            
                  #Indirect Effect
                  IndirectEffect := A*B
                  #Creating new parameter testing for the indirect effects
      
                  #Total Effect 
                  TotalEffect := C + (A*B)
                '
                mediationfit <- sem(model, data = f) 
                summary(mediationfit)
                
                # Modify the sem function to perform bootstrapping
                mediationfit <- sem(model, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
          
      #Interactions: BIF/ICA mediation model with SS*Sex and adjustment without age 
              model_int1 <- '
              # Direct effect
              CCBIFICAV1_w ~ C1*socsaf_mean + C2*Sex_num + C3*socsaf_Sex + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex_num + A3*socsaf_Sex 
            
              # Path B
              CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int1 <- sem(model_int1, data = f) 
                summary(mediationfit_int1)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int1 <- sem(model_int1, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int1, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                

      #Interactions: BIF/ICA mediation model with SS*Sex and adjustment with age 
                model_int2 <- '
              # Direct effect
              CCBIFICAV1_w ~ C1*socsaf_mean + C2*Sex_num + C3*socsaf_Sex + Age + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Sex_num + A3*socsaf_Sex + Age
            
              # Path B
              CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Age + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int2 <- sem(model_int2, data = f) 
                summary(mediationfit_int2)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int2 <- sem(model_int2, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int2, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
   
      #Interactions: BIF/ICA mediation model with SS*Age and adjustment without sex 
                model_int3 <- '
              # Direct effect
              CCBIFICAV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age
            
              # Path B
              CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int3 <- sem(model_int3, data = f) 
                summary(mediationfit_int3)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int3 <- sem(model_int3, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int3, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
   
      #Interactions: BIF/ICA mediation model with SS*Age and adjustment with sex 
                model_int4 <- '
              # Direct effect
              CCBIFICAV1_w ~ C1*socsaf_mean + C2*Age + C3*socsaf_Age + Sex + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*socsaf_Age + Sex
            
              # Path B
              CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + Sex + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Indirect effect of interaction
              IndirectEffect := A3*B
            
              # Total effect of interaction
              TotalEffect := C3 + (A3*B)
            '
                mediationfit_int4 <- sem(model_int4, data = f) 
                summary(mediationfit_int4)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int4 <- sem(model_int4, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int4, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
                
      #Interactions: BIF/ICA mediation model with SS + Age*Sex + adjustments
                model_int5 <- '
              # Direct effect
              CCBIFICAV1_w ~ C*socsaf_mean + Age_Sex + Age + Sex + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Path A
              mEMASVQTV_boxcox ~ A*socsaf_mean + Age_Sex + Age + Sex
            
              # Path B
              CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Indirect effect of interaction
              IndirectEffect := A*B
            
              # Total effect of interaction
              TotalEffect := C + (A*B)
            '
                mediationfit_int5 <- sem(model_int5, data = f) 
                summary(mediationfit_int5)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int5 <- sem(model_int5, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int5, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
 
      #Interactions: BIF/ICA mediation model with SS*Age*Sex + adjustments
              model_int6 <- '
              # Direct effect
              CCBIFICAV1_w ~ C1*socsaf_mean + C2*Age + C3*Sex_num + C4*socsaf_Sex_Age + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Path A
              mEMASVQTV_boxcox ~ A1*socsaf_mean + A2*Age + A3*Sex_num + A4*socsaf_Sex_Age
            
              # Path B
              CCBIFICAV1_w ~ B*mEMASVQTV_boxcox + BMIT1 + HIncome + MMBIFICAV1T1 + DBPT1  
            
              # Indirect effect of interaction
              IndirectEffect := A4*B
            
              # Total effect of interaction
              TotalEffect := C4 + (A4*B)
            '
                mediationfit_int6 <- sem(model_int6, data = f) 
                summary(mediationfit_int6)
                
                # Modify the sem function to perform bootstrapping
                mediationfit_int6 <- sem(model_int6, data = f, se = "bootstrap", bootstrap = 5000)
                
                # Display summary with confidence intervals
                summary(mediationfit_int6, fit.measures = TRUE, standardized = TRUE, ci = TRUE)
                
 
 
#Age-Sex stratified descriptives---- 
  #create Age stratified by decade (Age_stra)
    f$Age_stra <- ifelse(f$Age < 31, 1,
                         ifelse(f$Age < 41, 2,
                          ifelse(f$Age < 51, 3, 
                               ifelse(f$Age < 61, 4, 
                                       ifelse(f$Age < 71, 5, -999)))))
    table(f$Age_stra)
  
  f$Age_straf <- factor(f$Age_stra, levels = c(1,2,3,4,5), labels = c("21 - 30", "31 - 40", "41 - 50", "51 - 60", "61 - 70"))
  table(f$Age_straf)


#Calculate mean FCCA IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMFCCAV1T1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of FCCA IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMFCCAV1T1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of FCCA IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMFCCAV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid FCCA IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMFCCAV1T1)))


#Calculate mean BIF IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMBIFV1T1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BIF IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMBIFV1T1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BIF IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMBIFV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid BIF IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMBIFV1T1)))
  
#Calculate mean BIF IMT T1 (not winsorized) by sex 
  result <- aggregate(MMBIFV1T1 ~ Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BIF IMT T1 (not winsorized) by sex 
  result <- aggregate(MMBIFV1T1 ~ Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BIF IMT T1 (not winsorized) by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMBIFV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid BIF IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMBIFV1T1)))
  
#Calculate mean BIF IMT T1 (winsorized) by Age and sex 
  result <- aggregate(MMBIFV1T1_w ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BIF IMT T1 (winsorized) by Age and sex 
  result <- aggregate(MMBIFV1T1_w ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BIF IMT T1 (winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMBIFV1T1_w, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid BIF IMT T1 (winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMBIFV1T1_w)))
  
#Calculate mean ICA IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMICAV1T1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of ICA IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMICAV1T1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of ICA IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMICAV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid ICA IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMICAV1T1)))
  
  
#Calculate mean ICA IMT T1 (not winsorized) by sex 
  result <- aggregate(MMICAV1T1 ~ Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of ICA IMT T1 (not winsorized) by sex 
  result <- aggregate(MMICAV1T1 ~ Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of ICA IMT T1 (not winsorized) by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMICAV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid ICA IMT T1 (not winsorized) by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMICAV1T1)))
  
  
#Calculate mean ICA IMT T1 (winsorized) by Age and sex 
  result <- aggregate(MMICAV1T1_w ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of ICA IMT T1 (winsorized) by Age and sex 
  result <- aggregate(MMICAV1T1_w ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of ICA IMT T1 (winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMICAV1T1_w, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid ICA IMT T1 (winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMICAV1T1_w)))
  
#Calculate mean BIF/ICA IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMBIFICAV1T1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BIF/ICA IMT T1 (not winsorized) by Age and sex 
  result <- aggregate(MMBIFICAV1T1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BIF/ICA IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMBIFICAV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid BIF/ICA IMT T1 (not winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMBIFICAV1T1)))
  
#Calculate mean BIF/ICA IMT T1 (not winsorized) by sex 
  result <- aggregate(MMBIFICAV1T1 ~ Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BIF/ICA IMT T1 (not winsorized) by sex 
  result <- aggregate(MMBIFICAV1T1 ~ Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BIF/ICA IMT T1 (not winsorized) by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMBIFICAV1T1, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid BIF/ICA IMT T1 (not winsorized) by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMBIFICAV1T1)))
  
#Calculate mean BIF/ICA IMT T1 (winsorized) by Age and sex 
  result <- aggregate(MMBIFICAV1T1_w ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BIF/ICA IMT T1 (winsorized) by Age and sex 
  result <- aggregate(MMBIFICAV1T1_w ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BIF/ICA IMT T1 (winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(MMBIFICAV1T1_w, probs = 0.9, na.rm = TRUE)
    )
  
  
  #Calculate n valid BIF/ICA IMT T1 (winsorized) by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(MMBIFICAV1T1_w)))
  
#Calculate mean BMI @ T1 by Age and sex 
  result <- aggregate(BMIT1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BMI @ T1 by Age and sex 
  result <- aggregate(BMIT1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BMI @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(BMIT1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid BMI @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(BMIT1)))
  
#Calculate mean BMI @ T1 by sex 
  result <- aggregate(BMIT1 ~ Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of BMI @ T1 by sex 
  result <- aggregate(BMIT1 ~ Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of BMI @ T1 by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(
      quantile_90 = quantile(BMIT1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid BMI @ T1 by  sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(non_na_count = sum(!is.na(BMIT1)))
  
  
#Calculate mean SBP @ T1 by Age and sex 
  result <- aggregate(SBPT1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of SBP @ T1 by Age and sex 
  result <- aggregate(SBPT1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of SBP @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(SBPT1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid SBP @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(SBPT1)))
  
#Calculate mean SBP @ T1 by sex 
  result <- aggregate(SBPT1 ~ Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of SBP @ T1 by sex 
  result <- aggregate(SBPT1 ~ Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of SBP @ T1 by  sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(
      quantile_90 = quantile(SBPT1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid SBP @ T1 by sex 
  f %>%
    group_by(Sex) %>%
    dplyr::summarize(non_na_count = sum(!is.na(SBPT1)))
  
#Calculate mean DBP @ T1 by Age and sex 
  result <- aggregate(DBPT1 ~ Age_straf + Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of DBP @ T1 by Age and sex 
  result <- aggregate(DBPT1 ~ Age_straf + Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of DBP @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(DBPT1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid DBP @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(DBPT1)))
  
#Calculate mean DBP @ T1 by sex 
  result <- aggregate(DBPT1 ~ Sex, data = f, FUN = mean)
  
  # Print the result
  print(result)
  
  #Calculate standard deviation of DBP @ T1 by sex 
  result <- aggregate(DBPT1 ~ Sex, data = f, FUN = sd)
  
  # Print the result
  print(result)
  
  #Calculate 90th percentile of DBP @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(
      quantile_90 = quantile(DBPT1, probs = 0.9, na.rm = TRUE)
    )
  
  #Calculate n valid DBP @ T1 by Age and sex 
  f %>%
    group_by(Sex, Age_straf) %>%
    dplyr::summarize(non_na_count = sum(!is.na(DBPT1)))
  
#Calculate Lipid Med frequency by Age and sex 
  table(f$MedLipT1, f$Age_straf, f$Sex)

#Calculate Diabetes Med frequency by Age and sex 
  table(f$MedDT1, f$Age_straf, f$Sex)
  
#Calcuate freq of Diabetes med and or lipid med 

  table(f$MedDT1)
  table(f$MedLipT1)
  

  table(f$MedLipT1, f$MedDT1)
  #46 ppl take diabetes and/or lipid mgmt meds 
  #15.33% of N = 300 

#Creating Stage I Hypertension factor variable (140 > SBP > 129 OR 90 > DBP > 79)
  f$Hypt_S1 <- ifelse(
    is.na(f$SBPT1) | is.na(f$DBPT1), NA,
    ifelse((f$SBPT1 > 129 & f$SBPT1 < 140) | (f$DBPT1 > 79 & f$DBPT1 < 90), 1, 0)
  )
  
  table(f$Hypt_S1)
  table(f$Hypt_S1, f$Age_straf, f$Sex)
  
#Calculate Stage I Hypertension frequency by Age and sex 
  table(f$Hypt_S1, f$Age_straf, f$Sex)
  
#Creating Stage II Hypertension factor variable (SBP >= 140 OR DBP >= 90)
  f$Hypt_S2 <- ifelse(is.na(f$SBPT1) | is.na(f$DBPT1), NA, 
                      ifelse(f$SBPT1 > 139 | f$DBPT1 > 89, 1, 0))
  
  table(f$Hypt_S2)
  
  #Calculate Stage II Hypertension frequency by Age and sex 
  table(f$Hypt_S2, f$Age_straf, f$Sex)
  
#Creating three-level Hypertension factor variable (
    #0 = not hypertensive 
    #1 = stage 1 HT: 140 > SBP > 129 OR 90 > DBP > 79
    #2 = stage 2 HT: SBP >= 140 OR DBP >= 90 
  f$Hypt <- ifelse(
    is.na(f$SBPT1) | is.na(f$DBPT1), NA,
    ifelse((f$SBPT1 > 129 & f$SBPT1 < 140) | (f$DBPT1 > 79 & f$DBPT1 < 90), 1,
           ifelse(f$SBPT1 > 139 | f$DBPT1 > 89, 2, 0)))
  
  table(f$Hypt)
  table(f$Hypt, f$Age_straf, f$Sex)

