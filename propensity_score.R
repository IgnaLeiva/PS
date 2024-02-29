library(tidyverse)
library(cobalt)
# Fit model
dat <- readRDS("final_data.Rds")
#class
sapply(dat, class)
#factor

# Transform
#dat <- dat %>% 
#mutate_if(names(dat) %in% factor_vars, factor)

# Model
#Data with all variables include in the regression plus patid
#dat_regression <- subset(dat, select = -c(patid, indexdate, pracid, deathdate, enddate, died, time))
#dat_imp <- readRDS("data_imp.Rds")
#variables_model <- subset(dat, select = -c(patid))

#################################################################
##             Propensity Score.                               ##
##                                                             ##
#################################################################                          

#----------------------                        
# Logistic regression
#----------------------

vars <- c("age","gender", "eth5", "imd_person", "calendarperiod", "prior_gast_cancer", "recent_gerd", "recent_peptic", "bmi", "n_consult","pracid")
formula <- as.formula(paste0("ppi~", paste(vars, collapse = "+")))
                      
#----------------------                        
# Logistic regression
#----------------------


PS_model <- glm(formula,  family=binomial(), data= dat) 
# Display model coefficients
#model_out <- broom::tidy(PS_model, exponentiate=TRUE, conf.int=TRUE) 

# add PS in the dataset
dat$PS <- PS_model %>% predict(dat,  type = "response") 

#histogram
hist(dat$PS, xlab="Propensity score", ylab= "Number of individuals" , main="Propensity score distribution", col="darkgrey", breaks=50) 

#boxplot
boxplot(PS ~ ppi, data= dat, ylab="Propensity score", xlab="PPI prescription") 

# histo
ggplot(dat, aes(x=PS, color=factor(ppi),fill = factor(ppi))) +  
  geom_histogram(position="identity",bins = 50, alpha=0.5)+    
  labs(x = "Propensity score", y="Number of individuals" )

ggplot(dat, aes(x=PS, color=factor(ppi),fill = factor(ppi))) +  
  geom_density(position="identity",bins = 50, alpha=0.5)+    
  labs(x = "Propensity score", y="Number of individuals" )


#----------------------                        
# Twang package
#----------------------
library(twang)

# Data need to be numeric.
dat_numeric <- dat
dat_numeric <- dat_numeric |> 
  mutate(ppi= ifelse(ppi == 0, 0, 1),
         gender = ifelse(gender == "Female", 0, 1),
         eth5 = case_when(
           eth5 == "White" ~ 0,
           eth5 == "South-Asian" ~ 1,
           eth5 == "Black" ~ 2,
           eth5 == "Unknown" ~ 3,
           eth5 == "Mixed" ~ 4,
           eth5 == "Other" ~ 5),
         imd_person = case_when(
           imd_person == "Least Deprived (1)" ~ 0,
           imd_person == 2 ~ 1,
           imd_person == 3 ~ 2,
           imd_person == 4 ~ 3,
           imd_person == "Most Deprived (5)" ~ 4),
         calendarperiod = case_when(
           calendarperiod == "1991-1999" ~ 0,
           calendarperiod == "2000-2004" ~ 1,
           calendarperiod == "2005-2009" ~ 2,
           calendarperiod == "2010-2014" ~ 3,
           calendarperiod == "2015-2017"~ 4),
         prior_gast_cancer = ifelse(prior_gast_cancer == 0, 0, 1),
         recent_gerd = ifelse(recent_gerd == 0, 0, 1),
         recent_peptic = ifelse(recent_peptic == 0, 0, 1))


ps.gbm = ps(formula,
              data = as.data.frame(dat),
              n.trees=5000,
              
              interaction.depth=2,
              shrinkage=0.01,
              estimand = "ATE",
              stop.method=c("es.mean","ks.max"),
              n.minobsinnode = 10,
              n.keep = 1,
              n.grid = 25,
              ks.exact = NULL,
              verbose=FALSE)


dat$ps_m <- ps.gbm$ps$es.mean.ATE


#################################################################
##             Assessing the covariates                        ##
##                     SMD                                    ##
#################################################################

# Assessing the covariates

dat %>%
  group_by(ppi) %>%
  summarize(mean=mean(age), sd=sd(age))

(47.5-40.9)/(sqrt((16.6^2+14.8^2)/2))


# Tableone
library(tableone)  

vars<- c(factor_vars[2:8], "bmi", "n_consult")
SD_crude <- CreateTableOne(vars, data= dat, strata="ppi", test=FALSE)  
tbl_before <- data.frame(print(SD_crude, smd = TRUE))
tbl_before <- cbind(rownames(tbl_before), tbl_before)
colnames(tbl_before) <- (c('Variable','H2RA','PPI', 'ASD'))
tbl_before$Variable[tbl_before$Variable == 'n'] <- 'N'
tbl_before$Variable[tbl_before$Variable == 'gender = Male (%)'] <- 'Male (%)'
tbl_before$Variable[tbl_before$Variable == 'eth5 (%)'] <- 'Ethnicity (%)'
tbl_before$Variable[tbl_before$Variable == 'calendarperiod (%)'] <- 'Calendar period (%)'
tbl_before$Variable[tbl_before$Variable == 'prior_gast_cancer = 1 (%)'] <- 'Prior gastric cancer (%)'
tbl_before$Variable[tbl_before$Variable == 'recent_gerd = 1 (%)'] <- 'Recebt GERD (%)'
tbl_before$Variable[tbl_before$Variable == 'recent_peptic = 1 (%)'] <- 'Recent peptic ulcer (%)'
tbl_before$Variable[tbl_before$Variable == 'bmi (mean (SD))'] <- 'BMI (mean (SD))'
tbl_before$Variable[tbl_before$Variable == 'n_consult (mean (SD))'] <- 'No Consultations (mean (SD))'
tbl_final <- cbind(tbl_before, tbl_after)
colnames(tbl_final) <- (c('Variable','H2RA','PPI', 'ASD','H2RA.','PPI.', 'ASD.'))
gt(tbl_final)
# Cobalt package
  

covs <- subset(dat, select =vars ) 
bal.tab.unwt <- bal.tab(covs, 
        treat = dat$ppi, 
        binary="std", #mean differences for binary variables (i.e., difference in proportion) should be standardized
        stats = "mean.diffs",
        continuous="std", #mean differences for continuous variables should be standardized 
        s.d.denom="pooled",
        abs = T) 

bal.tab.unwt$Observations

#################################################################
##             Propensity score Weighting                      ##
##                                                             ##
#################################################################



dat$ATEwg <- ifelse(dat$ppi== 1, 1/dat$PS, 1/(1- dat$PS)) # as IPTW at Charire
dat$ATTwg <- ifelse(dat$ppi== 1, 1, dat$PS/(1- dat$PS)) 

dat %>%
  group_by(ppi) %>%
  summarize( min_ATEwg = min(ATEwg),  
             max_ATEwg = max(ATEwg), 
             mean_ATEwg = mean(ATEwg),  
             sd_ATEwg = sd(ATEwg),   
             median_ATEwg = median(ATEwg)) 

dat %>%   
  group_by(ppi) %>%
  summarize(min_ATTwg = min(ATTwg),  
            max_ATTwg = max(ATTwg), 
            mean_ATTwg = mean(ATTwg),  
            sd_ATTwg = sd(ATTwg),   
            median_ATTwg = median(ATTwg)) 

#-------------------------------------------------------
# Assessing covariate balance in the weighted sample
#-------------------------------------------------------
bal.tab(covs,
        treat = dat$ppi, 
        weights = dat$ATTwg, 
        binary="std", 
        continuous="std", 
        s.d.denom="pooled",
        abs = T)  


###########################################################################
###########################################################################
###                                                                     ###
###                         ESTIMATES                                   ###
###                                                                     ###
###########################################################################
###########################################################################

#################################################################
##             Cox regression                                  ##
##              plus PS                                        ##
#################################################################
library(survival)
cox.PS <- coxph(Surv(time, as.numeric(died)) ~ as.factor(ppi) + PS, data = dat)
broom::tidy(cox.PS, exponentiate=TRUE, conf.int=TRUE) 


# Assess PH assumption
sch.resid.cox.PS=cox.zph(cox.PS, transform = 'identity')
sch.resid.cox.PS
plot(sch.resid.cox.PS, col = "red")

#################################################################
##             Cox regression                                  ##
##              confounding control                            ##
#################################################################

cox.confounding <- coxph(Surv(time, as.numeric(died)) ~ as.factor(ppi)+ 
                               age + gender + eth5 + imd_person  +
                               prior_gast_cancer + recent_gerd +
                               recent_peptic + bmi + n_consult, data = dat)
broom::tidy(cox.confounding, exponentiate=TRUE, conf.int=TRUE) 

# Assess PH assumption
sch.resid.cox.confounding=cox.zph(cox.confounding, transform = 'identity')
sch.resid.cox.confounding
plot(sch.resid.cox.confounding, col = "red")
ggcoxzph(sch.resid.cox.confounding)

# ggcoxdiagnostics(cox.confounding, type = "schoenfeld", ox.scale = "time",
#                  title = "Diagnostic plot", subtitle = "Data comes from survey XYZ",
#                  font.subtitle = 9)
#################################################################
##             Using standarization                            ##
##                                                             ##
#################################################################

# Cohort where all got PPI
ppi_yes <- dat %>% 
  mutate(ppi= 1) 
# Cohort where all got H2RA
ppi_no <- dat %>% 
  mutate(ppi= 0) 

dat$ppi_1 <- cox.PS  %>% predict(ppi_yes, type = "expected") 
dat$ppi_0 <- cox.PS  %>% predict(ppi_no,  type = "expected")

dat %>%  
  summarize(N = n(),  
            mean_deaths_PPI = mean(ppi_1),  
            mean_deaths_H2RA = mean(ppi_0), 
            risk_ratio = mean_deaths_PPI/ mean_deaths_H2RA, 
            risk_diff = mean_deaths_PPI - mean_deaths_H2RA) 

#################################################################
##             Regression with Weighs                          ##
#################################################################

#-------------------------------
#. Normal package using weight
#-------------------------------
dat1 <- dat[, c("ppi", "age","gender","eth5","imd_person","calendarperiod",
                "prior_gast_cancer","recent_gerd","recent_peptic","bmi", "n_consult","time","died", "ATTwg")
                ]
cohort.cox.final2.w <- coxph(Surv(time, died) ~ ppi, weights = ATEwg, data = dat1)
broom::tidy(cohort.cox.final2.w, exponentiate = TRUE, conf.int = TRUE)  

# Assess PH assumption
bfit1.ph =cox.zph(cohort.cox.final2.w)
sch.resid.cohort.cox.final2.w
plot(bfit1.ph, var = 3,col = "red")
ggcoxzph(sch.resid.cohort.cox.final2.w)

if (pdfind) {  pdf(file = "figure2B.pdf", width = 5, height = 5) }
par(oma = c(2, 2, 0.5, 0.5), mar = c(2, 2, 0, 0))
plotcoxzph(x =  bfit1.ph[1], wd = 2)
mtext(side = 1, line = 2.5, text = "duration of therapy (days)", cex = 1.2)
mtext(side = 2, line = 2.2, text = expression(hat(beta)), cex = 1.2)
abline(a = 0, b = 0, lty = 3)
abline(lm(bfit1.ph$y[, 1] ~ bfit1.ph$x)$coefficients, lty = 3, col = "red", lwd = 2)
legend("bottomleft", legend = "bfb", bty = "n", inset = 0.08, cex = 1.5)


plot.cox.zph 
#-------------------------------
#. survey package using weight
#------------------------------- 
library(survey)
ATEwg_design <-svydesign(id=~1, weights=~ATEwg, data=dat)
#rpbc<-as.svrepdesign(ATTwg_design)

cox.ATE.rr <- svycoxph(Surv(time, as.numeric(died))~ as.factor(ppi), design=ATEwg_design)
broom::tidy(cox.ATE.rr, exponentiate=TRUE, conf.int=TRUE)

# glm.ATT.rr <- svyglm(died ~ ppi, family=quasibinomial(link="log"), design=ATTwg_design)
# broom::tidy(glm.ATT.rr, exponentiate=TRUE, conf.int=TRUE)

# Assess PH assumption
sch.resid.cox.ATT.rr =cox.zph(cox.ATT.rr, transform = 'identity')
sch.resid.cox.ATT.rr
plot(sch.resid.cox.ATT.rr, col = "red")
ggcoxzph(sch.resid.cox.ATT.rr)

#################################################################
##             Regression MAtchinf                         ##
#################################################################
library(MatchIt)
nnmatch <- matchit(ppi~PS, data = dat, distance="glm",
                   link="linear.logit", method= "nearest", ratio = 1 , caliper =0.1) 
summary(nnmatch, standardize=TRUE) 
ps_match_data <- match.data(nnmatch) 
tab_match <- CreateTableOne(vars, strata="ppi", ps_match_data, test=FALSE)  
print(tab_match, smd = TRUE) 

matched_model <- coxph(Surv(time, died) ~ ppi, data = ps_match_data)
broom::tidy(matched_model, exponentiate=TRUE, conf.int=TRUE) 

###########################################################################
###########################################################################
###                                                                     ###
###                         TABLE AFTER WEIGHTs                         ###
###                                                                     ###
###########################################################################
###########################################################################

## Matching weight

## Weighted data
library(survey)
library(tableone)
rhcSvy <- svydesign(ids = ~ 1, data = dat, weights = ~ ATTwg)

## Construct a table (This is a bit slow.)
tabWeighted <- svyCreateTableOne(vars = vars, strata = "ppi", data = rhcSvy, test = FALSE)
## Show table with SMD
tbl_after <- data.frame(print(tabWeighted, smd = TRUE))

### appling the ATTwg
bal.tab.wt <- bal.tab(covs, 
        treat = dat$ppi, 
        weights = dat$ATTwg, 
        #distance = dat$PS,
        binary="std", 
        continuous="std", 
        s.d.denom="pooled",
        abs = T)  

bal.tab.wt

#----------------------------------------
#. Table at Charite
#----------------------------------------
library(MatchIt) 
set.seed(4567) #set same seed for now, so we all get the same result!

#temp <- matchit(ppi~.,
       #         data= variables_model, caliper=.05, method="nearest") 
#temp <- matchit(ppi~.,
#                data= variables_model) 

# the nearest but with a caliper of .05. 
#Trade-off: how many people it is possible to match and how well the balance is. Check seminar.with a high caliper more matched
temp

matched.data <- dat[temp$weight==1,]
View(matched.data)

tablematched2 <- CreateTableOne(vars = vars, 
                                  strata = "ppi", 
                                  data = matched.data, 
                                  test=FALSE)
print(tablematched2, smd = TRUE)

covs2 <- subset(matched.data, select =vars ) 
bal.tab(covs2, treat = matched.data$ppi, stats = "mean.diffs", binary="std", continuous="std", s.d.denom="pooled")$Observations





#################################################################
##             IPTW       at charite                       ##
##                   and plots                                   ##
#################################################################

# same ATE

dat <- within(dat, {
  IPTW <- NA
  IPTW[ppi == "1"] <- (1/dat$PS)[ppi == "1"]
  IPTW[ppi == "0"] <- (1/(1-dat$PS))[ppi == "0"]
})

#---------
#.  plots
#---------
library(ggplot2)
## Unweighted cohort
ggplot(data = dat,
       mapping = aes(x = PS, fill = factor(ppi), weight = NULL)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  theme(legend.key = element_blank()) +
  labs(title = "Unweighted cohort")

#-------
#. after
#--------

## IPTW the distributions almost overlap, this is the result by weighting ¡, moving confounding
ggplot(data = dat,
       mapping = aes(x = PS, fill = factor(ppi), weight = ATTwg)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  theme(legend.key = element_blank()) + 
  labs(title = "IPTW")




