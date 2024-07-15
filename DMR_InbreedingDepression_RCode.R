###########################################################

# R Code supporting paper "Evidence for inbreeding depression in captive Damaraland mole-rats"
# Authors: Dave Seager and Jack Thorley
# Email: ds993@cam.ac.uk

###########################################################

# The analysis consists of 12 statistical models separated into two parts 
# Part 1 covers the effects on litter size, inter-birth interval (IBI), birth weight, and survival 
# Part 2 covers the effects on skeletal and body mass growth. 

# The models consist of the following:
# Part 1: 
  # 1- Litter Size ~ treatment 
  # 2- Interbirth interval (IBI) ~ treatment
  # 3- Birth weight ~ treatment 
  # 4- Early life survival ~ treatment   # (survival to 30 days)
  # 5- Early life survival ~ treatment*group size 
  # 6- Longer-term survival ~ treatment

# Part 2: 
  # 7 - Male body mass growth ~ treatment    ("bodymass.nlme.male")
  # 8 - Female body mass growth ~ treatment  ("bodymass.nlme.female")
  # 9 - Male body length growth ~ treatment    ("bodylength.nlme.male")
  # 10 - Female body length growth ~ treatment   ("bodylength.nlme.female")
  # 11 - Male teeth width growth ~ treatment   ("teethwidth.nlme.male")
  # 12 - Female teeth width growth ~ treatment   ("teethwidth.nlme.male")

# In all cases, treatment refers to whether the offspring were inbred or outbred

#-----------------------------------------------------------

# set working directory to the location of the data sets
setwd("INSERT FILE PATH")

# load required packages
lapply(c("coxme", "DHARMa", "glmmTMB", "ggplot2", "MASS", "nlme", "patchwork", "survminer", "tidyverse"), FUN = library, character.only = TRUE)

# load the required data sets
LitterSize_Dataset <- read.csv("LitterSize_Dataset.csv", header = T)
InterBirthInterval_Dataset <- read.csv("InterBirthInterval_Dataset.csv", header = T)
Pup_Weight_Dataset <- read.csv("PupWeights_Dataset.csv", header = T)
Survival_Dataset <- read.csv("Survival_Dataset.csv", header = T)
IBWeights <- read.csv("LifetimeWeights_Dataset.csv", header = T)
Skeletal_Dataset <-  read.csv("LifetimeSkeletal_Dataset.csv", header = T)

#=============================

# PART 1: LITTER SIZE, IBI, BIRTH WEIGHT AND EARLY SURVIVAL

#============================

# Model 1: Modelling Litter Size ~ treatment (n = 109)

# To model litter size we used a dataset of 109 litters for which litter size is known (n = 16 females).
## The response variable is litter size
### The fixed effects tested were treatment (inbred/outbred), group Size (number of individuals in group) and density (number of animals per metre of tunnel)
#### maternal ID (mother of the litter) is included as a random effect

model1 <- glmmTMB(LitterSize ~ Treatment + GroupSize + Density + (1|MaternalID), 
                  data = LitterSize_Dataset)
summary(model1)
residuals1 <- simulateResiduals(model1)
plot(residuals1)
# None of the fixed effects predict litter size 

#---------------------------------------------------------------------#

# Model 2: Modelling Inter-birth interval ~ treatment (n = 118)

# To model inter-birth Interval we used a dataset of 118 litters
## 1 litter removed due to unusually long time between litters (818 days)
## The response variable is inter-birth interval in days - gestation period (90 days)
### The fixed effects tested were treatment, group size (number of individuals in group) and density (number of animals per metre of tunnel)
#### maternal ID is included as a random effect

model2 <- glmmTMB(InterBirth_rel ~ Treatment + GroupSize + Density + (1|MaternalID),
                data = InterBirthInterval_Dataset, 
                family = Gamma(link = "log"))
summary(model2)
residuals2 <- simulateResiduals(model2)
plot(residuals2) 
# None of the fixed effects predict inter-birth interval but there is a trend for larger colonies to have a shorter interbirth interval 

#---------------------------------------------------------------------#

# Model 3: Modelling birth weight ~ treatment (n = 292)

# To model offspring mass at birth and its effect on survival we used a dataset of 292 individuals with a weight taken within 5 days of birth (n= 16 mothers, n = 97 litters)
## The response variable was offspring weight at birth in grams
### The fixed effects were Treatment, group size, Litter size, density and the difference in days between birth and the weight being taken (maximum 5)
#### maternal ID and Litter Reference are included as random effects

model3 <- glmmTMB(Weight ~ Treatment + GroupSize + LitterSize + Density + Wt_Diff + (1|Mother) + (1|LitterRef),
                data = Pup_Weight_Dataset)
summary(model3)
residuals3 <- simulateResiduals(model3)
plot(residuals3)

# Individuals from outbred colonies are heavier than those from inbred colonies
## Individuals in larger colonies are lighter
### Individuals in larger litters are lighter
#### Weight increases with the number of days that animals are weighed after birth

# plot birth weight ~ treatment to produce plot 1a
newdf <- data.frame(Treatment = unique(Pup_Weight_Dataset$Treatment), 
                    GroupSize = mean(Pup_Weight_Dataset$GroupSize), 
                    LitterSize = mean(Pup_Weight_Dataset$LitterSize), 
                    Density = mean(Pup_Weight_Dataset$Density), 
                    Wt_Diff = mean(Pup_Weight_Dataset$Wt_Diff), 
                    Mother = unique(Pup_Weight_Dataset$Mother)[1], 
                    LitterRef = unique(Pup_Weight_Dataset$LitterRef)[1])

newdf$pred <- predict(model3, newdata = newdf, re.form = NA)
newdf$SE <- predict(model3, newdata = newdf, re.form = NA, se.fit = TRUE)$se.fit
newdf$l95 <-  newdf$pred - 1.96*newdf$SE
newdf$u95 <-  newdf$pred + 1.96*newdf$SE

p1 <- ggplot(Pup_Weight_Dataset,aes(x=Treatment, y = Weight, colour = Treatment)) +
  geom_jitter(aes(x=Treatment,y=Weight), height=0.05,width=0.1,alpha=0.4, shape = 1) +
  geom_errorbar(data = newdf,  aes(y = pred, ymin=l95, ymax = u95), 
                colour = "black", width=0, linewidth = 1) +
  geom_point(data = newdf,  aes(y = pred), size = 3.5) +
  geom_point(data = newdf,  aes(y = pred), size = 3.5, shape = 1, colour = "black") +
  labs(x = "Treatment", y = "Pup weight (g)", tag = "(a)") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none", 
        plot.tag = element_text(size = 13)) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(4, 14, 1), 
                     labels = c("4", "", "6", "", "8", "", "10", "", "12", "", "14"))

#---------------------------------------------------------------------#
# Model 4: Modelling early life survival ~ birthweight (n = 292)

# To model the effect of offspring mass at birth on survival we used a dataset of 292 individuals with a weight taken within 5 days of birth (n = 16 monthers, n = 97 litters)
## The response variable was individual survival to 30 days
### The fixed effects were Treatment, group size, weight at birth  and density
#### maternal ID and Litter Reference are included as random effects

# To test how birthweight affects survival to 30 days
model4 <- glmmTMB(Sur30 ~ Treatment + GroupSize + Weight + Density + 
                  (1|Mother) + (1|LitterRef),
                data=Pup_Weight_Dataset, 
                family = "binomial") 
summary(model4)
residuals4 = simulateResiduals(model4)
plot(residuals4)

# Offspring weight predicts survival to 30 days with lighter individuals more likely to die early

# plot survival to 30 days ~ birth weight to produce plot 1b
newdf2 <- data.frame(Treatment = unique(Pup_Weight_Dataset$Treatment), 
                     GroupSize = mean(Pup_Weight_Dataset$GroupSize), 
                     Weight = seq(5, 14, 0.2), 
                     Density = mean(Pup_Weight_Dataset$Density), 
                     Mother = unique(Pup_Weight_Dataset$Mother)[1], 
                     LitterRef = unique(Pup_Weight_Dataset$LitterRef)[1])

# restrict the plot over the weight range of each treatment 
range(Pup_Weight_Dataset$Weight[Pup_Weight_Dataset$Treatment == "Inbred" & !is.na(Pup_Weight_Dataset$Sur30)])
range(Pup_Weight_Dataset$Weight[Pup_Weight_Dataset$Treatment == "Outbred" & !is.na(Pup_Weight_Dataset$Sur30)])
newdf2 <- newdf2 %>% 
  filter(!(Treatment == "Inbred" & !between(Weight, 5, 13))) %>% 
  filter(!(Treatment == "Outbred" & !between(Weight, 6, 14))) 

newdf2$pred <- predict(model4, newdata = newdf2, re.form = NA, type = "link")
newdf2$SE <- predict(model4, newdata = newdf2, re.form = NA, se.fit = TRUE, type = "link")$se.fit
inv_logit <- function(x) { exp(x)/(1+exp(x)) }
newdf2$l95 <-  inv_logit(newdf2$pred - 1.96*newdf2$SE)
newdf2$u95 <-  inv_logit(newdf2$pred + 1.96*newdf2$SE)
newdf2$pred <- inv_logit(newdf2$pred)

p2 <- ggplot(Pup_Weight_Dataset,aes(x=Weight, group = Treatment, colour = Treatment)) +
  geom_jitter(aes(y = Sur30), height=0.05, width = 0.1, alpha= 0.4, shape = 1) +
  geom_ribbon(data = newdf2, aes(ymin = l95, ymax = u95, fill = Treatment), 
              colour = NA, alpha = 0.15) +
  geom_line(data = newdf2, aes(y = pred), linewidth = 1) +
  labs(x = "Pup weight (g)", y = "Survival to 30 days", tag = "(b)") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none", 
        plot.tag = element_text(size = 13)) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_x_continuous(breaks = seq(6, 14, 1),
                     labels = c("6", "", "8", "", "10", "", "12", "", "14"))

#---------------------------------------------------------------

# Model 5: Model the interaction between status and group size to survival at 30 days (n = 328)

# To model the interaction between group size and survival to 30 days we used a dataset of 328 individuals 
## The response variable was survival to 30 days
### The fixed effects were Treatment, group size and density
#### maternal ID (n = 16) and Litter Reference (n = 108) are included as random effects
Survival_Dataset$Treatment
model5 <- glmmTMB(Sur30 ~ Density + Treatment * GroupSize + (1|Mother) + (1|LitterRef),
                 data = Survival_Dataset,
                 family = "binomial") 
summary(model5)
residuals5 <- simulateResiduals(model5)
plot(residuals5)

# There is a significant interaction between group size and survival to 30 days with inbred individuals in smaller colonies being less likely to survive

# plot the interaction effect of survival to 30 days ~ group size * treatment to produce plot 1c
newdf3 <- expand.grid(Treatment = unique(Survival_Dataset$Treatment), 
                      GroupSize = seq(2, 26, 1), 
                      Density = mean(Survival_Dataset$Density), 
                      Mother = unique(Survival_Dataset$Mother)[1], 
                      LitterRef = unique(Survival_Dataset$LitterRef)[1])

# restrict the plot over the weight range of each treatment 
range(Survival_Dataset$GroupSize[Survival_Dataset$Treatment == "Inbred" & !is.na(Survival_Dataset$Sur30)])
range(Survival_Dataset$GroupSize[Survival_Dataset$Treatment == "Outbred" & !is.na(Survival_Dataset$Sur30)])
newdf3 <- newdf3 %>% 
  filter(!(Treatment == "Inbred" & !between(GroupSize, 2, 26))) %>% 
  filter(!(Treatment == "Outbred" & !between(GroupSize, 2, 24))) 

newdf3$pred <- predict(model5, newdata = newdf3, re.form = NA, type = "link")
newdf3$SE <- predict(model5, newdata = newdf3, re.form = NA, se.fit = TRUE, type = "link")$se.fit
newdf3$l95 <-  inv_logit(newdf3$pred - 1.96*newdf3$SE)
newdf3$u95 <-  inv_logit(newdf3$pred + 1.96*newdf3$SE)
newdf3$pred <-  inv_logit(newdf3$pred)

p3 <- ggplot(Survival_Dataset,aes(x = GroupSize, group = Treatment, 
                                  colour = Treatment)) +
  geom_jitter(aes(y = Sur30), height=0.05, width = 0.1, alpha=0.4, shape = 1) +
  geom_ribbon(data = newdf3, aes(ymin = l95, ymax = u95, fill = Treatment), 
              colour = NA, alpha = 0.15) +
  geom_line(data = newdf3, aes(y = pred), linewidth = 1) +
  labs(x = "Group size", y = "Survival to 30 days", tag = "(c)") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none", 
        plot.tag = element_text(size = 13)) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_x_continuous(breaks = seq(2, 26, 4))

#---------------------------------------------------------------

# Model 6: Use Cox models to model the difference in survival between treatments

# To model the difference in individual survival between treatments  we used a dataset of 328 individuals 
## The fixed effects were Treatment and density
#### maternal ID and Litter Reference are included as random effects

model6 <- coxme(Surv(Mnths,Censorship_D) ~ Treatment + Density + GroupSize + 
                  (1|LitterRef) + (1|Mother), 
                data = Survival_Dataset)
summary(model6)

# Inbred individuals have significantly lower survival than outbred individuals

# plot the output from the cox survival model 1d
p4 <- ggsurvplot(
  fit = survfit(Surv(Mnths, Censorship_D) ~ Treatment + (1|LitterRef) + (1|Mother), 
                data = Survival_Dataset), 
  xlab = "Months", 
  ylab = "Survival",
  legend.labs = c("Inbred Colonies", "Outbred Colonies"),
  conf.int = TRUE)

p4 <- p4[[1]]
p4 <- p4 + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none", 
        plot.tag = element_text(size = 13)) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  labs(tag = "(d)") + 
  scale_x_continuous(breaks = seq(0, 48, 12), limits = c(0, 50)) # + 
  #annotate("text", x = 36, y = 0.15, label = "Inbred", colour = "indianred") + 
  #annotate("text", x = 36, y = 0.25, label = "Outbred", colour = "dodgerblue")

# Combine the plots together with patchwork
(p1 | p2) / (p3 | p4)


#=============================

# PART 2: GROWTH MODELLING

#=============================

# Overall the body weights data set consists of 9126 weights from 92 individuals (40 females and 52 males)  
# Before modelling body mass growth, we need to separate the two sexes. 

# Due to the high degree of sexual dimorphism the sexes should be modelled separately
Male   <- filter(IBWeights, Sex == "M") %>% 
  arrange(AnimalID, Age)
Female <- filter(IBWeights, Sex == "F") %>% 
  arrange(AnimalID, Age)

# To check the number of males and observations and provide a median and SD
length(unique(Male$AnimalID)) # n = 52
nrow(Male)
Male %>% 
  group_by(AnimalID) %>% 
  summarise(n = n())  %>% 
  summarise(median(n), sd(n)) %>% 
  data.frame()

# To check the number of females and observations and provide a median and SD
length(unique(Female$AnimalID)) # n = 40
nrow(Female)
Female %>% 
  group_by(AnimalID) %>% 
  summarise(n = n())  %>% 
  summarise(median(n), sd(n)) %>% 
  data.frame()

# We also need to include inbreeding as a dummy variable coded as 0/1
Male$Inbred <- if_else(Male$Treatment == "IB", 1, 0)
Female$Inbred <- if_else(Female$Treatment == "IB", 1, 0)

# Plot the raw data to see the spread of weights ~ age
# For males
ggplot(Male, aes(x = Age, y = Weight, 
                 group = as.factor(Treatment), 
                 col = as.factor(Treatment))) + 
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point() +
  geom_smooth(method = "gam", colour = "black", 
              formula = y ~ s(x, k = 6, bs = "tp")) + 
  facet_wrap(~AnimalID) +
  xlim(c(0, 1250)) +
  ylim(c(0, 250)) +
  theme_bw() + 
  theme(panel.grid = element_blank())

## And for females
ggplot(Female, aes(x = Age, y = Weight, 
                   group = as.factor(Treatment), 
                   col = as.factor(Treatment))) + 
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point() +
  geom_smooth(method = "gam", colour = "black", 
              formula = y ~ s(x, k = 6, bs = "tp")) + 
  facet_wrap(~AnimalID) +
  xlim(c(0, 1250)) +
  ylim(c(0, 250)) +
  theme_bw() + 
  theme(panel.grid = element_blank())


#---------------------------------------------------------------

## Modelling body mass growth with a monomolecular model 

#----------------------------------------------------------------

# For body mass, we use a biphasic formulation of a monomolecular growth curve,
# where the two phases of growth occur either side of 50 days (~ age of weaning)

## First fit the model without random effects just to get an idea of good starting parameters
## i.e. fit as NLS

# The biphasic model formula: 
biphasic_inbred <- function(age,t0, Asym, AsymI, k1, k1I, k2, k2I, t1 = 50, Inbred) {
  
  val1 <- (Asym + AsymI*Inbred)*(1 - exp((k1 + k1I*Inbred)*(t0 - age)))
  val2 <- (Asym + AsymI*Inbred)*(1 - exp((k1 + k1I*Inbred)*(t0 - t1) + (k2 + k2I*Inbred)*(t1 - age)))
  val1[age > t1] <- val2[age > t1]
  
  return(val1)
}

# Male body mass NLS
bodymass.nls.male <- nls(Weight ~ biphasic_inbred(age=Age,t0=t0,Asym=Asym, AsymI=AsymI, 
                                                  k1=k1,k1I=k1I,k2=k2, k2I=k2I, Inbred = Inbred),
                         data = Male,
                         start = c(Asym = 150, AsymI = -30, k1 = 0.001, k1I = 0, k2 = 0.001, k2I = 0, 
                                   t0 = 0),
                         na.action = na.pass) 
summary(bodymass.nls.male)   

# Female body mass NLS
bodymass.nls.female <-nls(Weight ~ biphasic_inbred(age=Age,t0=t0,Asym=Asym, AsymI=AsymI, 
                                                   k1=k1,k1I=k1I,k2=k2, k2I=k2I, Inbred = Inbred),
                          data = Female,
                          start = c(Asym = 150, AsymI = -30, k1 = 0.001, 
                                    k1I = 0.0001, k2 = 0.001, k2I = 0.0001, 
                                    t0 = 0),
                          na.action = na.pass) 
summary(bodymass.nls.female) 

# Fit as NLMM
# Male body mass NLMM (model 7)
bodymass.nlme.male <- nlme(Weight ~ biphasic_inbred(age=Age,Asym=Asym,AsymI=AsymI, 
                                                    k1=k1,k1I=k1I,k2=k2, k2I=k2I, t0=t0, Inbred = Inbred),
                           data = Male,
                           fixed = (Asym+AsymI+k1+k1I+k2+k2I+t0~1),
                           random = list(AnimalID = pdDiag(Asym+k1+k2 +t0~ 1)),
                           correlation = corAR1(form = ~ 1|AnimalID),
                           start = c(Asym = 250, AsymI = -40, k1 = 0.002, k1I = 0, 
                                     k2 = 0.001, k2I = 0, t0 = 0),
                           na.action = na.pass)
summary(bodymass.nlme.male) 
plot(bodymass.nlme.male, resid(.) ~ Age) 

# Female body mass NLMM (model 8)
bodymass.nlme.female <- nlme(Weight ~ biphasic_inbred(age=Age,t0=t0,Asym=Asym, AsymI=AsymI, 
                                                      k1=k1,k1I=k1I,k2=k2, 
                                                      k2I=k2I, Inbred = Inbred),
                             data = Female,
                             fixed = (Asym+AsymI+k1+k1I+k2+k2I+t0~1),
                             random = list(AnimalID = pdDiag(Asym+k1+k2 + t0~ 1)),
                             correlation = corAR1(form = ~ 1|AnimalID),
                             weights = varPower(1),
                             start = c(Asym = 140, AsymI = -30, k1 = 0.001, k1I = 0, 
                                       k2 = 0.001, k2I = 0,  t0 = 0),
                             na.action = na.pass)
summary(bodymass.nlme.female) 
plot(bodymass.nlme.female, resid(.) ~ Age) 

# Predict and plot the best fitting model in either case: 
# range(Male$Age) ; range(Female$Age)
newdf <- expand_grid(Age = 1:1200, Inbred = c(0, 1))

male_pop <- bind_cols(newdf, PredWeight = predict(bodymass.nlme.male, newdata = newdf, level = 0)) %>% 
  mutate(Sex = "Male")
female_pop <- bind_cols(newdf, PredWeight = predict(bodymass.nlme.female, newdata = newdf, level = 0)) %>% 
  mutate(Sex = "Female")

# Get the confidence intervals (population prediction intervals in the parlance of: 
# "https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme)"

# now get the bootstrapped values for the mean
nresamp <- 1000

## pick new parameter values by sampling from multivariate normal distribution based on fit
pars_male <- mvrnorm(nresamp, 
                     mu = fixef(bodymass.nlme.male), 
                     Sigma = vcov(bodymass.nlme.male))
pars_female <- mvrnorm(nresamp, 
                       mu = fixef(bodymass.nlme.female), 
                       Sigma = vcov(bodymass.nlme.female))

# utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}

set.seed(101)
yvals1 <- matrix(nrow = length(male_pop$Age), ncol = 1000)
yvals2 <- matrix(nrow = length(female_pop$Age), ncol = 1000)
# 1000 cols for each Age. So I will predict for 1000 draws from the param distribution 
# (assuming multivariate normal dist)

biphasic_inbred <- function(age,t0, Asym, AsymI, k1, k1I, k2, k2I, t1 = 50, Inbred) {
  
  val1 <- (Asym + AsymI*Inbred)*(1 - exp((k1 + k1I*Inbred)*(t0 - age)))
  val2 <- (Asym + AsymI*Inbred)*(1 - exp((k1 + k1I*Inbred)*(t0 - t1) + (k2 + k2I*Inbred)*(t1 - age)))
  val1[age > t1] <- val2[age > t1]
  
  return(val1)
}

# male CI prediction
## There are a thousand parameter draws. 
### Each will be predicted against each age (to generate a 95% CI for each age from 0 to 1200 days)
for (i in 1:length(male_pop$Age)) {
  print(i)
  yvals1[i,] <- biphasic_inbred(age = male_pop$Age[i], 
                                Inbred = male_pop$Inbred[i],
                                Asym = pars_male[,1], 
                                AsymI = pars_male[,2], 
                                k1 = pars_male[,3], 
                                k1I = pars_male[,4], 
                                k2 = pars_male[,5], 
                                k2I = pars_male[,6], 
                                t0 = pars_male[,7])
  
  
}

ci.males <- get_CI(yvals1)
male_pop <- cbind(male_pop, ci.males)  

# Female prediction
for (i in 1:length(female_pop$Age)) {
  print(i)
  yvals2[i,] <- biphasic_inbred(age = female_pop$Age[i], 
                                Inbred = female_pop$Inbred[i],
                                Asym = pars_female[,1], 
                                AsymI = pars_female[,2], 
                                k1 = pars_female[,3], 
                                k1I = pars_female[,4], 
                                k2 = pars_female[,5], 
                                k2I = pars_female[,6], 
                                t0 = pars_female[,7])
}

ci.females <- get_CI(yvals2)
female_pop <- cbind(female_pop, ci.females)  

# Join the two sexes in a single data frame: 
twosex_pop <- bind_rows(male_pop, female_pop) %>% 
  mutate(Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age/365.35)

# Location of the 95% aysymptotic mass
male_outbred_95 <- 0.95*220.21132                  # male outbred
male_inbred_95  <- 0.95*(220.21132   - 50.46211)  # male inbred
female_outbred_95 <- 0.95*151.11225                # female outbred
female_inbred_95  <- 0.95*(151.11225 - 36.06486)   # female inbred 

male_outbred_95_age <- filter(twosex_pop, Sex == "Male" & Inbred == "Outbred") %>% 
  mutate(diffweight = abs(PredWeight - male_outbred_95)) %>% 
  filter(diffweight == min(diffweight)) %>% 
  .$Age # > 1200 days

male_inbred_95_age <- filter(twosex_pop, Sex == "Male" & Inbred == "Inbred") %>% 
  mutate(diffweight = abs(PredWeight - male_inbred_95)) %>% 
  filter(diffweight == min(diffweight)) %>% 
  .$Age # > 1200 days 

female_outbred_95_age <- filter(twosex_pop, Sex == "Female" & Inbred == "Outbred") %>% 
  mutate(diffweight = abs(PredWeight - female_outbred_95)) %>% 
  filter(diffweight == min(diffweight)) %>% 
  .$Age # > 1200 days

female_inbred_95_age <- filter(twosex_pop, Sex == "Female" & Inbred == "Inbred") %>% 
  mutate(diffweight = abs(PredWeight - female_inbred_95)) %>% 
  filter(diffweight == min(diffweight)) %>% 
  .$Age # 930 days

# The body mass plot with both
p_bodymass <- ggplot(twosex_pop, aes(x = Age, y = PredWeight, 
                                       group = as.factor(Inbred), 
                                       colour = as.factor(Inbred), 
                                       fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") +  
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = c(0.12, 0.88), 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        legend.background = element_blank()) + 
  labs(y = "Body mass (g)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) 

# Estimate the growth rate at each phase of life (g/day)
# Can do this easily by taking the successive deltas
twosex_pop <- twosex_pop %>% 
  group_by(Sex, Inbred) %>% 
  mutate(Delta = PredWeight - lag(PredWeight)) %>% 
  data.frame()

# To plot this we first plot by sex and then grob in the inset IGR
p_bodymass_males <- ggplot(filter(twosex_pop, Sex == "Male"),
                             
                             aes(x = Age, y = PredWeight, 
                                 group = as.factor(Inbred), 
                                 colour = as.factor(Inbred), 
                                 fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") +  
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none", 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        legend.background = element_blank()) + 
  labs(y = "Body mass (g)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 210))

p_malebmgrowthrate <- ggplot(data = filter(twosex_pop, Age <= 50 & Sex == "Male"), 
                           aes(x = Age,
                               y = Delta, 
                               group = as.factor(Inbred), 
                               colour = as.factor(Inbred), 
                               fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") + 
  geom_path(linewidth = 0.6) + 
  geom_path(data = filter(twosex_pop, Age > 50 & Sex == "Male"), linewidth = 0.6) + 
  facet_wrap(~ Sex) + 
  theme_classic() + 
  theme(axis.text.x = element_blank() , 
        axis.text.y = element_text(size = 8, colour = "black"), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 0.5))

p_bodymass_males <- p_bodymass_males + 
  annotation_custom(ggplotGrob(p_malebmgrowthrate), 
                    xmin = 600, xmax = 1200, 
                    ymin = 0, ymax = 90)


p_bodymass_females <- ggplot(filter(twosex_pop, Sex == "Female"),
                               
                               aes(x = Age, y = PredWeight, 
                                   group = as.factor(Inbred), 
                                   colour = as.factor(Inbred), 
                                   fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") +  
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = c(0.3, 0.88), 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        legend.background = element_blank()) + 
  labs(y = "Body mass (g)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 210))

p_femalebmgrowthrate <- ggplot(data = filter(twosex_pop, Age <= 50 & Sex == "Female"), 
                             aes(x = Age,
                                 y = Delta, 
                                 group = as.factor(Inbred), 
                                 colour = as.factor(Inbred), 
                                 fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") + 
  geom_path(linewidth = 0.6) + 
  geom_path(data = filter(twosex_pop, Age > 50 & Sex == "Female"), linewidth = 0.6) + 
  facet_wrap(~ Sex) + 
  theme_classic() + 
  theme(axis.text.x = element_blank() , 
        axis.text.y = element_text(size = 8, colour = "black"), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 0.5))

p_bodymass_females <- p_bodymass_females + 
  annotation_custom(ggplotGrob(p_femalebmgrowthrate), 
                    xmin = 600, xmax = 1200, 
                    ymin = 0, ymax = 90)

# The final plot for body mass
p_bodymass2 <- (p_bodymass_females + 
                    theme(plot.margin = margin(r = 0), 
                          axis.title.x = element_blank(), 
                          axis.text.x = element_blank())) | (p_bodymass_males + 
                                                               theme(axis.text.y = element_blank(), 
                                                                     axis.ticks.y = element_blank(), 
                                                                     axis.title.y = element_blank(), 
                                                                     plot.margin = margin(l = 0), 
                                                                     axis.title.x = element_blank(), 
                                                                     axis.text.x = element_blank()))

#--------------------------------
### Generate body mass predictions for each individual to see how the model captures individual variation

newdf_male <- Male %>% 
  group_by(AnimalID, Inbred) %>%
  tidyr::expand(Age = seq(1, max(Age), 1)) %>% 
  data.frame()

newdf_female <- Female %>% 
  group_by(AnimalID, Inbred) %>%
  tidyr::expand(Age = seq(1, max(Age), 1)) %>% 
  data.frame()

male_ind <- bind_cols(newdf_male, PredWeight = predict(bodymass.nlme.male, 
                                                       newdata = newdf_male)) %>% 
  mutate(Sex = "Male", 
         Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age/365.35)

female_ind <- bind_cols(newdf_female, PredWeight = predict(bodymass.nlme.female, 
                                                           newdata = newdf_female)) %>% 
  mutate(Sex = "Female", 
         Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age/365.35)

# Male individual curves
Male_toplot <- mutate(Male, Inbred = if_else(Inbred == 0, "Outbred", "Inbred"))

p_bodymass_maleind <- ggplot(Male_toplot, aes(x = Age, y = Weight, group = Inbred, col = Inbred)) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point(alpha = 0.1) +
  geom_line(data = male_ind, aes(x = Age, y = PredWeight), 
            colour = "black", linewidth = 0.8) + 
  facet_wrap(~AnimalID, ncol = 7) +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(colour = "black", angle= 90, size = 9), 
        axis.text.y = element_text(colour = "black"), 
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x = "Age (days)", y = "Body mass (g)", 
       subtitle = "Male growth curves", 
       tag = "B") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Female individual curves
Female_toplot <- mutate(Female, Inbred = if_else(Inbred == 0, "Outbred", "Inbred"))

p_bodymass_femaleind <- ggplot(Female_toplot, aes(x = Age, y = Weight, group = Inbred, col = Inbred)) +
  geom_hline(yintercept = 10, linetype = 2) + 
  geom_point(alpha = 0.1) +
  geom_line(data = female_ind, aes(x = Age, y = PredWeight), 
            colour = "black", linewidth = 0.8) + 
  facet_wrap(~AnimalID) +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(colour = "black", angle= 90, size = 9), 
        axis.text.y = element_text(colour = "black"), 
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x = "Age (days)", y = "Body mass (g)", 
       subtitle = "Female growth curves", 
       tag = "A") + 
  scale_colour_manual(values = c("indianred", "dodgerblue"))  + 
  scale_x_continuous(breaks = seq(0, 1250, 250))

p_bodymassindcurves <- (p_bodymass_femaleind / p_bodymass_maleind) + 
  plot_layout(heights = c(6,8.5))

#---------------------------------------------------------------

## Modelling body length and teeth width a monomolecular model 

#---------------------------------------------------------------

# Here we use a monophasic form of the monomolecular model as specified directly in nls/nlme formula. 

# As for body mass we want to model males and females separately
MaleSkel <- filter(Skeletal_Dataset, Sex == "M") %>% 
  arrange(AnimalID, Age_days)
FemaleSkel <- filter(Skeletal_Dataset, Sex == "F") %>% 
  arrange(AnimalID, Age_days) %>% 
  mutate(UpTeethWidth = if_else(UpTeethWidth < 2, NA, UpTeethWidth))  # obvious outlier 

# and get the number of individuals and repeats/individual
length(unique(MaleSkel$AnimalID)) # n = 52
nrow(MaleSkel)
MaleSkel %>% 
  group_by(AnimalID) %>% 
  summarise(n = n())  %>% 
  summarise(median(n), sd(n)) %>% 
  data.frame()

length(unique(FemaleSkel$AnimalID)) # n = 39
nrow(FemaleSkel)
FemaleSkel %>% 
  group_by(AnimalID) %>% 
  summarise(n = n())  %>% 
  summarise(median(n), sd(n)) %>% 
  data.frame()

# I will include inbreeding as a dummy var coded as 0/1
MaleSkel$Inbred <- if_else(MaleSkel$Treatment == "IB", 1, 0)
FemaleSkel$Inbred <- if_else(FemaleSkel$Treatment == "IB", 1, 0)

# Plot the raw data to see the spread of skeletal measures ~ age
# Male body length (all individuals)
ggplot(MaleSkel, aes(x = Age_days, y = BodyLength, 
                     group = as.factor(Treatment), 
                     col = as.factor(Treatment))) + 
  geom_point() +
  geom_smooth() + 
  labs(x = "Age (days)", y = "Body Length (cm)", colour = "Inbreeding", 
       title = "Male body length") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Female body length (all individuals)
ggplot(FemaleSkel, aes(x = Age_days, y = BodyLength, 
                                         group = as.factor(Treatment), 
                                         col = as.factor(Treatment))) + 
  geom_point() +
  geom_smooth() + 
  labs(x = "Age (days)", y = "Body Length (cm)", colour = "Inbreeding", 
       title = "Female body length") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Male teeth width (all individuals)
ggplot(MaleSkel, aes(x = Age_days, y = UpTeethWidth, 
                     group = as.factor(Treatment), 
                     col = as.factor(Treatment))) + 
  geom_point() +
  geom_smooth() + 
  labs(x = "Age (days)", y = "Teeth width (mm)", colour = "Inbreeding", 
       title = "Male teeth width") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Female teeth width (all individuals)
ggplot(FemaleSkel, aes(x = Age_days, y = UpTeethWidth, 
                       group = as.factor(Treatment), 
                       col = as.factor(Treatment))) + 
  geom_point() +
  geom_smooth() + 
  labs(x = "Age (days)", y = "Teeth width (mm)", colour = "Inbreeding", 
       title = "Female teeth width") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Also plot the raw data within individuals to get an idea of the spread 
# Male body length (by individual)
ggplot(MaleSkel, aes(x = Age_days, y = BodyLength, 
           group = as.factor(Treatment), 
           col = as.factor(Treatment))) + 
  geom_point() +
  geom_line(colour = "black") + 
  labs(x = "Age (days)", y = "Body Length (cm)", 
       colour = "Inbreeding", title = "Male growth curves") + 
  facet_wrap(~AnimalID, ncol = 6) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Female body length (by individual)
ggplot(FemaleSkel, aes(x = Age_days, y = BodyLength, 
                       group = as.factor(Treatment), 
                       col = as.factor(Treatment))) + 
  geom_point() +
  geom_line(colour = "black") + 
  labs(x = "Age (days)", y = "Body Length (cm)", 
       colour = "Inbreeding", title = "Female growth curves", tag = "A") + 
  facet_wrap(~AnimalID, ncol = 6) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250))

# Male teeth width (by individual)
ggplot(MaleSkel, aes(x = Age_days, y = UpTeethWidth, 
                                     group = as.factor(Treatment), 
                                     col = as.factor(Treatment))) + 
  geom_point() +
  geom_line(colour = "black") + 
  labs(x = "Age (days)", y = "Incisor width (mm)", 
       colour = "Inbreeding", title = "Male growth curves") + 
  facet_wrap(~AnimalID, ncol = 6) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250)) + 
  scale_y_continuous(breaks = seq(2, 8, 2))

# Female teeth width (by individual)
ggplot(FemaleSkel, aes(x = Age_days, y = UpTeethWidth, 
                                         group = as.factor(Treatment), 
                                         col = as.factor(Treatment))) + 
  geom_point() +
  geom_line(colour = "black") + 
  labs(x = "Age (days)", y = "Incisor width (mm)", 
       colour = "Inbreeding", title = "Female growth curves") + 
  facet_wrap(~AnimalID, ncol = 6) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 1250, 250)) + 
  scale_y_continuous(breaks = seq(2, 8, 2))

#---------------------------------------------------------------

## Modelling body length growth 

#----------------------------------------------------------------

# Run the model as an NLS first
# Male Body Length NLS
bodylength.nls.male <- nls(BodyLength ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                           data = MaleSkel,
                           start = c(Asym = 19, AsymI = -1, 
                                     k = 0.001, kI = 0,
                                     t0 = 0),
                           na.action = na.omit)   # shouldn't need na.action
summary(bodylength.nls.male)   

## Female body length NLS
bodylength.nls.female <- nls(BodyLength ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                             data = FemaleSkel,
                             start = c(Asym = 16, AsymI = -1, 
                                       k = 0.001, kI = 0,
                                       t0 = 0),
                             na.action = na.omit)  # shouldn't need na.action
summary(bodylength.nls.female)   

# Now as a full NLMM 
# Male Body Length NLMM (model 9)
bodylength.nlme.male <- nlme(BodyLength ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                             data = MaleSkel,
                             fixed = (Asym+AsymI+k+kI+t0~1),
                             random = list(AnimalID = pdDiag(Asym+k+t0~ 1)),
                             correlation = corAR1(form = ~ 1|AnimalID),
                             start = c(Asym = 19, AsymI = -1, 
                                       k = 0.003, kI = 0,
                                       t0 = -200))
summary(bodylength.nlme.male) 
plot(bodylength.nlme.male, resid(.) ~ Age_days) 

# Male Body Length NLMM (model 10)
bodylength.nlme.female <- nlme(BodyLength ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                               data = FemaleSkel,
                               fixed = (Asym+AsymI+k+kI+t0~1),
                               random = list(AnimalID = pdDiag(Asym+k+t0~ 1)),
                               correlation = corAR1(form = ~ 1|AnimalID),
                               start = c(Asym = 17, AsymI = -0.5, 
                                         k = 0.003, kI = 0,
                                         t0 = -200))
summary(bodylength.nlme.female) 
plot(bodylength.nlme.female, resid(.) ~ Age_days) 

# Predict and plot the best fitting body length model in either case: 
# range(Male$Age) ; range(Female$Age)
newdf <- expand_grid(Age_days = 1:1200, Inbred = c(0, 1))

male_bodylength <- bind_cols(newdf, 
                             PredLength = predict(bodylength.nlme.male, 
                                                  newdata = newdf, 
                                                  level = 0)) %>% 
  mutate(Sex = "Male")

female_bodylength <- bind_cols(newdf, 
                               PredLength = predict(bodylength.nlme.female, 
                                                    newdata = newdf, 
                                                    level = 0)) %>% 
  mutate(Sex = "Female")

# Repeat steps to get confidence intervals
pars_male <- mvrnorm(nresamp, 
                     mu = fixef(bodylength.nlme.male), 
                     Sigma = vcov(bodylength.nlme.male))
pars_female <- mvrnorm(nresamp, 
                       mu = fixef(bodylength.nlme.female), 
                       Sigma = vcov(bodylength.nlme.female))

set.seed(101)
yvals1 <- matrix(nrow = length(male_bodylength$Age_days), ncol = 1000)
yvals2 <- matrix(nrow = length(female_bodylength$Age_days), ncol = 1000)
# 1000 cols for each Age. So I will predict for 1000 draws from the param distribution 
# (assuming multivariate normal dist)

# this time we need to predict over the monophasic form of the monomolecular curve
mono_inbred <- function(age,t0, Asym, AsymI, k, kI, t1 = 50, Inbred) {
  val1 <- (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - age)))
  return(val1)
}

# Male body length CI prediction
# There are a thousand parameter draws. 
# Each will be predicted against each age (to generate a 95% CI for each age from 0 to 1200 days)
for (i in 1:length(male_bodylength$Age_days)) {
  print(i)
  yvals1[i,] <- mono_inbred(age = male_bodylength$Age_days[i], 
                            Inbred = male_bodylength$Inbred[i],
                            Asym = pars_male[,1], 
                            AsymI = pars_male[,2], 
                            k = pars_male[,3], 
                            kI = pars_male[,4], 
                            t0 = pars_male[,5])
}

ci.males <- get_CI(yvals1)
male_bodylength <- cbind(male_bodylength, ci.males)  

# Female body length CI prediction
for (i in 1:length(female_bodylength$Age_days)) {
  print(i)
  yvals2[i,] <- mono_inbred(age = female_bodylength$Age_days[i], 
                            Inbred = female_bodylength$Inbred[i],
                            Asym = pars_female[,1], 
                            AsymI = pars_female[,2], 
                            k = pars_female[,3], 
                            kI = pars_female[,4], 
                            t0 = pars_female[,5])
}

ci.females <- get_CI(yvals2)
female_bodylength <- cbind(female_bodylength, ci.females)  

# Join the two sexes in a single data frame: 
twosex_bodylength <- bind_rows(male_bodylength, female_bodylength) %>% 
  mutate(Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age_days/365.35)

# Location of the 95% aysymptotic mass
male_bl_outbred_95 <- 0.95*19.21917                # male outbred
male_bl_inbred_95  <- 0.95*(19.21917 - 0.91626)  # male inbred
female_bl_outbred_95 <- 0.95*17.31919                # female outbred
female_bl_inbred_95  <- 0.95*(17.31919 -0.52058)   # female inbred 
male_bl_outbred_95 <- 0.95*19.21917                # male outbred

male_bl_outbred_95_age <- filter(twosex_bodylength, Sex == "Male" & Inbred == "Outbred") %>% 
  mutate(difflength = abs(PredLength - male_bl_outbred_95)) %>% 
  filter(difflength == min(difflength)) %>% 
  .$Age_days # 799

male_bl_inbred_95_age <- filter(twosex_bodylength, Sex == "Male" & Inbred == "Inbred") %>% 
  mutate(difflength = abs(PredLength - male_bl_inbred_95)) %>% 
  filter(difflength == min(difflength)) %>% 
  .$Age_days # 719

female_bl_outbred_95_age <- filter(twosex_bodylength, Sex == "Female" & Inbred == "Outbred") %>% 
  mutate(difflength = abs(PredLength - female_bl_outbred_95)) %>% 
  filter(difflength == min(difflength)) %>% 
  .$Age_days  # 690

female_bl_inbred_95_age <- filter(twosex_bodylength, Sex == "Female" & Inbred == "Inbred") %>% 
  mutate(difflength = abs(PredLength - female_bl_inbred_95)) %>% 
  filter(difflength == min(difflength)) %>% 
  .$Age_days # 658

# Make the body length plot
p_bodylength <- ggplot(twosex_bodylength, aes(x = Age_days, 
                                              y = PredWeight, 
                                              group = as.factor(Inbred), 
                                              colour = as.factor(Inbred), 
                                              fill = as.factor(Inbred))) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none") + 
  labs(y = "Body length (cm)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(9, 19, 2))

# Estimate the growth rate at each phase of life (g/day)
# Can do this easily by taking the successive deltas
twosex_bodylength <- twosex_bodylength %>% 
  group_by(Sex, Inbred) %>% 
  mutate(Delta = PredLength - lag(PredLength)) %>% 
  data.frame()

# To plot this we first plot by sex and then grob in the inset IGR
p_bodylength_males <- ggplot(filter(twosex_bodylength, Sex == "Male"),
                             aes(x = Age_days, y = PredLength, 
                                 group = as.factor(Inbred), 
                                 colour = as.factor(Inbred), 
                                 fill = as.factor(Inbred))) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none") + 
  labs(y = "Body length (cm)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(9, 19, 2), limits = c(8.5, 19.5))

p_maleblgrowthrate <- ggplot(data = filter(twosex_bodylength, Sex == "Male"), 
                             aes(x = Age_days,
                                 y = Delta, 
                                 group = as.factor(Inbred), 
                                 colour = as.factor(Inbred), 
                                 fill = as.factor(Inbred))) + 
  geom_path(linewidth = 0.6) + 
  facet_wrap(~ Sex) + 
  theme_classic() + 
  theme(axis.text.x = element_blank() , 
        axis.text.y = element_text(size = 8, colour = "black"), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 0.025))

p_bodylength_males <- p_bodylength_males + 
  annotation_custom(ggplotGrob(p_maleblgrowthrate), 
                    xmin = 510, xmax = 1200, 
                    ymin = 8.5, ymax = 14)

p_bodylength_females <- ggplot(filter(twosex_bodylength, Sex == "Female"),
                                 aes(x = Age_days, y = PredLength, 
                                     group = as.factor(Inbred), 
                                     colour = as.factor(Inbred), 
                                     fill = as.factor(Inbred))) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none") + 
  labs(y = "Body length (cm)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(9, 19, 2), limits = c(8.5, 19.5))

p_femaleblgrowthrate <- ggplot(data = filter(twosex_bodylength, Sex == "Female"), 
                               aes(x = Age_days,
                                   y = Delta, 
                                   group = as.factor(Inbred), 
                                   colour = as.factor(Inbred), 
                                   fill = as.factor(Inbred))) + 
  geom_path(linewidth = 0.6) + 
  facet_wrap(~ Sex) + 
  theme_classic() + 
  theme(axis.text.x = element_blank() , 
        axis.text.y = element_text(size = 8, colour = "black"), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 0.025))

p_bodylength_females <- p_bodylength_females + 
  annotation_custom(ggplotGrob(p_femaleblgrowthrate), 
                    xmin = 510, xmax = 1200, 
                    ymin = 8.5, ymax = 14)

p_bodylength2 <- (p_bodylength_females + 
                      theme(plot.margin = margin(r = 0), 
                            strip.background = element_blank(),
                            strip.text.x = element_blank(), 
                            axis.title.x = element_blank(), 
                            axis.text.x = element_blank())) | (p_bodylength_males + 
                                                                 theme(axis.text.y = element_blank(), 
                                                                       axis.ticks.y = element_blank(), 
                                                                       axis.title.y = element_blank(), 
                                                                       plot.margin = margin(l = 0), 
                                                                       strip.background = element_blank(),
                                                                       strip.text.x = element_blank(), 
                                                                       axis.title.x = element_blank(), 
                                                                       axis.text.x = element_blank()))

#---------------------------------------------------------------

## Modelling teeth width growth 

#----------------------------------------------------------------

# Run the model as an NLS first
# Male Teeth Width NLS
teethwidth.nls.male <- nls(UpTeethWidth ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                           data = MaleSkel,
                           start = c(Asym = 7, AsymI = -1, 
                                     k = 0.001, kI = 0,
                                     t0 = 0),
                           na.action = na.omit)   
summary(teethwidth.nls.male)   

## Females body length NLS
teethwidth.nls.female <-  nls(UpTeethWidth ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                              data = FemaleSkel,
                              start = c(Asym = 5.5, AsymI = 0, 
                                        k = 0.003, kI = 0,
                                        t0 = 0),
                              na.action = na.omit) # one NA
summary(teethwidth.nls.female)   

# Now as a full NLMM 
# Male Teeth Width NLMM (Model 11)
teethwidth.nlme.male <- nlme(UpTeethWidth ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                             data = MaleSkel,
                             fixed = (Asym+AsymI+k+kI+t0~1),
                             random = list(AnimalID = pdDiag(Asym+k+t0~ 1)),
                             correlation = corAR1(form = ~ 1|AnimalID),
                             start = c(Asym = 7, AsymI = -1, 
                                       k = 0.001, kI = 0,
                                       t0 = 0),
                             na.action = na.omit)   # shouldn't need na.action
summary(teethwidth.nlme.male)   

## Females body length NLMM (Model 12)
teethwidth.nlme.female <-  nlme(UpTeethWidth ~ (Asym + AsymI*Inbred)*(1 - exp((k + kI*Inbred)*(t0 - Age_days))), 
                                data = FemaleSkel,
                                fixed = (Asym+AsymI+k+kI+t0~1),
                                random = list(AnimalID = pdDiag(Asym+k+t0~ 1)),
                                correlation = corAR1(form = ~ 1|AnimalID),
                                start = c(Asym = 5.5, AsymI = 0, 
                                          k = 0.003, kI = 0,
                                          t0 = 0),
                                na.action = na.omit) # one NA
summary(teethwidth.nlme.female)   

# Predict and plot the best fitting body length model in either case: 
# range(Male$Age) ; range(Female$Age)
newdf <- expand_grid(Age_days = 1:1200, Inbred = c(0, 1))

male_teethwidth <- bind_cols(newdf, 
                             PredWidth = predict(teethwidth.nlme.male, 
                                                 newdata = newdf, 
                                                 level = 0)) %>% 
  mutate(Sex = "Male")

female_teethwidth <- bind_cols(newdf, 
                               PredWidth = predict(teethwidth.nlme.female, 
                                                   newdata = newdf, 
                                                   level = 0)) %>% 
  mutate(Sex = "Female")

# Get the confidence intervals
pars_male2 <- mvrnorm(nresamp, 
                      mu = fixef(teethwidth.nlme.male), 
                      Sigma = vcov(teethwidth.nlme.male))
pars_female2 <- mvrnorm(nresamp, 
                        mu = fixef(teethwidth.nlme.female), 
                        Sigma = vcov(teethwidth.nlme.female))

set.seed(101)
yvals1 <- matrix(nrow = length(male_teethwidth$Age_days), ncol = 1000)
yvals2 <- matrix(nrow = length(female_teethwidth$Age_days), ncol = 1000)

# Male teeth width length CI prediction
# There are a thousand parameter draws. 
# Each will be predicted against each age (to generate a 95% CI for each age from 0 to 1200 days)
for (i in 1:length(male_teethwidth$Age_days)) {
  print(i)
  yvals1[i,] <- mono_inbred(age = male_teethwidth$Age_days[i], 
                            Inbred = male_teethwidth$Inbred[i],
                            Asym = pars_male2[,1], 
                            AsymI = pars_male2[,2], 
                            k = pars_male2[,3], 
                            kI = pars_male2[,4], 
                            t0 = pars_male2[,5])
}

ci.males <- get_CI(yvals1)
male_teethwidth <- cbind(male_teethwidth, ci.males)  

# Female body length CI prediction
for (i in 1:length(female_teethwidth$Age_days)) {
  print(i)
  yvals2[i,] <- mono_inbred(age = female_teethwidth$Age_days[i], 
                            Inbred = female_teethwidth$Inbred[i],
                            Asym = pars_female2[,1], 
                            AsymI = pars_female2[,2], 
                            k = pars_female2[,3], 
                            kI = pars_female2[,4], 
                            t0 = pars_female2[,5])
}

ci.females <- get_CI(yvals2)
female_teethwidth <- cbind(female_teethwidth, ci.females)  

# Join the two sexes in a single data frame: 
twosex_teethwidth <- bind_rows(male_teethwidth, female_teethwidth) %>% 
  mutate(Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age_days/365.35)

# Location of the 95% aysymptotic mass
male_tw_outbred_95 <- 0.95*7.38076                  # male outbred
male_tw_inbred_95  <- 0.95*(7.38076 - 0.55119)  # male inbred
female_tw_outbred_95 <- 0.95*6.28108                # female outbred
female_tw_inbred_95  <- 0.95*(6.28108 -0.36013)   # female inbred 

male_tw_outbred_95_age <- filter(twosex_teethwidth, Sex == "Male" & Inbred == "Outbred") %>% 
  mutate(diffwidth = abs(PredWidth - male_tw_outbred_95)) %>% 
  filter(diffwidth == min(diffwidth)) %>% 
  .$Age_days # 961

male_tw_inbred_95_age <- filter(twosex_teethwidth, Sex == "Male" & Inbred == "Inbred") %>% 
  mutate(diffwidth = abs(PredWidth - male_tw_inbred_95)) %>% 
  filter(diffwidth == min(diffwidth)) %>% 
  .$Age_days # 832

female_tw_outbred_95_age <- filter(twosex_teethwidth, Sex == "Female" & Inbred == "Outbred") %>% 
  mutate(diffwidth = abs(PredWidth - female_tw_outbred_95)) %>% 
  filter(diffwidth == min(diffwidth)) %>% 
  .$Age_days  # 758

female_tw_inbred_95_age <- filter(twosex_teethwidth, Sex == "Female" & Inbred == "Inbred") %>% 
  mutate(diffwidth = abs(PredWidth - female_tw_inbred_95)) %>% 
  filter(diffwidth == min(diffwidth)) %>% 
  .$Age_days # 729

# Make the plot
p_teethwidth <- ggplot(twosex_teethwidth, aes(x = Age_days, 
                                                         y = PredWeight, 
                                                         group = as.factor(Inbred), 
                                                         colour = as.factor(Inbred), 
                                                         fill = as.factor(Inbred))) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none") + 
  labs(y = "Incisor width (mm)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(1, 8, 1))

# Estimate the growth rate at each phase of life (g/day)
# Can do this easily by taking the successive deltas
twosex_teethwidth <- twosex_teethwidth %>% 
  group_by(Sex, Inbred) %>% 
  mutate(Delta = PredWidth - lag(PredWidth)) %>% 
  data.frame()

# To plot this we first plot by sex and then grob in the inset IGR
p_teethwidth_females <- ggplot(filter(twosex_teethwidth, Sex == "Female"),
                                 aes(x = Age_days, y = PredWidth, 
                                     group = as.factor(Inbred), 
                                     colour = as.factor(Inbred), 
                                     fill = as.factor(Inbred))) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none") + 
  labs(y = "Incisor width (mm)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(1, 8, 1), 
                     labels = c("", 2, "", 4, "", 6, "", 8), 
                     limits = c(1, 8))

p_femaletwgrowthrate <- ggplot(data = filter(twosex_teethwidth, Sex == "Female"), 
                               aes(x = Age_days,
                                   y = Delta, 
                                   group = as.factor(Inbred), 
                                   colour = as.factor(Inbred), 
                                   fill = as.factor(Inbred))) + 
  geom_path(linewidth = 0.6) + 
  facet_wrap(~ Sex) + 
  theme_classic() + 
  theme(axis.text.x = element_blank() , 
        axis.text.y = element_text(size = 8, colour = "black"), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 0.017))

p_teethwidth_females <- p_teethwidth_females + 
  annotation_custom(ggplotGrob(p_femaletwgrowthrate), 
                    xmin = 510, xmax = 1200, 
                    ymin = 1, ymax = 4.7)

p_teethwidth_males <- ggplot(filter(twosex_teethwidth, Sex == "Male"),
                               aes(x = Age_days, y = PredWidth, 
                                   group = as.factor(Inbred), 
                                   colour = as.factor(Inbred), 
                                   fill = as.factor(Inbred))) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) + 
  geom_path(linewidth = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none") + 
  labs(y = "Incisor width (mm)", x = " Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) + 
  scale_y_continuous(breaks = seq(1, 8, 1), 
                     labels = c("", 2, "", 4, "", 6, "", 8), 
                     limits = c(1, 8))

p_maletwgrowthrate <- ggplot(data = filter(twosex_teethwidth, Sex == "Male"), 
                             aes(x = Age_days,
                                 y = Delta, 
                                 group = as.factor(Inbred), 
                                 colour = as.factor(Inbred), 
                                 fill = as.factor(Inbred))) + 
  geom_path(linewidth = 0.6) + 
  facet_wrap(~ Sex) + 
  theme_classic() + 
  theme(axis.text.x = element_blank() , 
        axis.text.y = element_text(size = 8, colour = "black"), 
        axis.title = element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  ylim(c(0, 0.017))

p_teethwidth_males <- p_teethwidth_males + 
  annotation_custom(ggplotGrob(p_maletwgrowthrate), 
                    xmin = 510, xmax = 1200, 
                    ymin = 1, ymax = 4.7)

p_teethwidth2 <- (p_teethwidth_females + 
                      theme(plot.margin = margin(r = 0), 
                            strip.background = element_blank(),
                            strip.text.x = element_blank(), 
                            axis.text.x = element_text(hjust = 0.5, size = 9.5))) | 
  (p_teethwidth_males + 
     theme(axis.text.y = element_blank(), 
           axis.ticks.y = element_blank(), 
           axis.title.y = element_blank(), 
           plot.margin = margin(l = 0), 
           strip.background = element_blank(),
           strip.text.x = element_blank(), 
           axis.text.x = element_text(hjust = 0.5, size = 9.5)))

#---------------------------------------------

# Lastly, we want to combine of the growth-related plots into a single figure

p_all <- (p_bodymass2 + 
            theme(plot.margin = margin(t = 5))) / 
  (p_bodylength2  + 
     theme(plot.margin = margin(t = 5))) / 
  (p_teethwidth2 + 
     theme(plot.margin = margin(t = 5)))

############################## END #########################################