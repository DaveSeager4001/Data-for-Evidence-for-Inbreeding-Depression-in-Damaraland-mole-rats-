

lapply(c("coxme", "glmmTMB", "survminer", "patchwork", "tidyverse","DHARMa","nlme","ggplot2","MASS"), FUN = library, character.only = TRUE)


# Code supporting paper "Inbreeding depression in Damaraland Mole-rats"

# Model 1: Modelling Litter Size ~ treatment (n = 109)

# To model litter size we used a dataset of 109 litters for which litter size is known.
## The response variable is litter size
### The fixed effects tested were treatment (inbred/outbred), colony Size (number of individuals in group) and density (number of animals per metre of tunnel)
#### maternal ID (mother of the litter) is included as a random effect

model1<-glmmTMB(LitterSize~Treatment + ColonySize + Density + (1|MaternalID),data=LitterSize_Dataset)
model1<-glmmTMB(LitterSize~Treatment + Density + (1|MaternalID),data=LitterSize_Dataset)
summary(model1)
residuals1 = simulateResiduals(model1)
plot(residuals1)
cor(LitterSize_Dataset %>% select(GroupSize,Density,Mat.Age.years))
cor.test(LitterSize_Dataset$Density,LitterSize_Dataset$GroupSize)
# None of the fixed effects predict litter size 

#---------------------------------------------------------------------#

# Model 2: Modelling Inter-birth interval ~ treatment (n = 117)

# To model inter-birth Interval we used a dataset of 118 litters
## 1 litter removed due to unusually long time between litters (818 days)
## The response variable is inter-birth interval in days - gestation period (90 days)
### The fixed effects tested were treatment, colony size (number of individuals in group) and density (number of animals per metre of tunnel)
#### maternal ID is included as a random effect

model2<-glmmTMB(InterBirth_rel~Treatment + GroupSize + Density + (1|MaternalID),
                data=InterBirthInterval_Dataset1, 
                family = Gamma(link = "log"))
summary(model2)
residuals2 = simulateResiduals(model2)
plot(residuals2)
# None of the fixed effects predict inter-birth interval but there is a trend for larger colonies to have a shorter interbirth interval 

#---------------------------------------------------------------------#

# Model 3: Modelling birth weight ~ treatment (n = 292)

# To model offspring mass at birth and its effect on survival we used a dataset of 292 individuals with a weight taken within 5 days of birth
## The response variable was offspring weight at birth in grams
### The fixed effects were Treatment, colony size, Litter size, density and the difference in days between birth and the weight being taken (maximum 5)
#### maternal ID and Litter Reference are included as random effects

model3<-glmmTMB(Weight~Treatment + GroupSize + LitterSize + Density + Wt_Diff + (1|Mother) + (1|LitterRef),data=Pup_Weight_Dataset)
summary(model3)
residuals3 = simulateResiduals(model3)
plot(residuals3)
cor.test(Pup_Weight_Dataset$GroupSize,Pup_Weight_Dataset$Mat.Litter)

# Individuals from outbred colonies are heavier than those from inbred colonies
## Individuals in larger colonies are lighter
### Individuals in larger litters are lighter
#### Weight increases with the number of days that animals are weighed after birth


# plot birth weight ~ treatment to produce plot 1a
newdf <- data.frame(Treatment = unique(Pup_Weight_Dataset$Treatment), 
                    ColonySize = mean(Pup_Weight_Dataset$ColonySize), 
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

# To model the effect of offspring mass at birth on survival we used a dataset of 292 individuals with a weight taken within 5 days of birth
## The response variable was individual survival to 30 days
### The fixed effects were Treatment, colony size, weight at birth  and density
#### maternal ID and Litter Reference are included as random effects

# To test how birthweight affects survival to 30 days
model4<-glmmTMB(Sur30 ~ Treatment + ColonySize + Weight + Density + 
                  (1|Mother) + (1|LitterRef),
                data=Pup_Weight_Dataset, 
                family = binomial) 
summary(model4)
residuals4 = simulateResiduals(model4)
plot(residuals4)

# Offspring weight predicts survival to 30 days with lighter individuals more likely to die early

# plot survival to 30 days ~ birth weight to produce plot 1b
newdf2 <- data.frame(Treatment = unique(Pup_Weight_Dataset$Treatment), 
                     ColonySize = mean(Pup_Weight_Dataset$ColonySize), 
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
newdf2$l95 <-  boot::inv.logit(newdf2$pred - 1.96*newdf2$SE)
newdf2$u95 <-  boot::inv.logit(newdf2$pred + 1.96*newdf2$SE)
newdf2$pred <-  boot::inv.logit(newdf2$pred)

p2 <- ggplot(Pup_Weight_Dataset,aes(x=Weight, y = Sur30, group = Treatment, colour = Treatment)) +
  geom_jitter(height=0.05, width = 0.1, alpha= 0.4, shape = 1) +
  geom_line(data = newdf2, aes(y = l95), linewidth = 0.5, linetype = 2) +
  geom_line(data = newdf2, aes(y = u95), linewidth = 0.5, linetype = 2) +
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
  scale_x_continuous(breaks = seq(6, 14, 1),
                     labels = c("6", "", "8", "", "10", "", "12", "", "14"))

#---------------------------------------------------------------

# Model 5: Model the interaction between status and colony size to survival at 30 days (n = 328)

# To model the interaction between colony size and survival to 30 days we used a dataset of 328 individuals 
## The response variable was survival to 30 days
### The fixed effects were Treatment, colony size and density
#### maternal ID and Litter Reference are included as random effects

model5 <-glmmTMB(Sur30 ~ Density + Treatment * ColonySize + 
                   (1|Mother) + (1|LitterRef),
                 data=Survival_Dataset,family = "binomial") 
summary(model5)
residuals5 = simulateResiduals(model5)
plot(residuals5)

# There is a significant interaction between colony size and survival to 30 days with inbred individuals in smaller colonies being less likely to survive

# plot the interaction effect of survival to 30 days ~ colony size * treatment to produce plot 1c
newdf3 <- expand.grid(Treatment = unique(Survival_Dataset$Treatment), 
                      ColonySize = seq(2, 26, 1), 
                      Density = mean(Survival_Dataset$Density), 
                      Mother = unique(Survival_Dataset$Mother)[1], 
                      LitterRef = unique(Survival_Dataset$LitterRef)[1])

# restrict the plot over the weight range of each treatment 
range(Survival_Dataset$ColonySize[Survival_Dataset$Treatment == "Inbred" & !is.na(Survival_Dataset$Sur30)])
range(Survival_Dataset$ColonySize[Survival_Dataset$Treatment == "Outbred" & !is.na(Survival_Dataset$Sur30)])
newdf3 <- newdf3 %>% 
  filter(!(Treatment == "Inbred" & !between(ColonySize, 2, 26))) %>% 
  filter(!(Treatment == "Outbred" & !between(ColonySize, 2, 24))) 

newdf3$pred <- predict(model5, newdata = newdf3, re.form = NA, type = "link")
newdf3$SE <- predict(model5, newdata = newdf3, re.form = NA, se.fit = TRUE, type = "link")$se.fit
newdf3$l95 <-  boot::inv.logit(newdf3$pred - 1.96*newdf3$SE)
newdf3$u95 <-  boot::inv.logit(newdf3$pred + 1.96*newdf3$SE)
newdf3$pred <-  boot::inv.logit(newdf3$pred)

p3 <- ggplot(Survival_Dataset,aes(x=ColonySize, y = Sur30, group = Treatment, 
                          colour = Treatment)) +
  geom_jitter(height=0.05, width = 0.1, alpha=0.4, shape = 1) +
  geom_line(data = newdf3, aes(y = l95), linewidth = 0.5, linetype = 2) +
  geom_line(data = newdf3, aes(y = u95), linewidth = 0.5, linetype = 2) +
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
  scale_x_continuous(breaks = seq(2, 26, 4))

#---------------------------------------------------------------

# Model 6: Use Cox models to model the difference in survival between treatments

# To model the difference in individual survival between treatments  we used a dataset of 328 individuals 
## The fixed effects were Treatment and density
#### maternal ID and Litter Reference are included as random effects

model6 <- coxme(Surv(Mnths,Censorship_D) ~ Treatment + Density + ColonySize + 
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
  scale_x_continuous(breaks = seq(0, 48, 12), limits = c(0, 50)) #+ 
#annotate("text", x = 36, y = 0.15, label = "Inbred", colour = "indianred") + 
#annotate("text", x = 36, y = 0.25, label = "Outbred", colour = "dodgerblue")

# Combine the plots together
(p1 | p2) / (p3 | p4)

######################### END #########################################

# Model 7: Compare the weight curves of inbred and outbred individuals

# The datset InbreedingDepresssion_Weights_Dataset has 9126 weights from 92 individuals (40 females and 52 males)  

# To create a duplicate dataset
IBWeights<-InbreedingDepression_Weights_Dataset

# Due to the high degree of sexual dimorphism the sexes should be modelled sperately
## To create seperate datasets for male and female weights
Female <- filter(IBWeights, Sex == "F") %>% 
  arrange(AnimalID, Age)
Male   <- filter(IBWeights, Sex == "M") %>% 
  arrange(AnimalID, Age)

# To check the number of females and observations and provide a median and SD
length(unique(Female$AnimalID)) # n = 40
nrow(Female)
Female %>% 
  group_by(AnimalID) %>% 
  summarise(n = n())  %>% 
  summarise(median(n), sd(n))

# To check the number of males and observations and provide a median and SD
length(unique(Male$AnimalID)) # n = 52
nrow(Male)
Male %>% 
  group_by(AnimalID) %>% 
  summarise(n = n())  %>% 
  summarise(median(n), sd(n))


# I will include inbreeding as a dummy var coded as 0/1
Male$Inbred <- if_else(Male$Status == "IB", 1, 0)
Female$Inbred <- if_else(Female$Status == "IB", 1, 0)

#------------------------

# Plot the raw data to see the spread of weights ~ age
## For females
ggplot(Female, aes(x = Age, y = Weight, 
                   group = as.factor(Status), 
                   col = as.factor(Status))) + 
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point() +
  geom_smooth(method = "gam", colour = "black", 
              formula = y ~ s(x, k = 6, bs = "tp")) + 
  facet_wrap(~AnimalID) +
  xlim(c(0, 1250)) +
  ylim(c(0, 250)) +
  theme_bw() + 
  theme(panel.grid = element_blank())

# And for males
ggplot(Male, aes(x = Age, y = Weight, 
                 group = as.factor(Status), 
                 col = as.factor(Status))) + 
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point() +
  geom_smooth(method = "gam", colour = "black", 
              formula = y ~ s(x, k = 6, bs = "tp")) + 
  facet_wrap(~AnimalID) +
  xlim(c(0, 1250)) +
  ylim(c(0, 250)) +
  theme_bw() + 
  theme(panel.grid = element_blank())

## Modelling growth with a monomolecular model 

## First fit the model without random effects just to get an idea of good starting parameters
# For that still need just one value per day for each  

# The model formula: 
biphasic_inbred <- function(age,t0, Asym, AsymI, k1, k1I, k2, k2I, t1 = 50, Inbred) {
  
  val1 <- (Asym + AsymI*Inbred)*(1 - exp((k1 + k1I*Inbred)*(t0 - age)))
  val2 <- (Asym + AsymI*Inbred)*(1 - exp((k1 + k1I*Inbred)*(t0 - t1) + (k2 + k2I*Inbred)*(t1 - age)))
  val1[age > t1] <- val2[age > t1]
  
  return(val1)
}

# As before, run the model as an nls first
## Males 
biphasic.nls.male <- nls(Weight ~ biphasic_inbred(age=Age,t0=t0,Asym=Asym, AsymI=AsymI, 
                                                  k1=k1,k1I=k1I,k2=k2, k2I=k2I, Inbred = Inbred),
                         data = Male,
                         start = c(Asym = 150, AsymI = -30, k1 = 0.001, k1I = 0, k2 = 0.001, k2I = 0, 
                                   t0 = 0),
                         na.action = na.pass) 
summary(biphasic.nls.male)   

## Females 
biphasic.nls.female <-nls(Weight ~ biphasic_inbred(age=Age,t0=t0,Asym=Asym, AsymI=AsymI, 
                                                   k1=k1,k1I=k1I,k2=k2, k2I=k2I, Inbred = Inbred),
                          data = Female,
                          start = c(Asym = 150, AsymI = -30, k1 = 0.001, k1I = 0.0001, k2 = 0.001, k2I = 0.0001, 
                                    t0 = 0),
                          na.action = na.pass) 
summary(biphasic.nls.female) 


# Fit as nlme
## Males
biphasic.nlme.male <- nlme(Weight ~ biphasic_inbred(age=Age,Asym=Asym,AsymI=AsymI, 
                                                    k1=k1,k1I=k1I,k2=k2, k2I=k2I, t0=t0, Inbred = Inbred),
                           data = Male,
                           fixed = (Asym+AsymI+k1+k1I+k2+k2I+t0~1),
                           random = list(AnimalID = pdDiag(Asym+k1+k2 +t0~ 1)),
                           correlation = corAR1(form = ~ 1|AnimalID),
                           start = c(Asym = 250, AsymI = -40, k1 = 0.002, k1I = 0, 
                                     k2 = 0.001, k2I = 0, t0 = 0),
                           na.action = na.pass)
summary(biphasic.nlme.male) 
plot(biphasic.nlme.male, resid(.) ~ Age) 

## Females
biphasic.nlme.female <- nlme(Weight ~ biphasic_inbred(age=Age,t0=t0,Asym=Asym, AsymI=AsymI, 
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
summary(biphasic.nlme.female) 
plot(biphasic.nlme.female, resid(.) ~ Age) 


# Predict and plot the best fitting model in either case: 
# range(Male$Age) ; range(Female$Age)
newdf <- expand_grid(Age = 1:1500, Inbred = c(0, 1))

male_pop <- bind_cols(newdf, PredWeight = predict(biphasic.nlme.male, newdata = newdf, level = 0)) %>% 
  mutate(Sex = "Male")
female_pop <- bind_cols(newdf, PredWeight = predict(biphasic.nlme.female, newdata = newdf, level = 0)) %>% 
  mutate(Sex = "Female")

# Get the confidence intervals (population prediction intervals in the parlance of: 
# "https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme)"

# now get the bootstrapped values for the mean
nresamp <- 1000

## pick new parameter values by sampling from multivariate normal distribution based on fit
pars_male <- mvrnorm(nresamp, 
                     mu = fixef(biphasic.nlme.male), 
                     Sigma = vcov(biphasic.nlme.male))
pars_female <- mvrnorm(nresamp, 
                       mu = fixef(biphasic.nlme.female), 
                       Sigma = vcov(biphasic.nlme.female))

## utility function
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
# There are a thousand parameter draws. 
# Each will be predicted against each age (to generate a 95% CI for each age from 0 to 1200 days)

for (i in 1:length(male_pop$Age))
{
  print(i)
  yvals2[i,] <- biphasic_inbred(age = male_pop$Age[i], 
                                Inbred = male_pop$Inbred[i],
                                Asym = pars_male[,1], 
                                AsymI = pars_male[,2], 
                                k1 = pars_male[,3], 
                                k1I = pars_male[,4], 
                                k2 = pars_male[,5], 
                                k2I = pars_male[,6], 
                                t0 = pars_male[,7])
}

ci.males <- get_CI(yvals2)
male_pop <- cbind(male_pop, ci.males)  

# Female prediction
for (i in 1:length(female_pop$Age))
{
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

# add a rug to show where the data lies

# Produce the plot to this point 
p_inbreeding <- ggplot(twosex_pop, aes(x = Age, y = PredWeight, 
                                       group = as.factor(Inbred), 
                                       colour = as.factor(Inbred), 
                                       fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, colour = NA) + 
  geom_path(lwd = 1) + 
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


#--------------------------------


# Get the age at which 95% asymptotic mass reached
# Suggesting they haven't yet reached their asymptotic body mass
twosex_pop %>% 
  mutate(PopAsymp = case_when(Sex == "Male" & Inbred == "Outbred" ~ fixef(biphasic.nlme.male)[1], 
                              Sex == "Male" & Inbred == "Inbred" ~ fixef(biphasic.nlme.male)[1] + fixef(biphasic.nlme.male)[2], 
                              Sex == "Female" & Inbred == "Outbred" ~ fixef(biphasic.nlme.female)[1], 
                              Sex == "Female" & Inbred == "Inbred" ~ fixef(biphasic.nlme.female)[1] + fixef(biphasic.nlme.female)[2])) %>% 
  mutate(PopAsymp_95 = 0.95*PopAsymp) %>% 
  group_by(Sex, Inbred) %>% 
  mutate(MassDiff = abs(PopAsymp_95 - PredWeight)) %>%
  filter(MassDiff == min(MassDiff)) %>% 
  data.frame()
# Estimate the growth rate at each phase of life (g/day)
# Can do this easily by taking the successive deltas

twosex_pop <- twosex_pop %>% 
  group_by(Sex, Inbred) %>% 
  mutate(Delta = PredWeight - lag(PredWeight)) %>% 
  data.frame()

p_growthrates <- ggplot(data = filter(twosex_pop, Age <= 50), aes(x = Age, y = Delta, 
                                                                  group = as.factor(Inbred), 
                                                                  colour = as.factor(Inbred), 
                                                                  fill = as.factor(Inbred))) + 
  geom_vline(xintercept = 50, linetype = 2, colour = "gray") + 
  geom_path(lwd = 1) + 
  geom_path(data = filter(twosex_pop, Age > 50), lwd = 1) + 
  facet_wrap(~ Sex) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10.5, colour = "black") , 
        axis.text.x = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 11, colour = "black"), 
        legend.position = "none", 
        legend.text = element_text(size = 10), 
        legend.title = element_blank()) + 
  labs(y = "Daily growth rate (g/day)", x = "Age (days)") + 
  scale_colour_manual(values = c("indianred", "dodgerblue")) + 
  scale_fill_manual(values = c("indianred", "dodgerblue")) 


p_inbreeding2 <- p_inbreeding + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position = c(0.12, 0.88))

p_inbreeding2 / p_growthrates

#--------------------------------

### Generate predictions for each individual to see how the model captures individual variation

newdf_female <- Female %>% 
  group_by(AnimalID, Inbred) %>%
  expand(Age = seq(1, max(Age), 1)) %>% 
  data.frame()

newdf_male <- Male %>% 
  group_by(AnimalID, Inbred) %>%
  expand(Age = seq(1, max(Age), 1)) %>% 
  data.frame()

female_ind <- bind_cols(newdf_female, PredWeight = predict(biphasic.nlme.female, 
                                                           newdata = newdf_female)) %>% 
  mutate(Sex = "Female", 
         Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age/365.35)


male_ind <- bind_cols(newdf_male, PredWeight = predict(biphasic.nlme.male, 
                                                       newdata = newdf_male)) %>% 
  mutate(Sex = "Male", 
         Inbred = if_else(Inbred == 0, "Outbred", "Inbred"), 
         Age_yrs = Age/365.35)


# Female individual curves
Female_toplot <- mutate(Female, Inbred = if_else(Inbred == 0, "Outbred", "Inbred"))

p_females <- ggplot(Female_toplot, aes(x = Age, y = Weight, group = Inbred, col = Inbred)) +
  geom_hline(yintercept = 10, linetype = 2) + 
  geom_point(alpha = 0.1) +
  geom_line(data = female_ind, aes(x = Age, y = PredWeight), 
            colour = "black", lwd = 0.8) + 
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
  scale_x_continuous(breaks = seq(0, 1500, 250))


# Male individual curves
Male_toplot <- mutate(Male, Inbred = if_else(Inbred == 0, "Outbred", "Inbred"))

p_males <- ggplot(Male_toplot, aes(x = Age, y = Weight, group = Inbred, col = Inbred)) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point(alpha = 0.1) +
  geom_line(data = male_ind, aes(x = Age, y = PredWeight), 
            colour = "black", lwd = 0.8) + 
  facet_wrap(~AnimalID) +
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
  scale_x_continuous(breaks = seq(0, 1500, 250))

