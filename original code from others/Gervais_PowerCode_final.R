# retrieved from https://osf.io/nipv5/?view_only=533e0b6a3003477fbb7d7dd147b91ef5

library(pwr)
library(sciplot)
library(psych)
library(brglm)
library(ggplot2)
library(boot)

setwd("/Users/willgervais/Documents/Google Drive/Secrets/Research/Ongoing Projects/Power/Data & Code")


## Get Data

Full <- read.table("PowerData.csv", header=T, sep=",")
summary(Full)

nrow(Full)


Pow <- subset(Full, Choice < 1000)
summary(Pow)

nrow(Pow)





#############################################
####### Confirmatory Analyses: ##############
#### Tests of Registered Hypotheses #########
#############################################


##########################################
###### Hyp 1: Candidate Choices  #########
##### Will Differ Across Conditions ######
##########################################


## Does introducing power information change preferences?


## 2 x 3 chi square

X23 <- table(Pow$Choice, Pow$Condition)
chisq.test(X23)



#####################################
#### Pairwise Comparisons ###########
#####################################


#####################################
#### Hyp 2: Order of Preference #####
## For High Power Candidate Pattern #
## Findings < SS <<< Consequences ###
#####################################

# Findings Vs. Sample size

FS <- subset(Pow, Condition != "Consequences")
summary(FS)

FS.log <- glm(Choice ~ CondNum, data= FS, family= binomial)
summary(FS.log)
exp(coef(FS.log))
exp(confint(FS.log))

# Findings Vs. Consequences

FC <- subset(Pow, Condition != "Sample Size")
FC$CondNum <- FC$CondNum/2
summary(FC)

FC.log <- brglm(Choice ~ CondNum, data= FC, family= binomial)
summary(FC.log)
exp(coef(FC.log))
exp(confint(FC.log))


# SampleSize Vs. Consequences

SC <- subset(Pow, Condition != "Findings")
SC$CondNum <- SC$CondNum-1
summary(SC)

SC.log <- brglm(Choice ~ CondNum, data= SC, family= binomial)
summary(SC.log)
exp(coef(SC.log))
exp(confint(SC.log))

###########################################
############ Check preferences ############
########## within each condition ##########
###########################################


###########################################
##### Hyp 3: PPTS will significantly ######
#### prefer low-power candidate in the ####
## Findings and Sample Size Conditions. ###
###### This preference reverses in ########
###### the Consequences condition #########
########################################### 

Choices <- table(Pow$Choice, Pow$Condition)
prop.table(Choices, 2)
Choices


# dif from chance wintin each condition???


Findings <- subset(Pow, Condition == "Findings")
SS <- subset(Pow, Condition == "Sample Size")
Cons <- subset(Pow, Condition == "Consequences")


#Findings

F0 <- subset(Findings, Choice == "0")
F1 <- subset(Findings, Choice == "1")

F0b <- describe(F0$Choice)
F1b <- describe(F1$Choice)

binom.test(F1b$n, F1b$n + F0b$n, p = .5)

#Sample Size

SS0 <- subset(SS, Choice == "0")
SS1 <- subset(SS, Choice == "1")

SS0b <- describe(SS0$Choice)
SS1b <- describe(SS1$Choice)

binom.test(SS1b$n, SS1b$n + SS0b$n, p = .5)


#Consequences

Cons0 <- subset(Cons, Choice == "0")
Cons1 <- subset(Cons, Choice == "1")

Cons0b <- describe(Cons0$Choice)
Cons1b <- describe(Cons1$Choice)

binom.test(Cons1b$n, Cons1b$n + Cons0b$n, p = .5)












########################################
########################################
### Exploratory and Descriptive Stuff ##
########################################
########################################



########################
##### Describe #########
## research Practices###
########################


##Low

describe(Pow$NyearLow)

describe(Pow$DtypicalLow)

describe(Pow$NconditionLow)

describe(Pow$PcorrectLow)


##High

describe(Pow$NyearHigh)

describe(Pow$DtypicalHigh)

describe(Pow$NconditionHigh)

describe(Pow$PcorrectHigh)



##############################################
#### Calculate power, false pos rate, and ####
### replication rates from provided info #####
##############################################


##Using highest provided estimate 


powerHigh <- power.t.test(d=Pow$DtypicalHigh, power= NULL, n= Pow$NconditionHigh, sig.level=.05, alternative="two.sided")

Pow$N2.5 <- Pow$NconditionHigh * 2.5


powerHigh2.5 <- power.t.test(d=Pow$DtypicalHigh, power= NULL, n= Pow$N2.5, sig.level=.05, alternative="two.sided")


Pow$PowerTypicalHigh <- powerHigh$power


Pow$PowerTypicalHigh2.5 <- powerHigh2.5$power


rightHigh <- Pow$PcorrectHigh * Pow$PowerTypicalHigh

wrongHigh <- (1-Pow$PcorrectHigh)* .05

sigrateHigh <- rightHigh + wrongHigh


Pow$fprateHigh <- wrongHigh/sigrateHigh

Pow$reprateHigh <- ((1-Pow$fprateHigh) * Pow$PowerTypicalHigh) + (Pow$fprateHigh * .05)

Pow$reprateHigh2.5 <- ((1-Pow$fprateHigh) * Pow$PowerTypicalHigh2.5) + (Pow$fprateHigh * .05)





Pow$Power50 <- ifelse(Pow$PowerTypicalHigh <= .5, 1, 0)
Pow$Power80 <- ifelse(Pow$PowerTypicalHigh <= .8, 1, 0)
Pow$reprate50 <- ifelse(Pow$reprateHigh <= .5, 1, 0)

Pow$reprate50.2.5 <- ifelse(Pow$reprateHigh2.5 <= .5, 1, 0)

myvars <- c("PowerTypicalHigh" ,"fprateHigh", "reprateHigh", "reprateHigh2.5", "Power50", "Power80", "reprate50","reprate50.2.5")
PracticesHigh <- Pow[myvars]


describe(PracticesHigh)

NHigh <- describe(Pow$NconditionHigh)
PHigh <- describe(Pow$PowerTypicalHigh)
fHigh <- describe(Pow$fprateHigh)
rHigh <- describe(Pow$reprateHigh)


par(mfrow=c(2,2))

n <- density(Pow$NconditionHigh, na.rm=TRUE)
plot(n, main="N per conditon", xlim= c(0, 150))
polygon(n, col="cornflowerblue", border="black")
abline(v=NHigh$median, lty=3, col="green")


pwr <- density(Pow$PowerTypicalHigh, na.rm=TRUE)
plot(pwr, main="Power", xlim= c(0, 1))
polygon(pwr, col="cornflowerblue", border="black")
abline(v=PHigh$median, lty=3, col="green")



fp <- density(Pow$fprateHigh, na.rm=TRUE)
plot(fp, main="False Positive Rate", xlim= c(0, 1))
polygon(fp, col="cornflowerblue", border="black")
abline(v=fHigh$median, lty=3, col="green")


rep <- density(Pow$reprateHigh, na.rm=TRUE)
plot(rep, main="Replication Rate", xlim= c(0, 1))
polygon(rep, col="cornflowerblue", border="black")
abline(v=rHigh$median, lty=3, col="green")


######################################
######### Exploratory Analyses: ######
#### to see if institution, rank #####
## or practices (N & Power) predict ##
###### choice of candidates ##########
######################################

## Does the effect differ by academic rank?


Pow$N10 <- Pow$NconditionHigh/10
Pow$Pow10 <- Pow$PowerTypicalHigh/10

Findings <- subset(Pow, Condition == "Findings")
SS <- subset(Pow, Condition == "Sample Size")
Cons <- subset(Pow, Condition == "Consequences")



FRank <- table(Findings$Choice, Findings$RankID)
chisq.test(FRank)


SSRank <- table(SS$Choice, SS$RankID)
chisq.test(SSRank)

CRank <- table(Cons$Choice, Cons$RankID)
chisq.test(CRank)



## Does the effect differ by institution?


FInst <- table(Findings$Choice, Findings$InstID)
chisq.test(FInst)


SSInst <- table(SS$Choice, SS$InstID)
chisq.test(SSInst)

CInst <- table(Cons$Choice, Cons$InstID)
chisq.test(CInst)


## Does a participant's typical sample size relate to preferences when sample size is given? Consequences?

 




NS.log <- glm(Choice ~ N10, data=SS, family = binomial)
summary(NS.log)
exp(coef(NS.log))
exp(confint(NS.log))


NC.log <- glm(Choice ~ N10, data=Cons, family = binomial)
summary(NC.log)
exp(coef(NC.log))
exp(confint(NC.log))


## Does participant typical power relate to decisions?


#scale power to .10 increments



PS.log <- glm(Choice ~ Pow10, data=SS, family = binomial)
summary(PS.log)
exp(coef(PS.log))
exp(confint(PS.log))


PC.log <- brglm(Choice ~ Pow10, data=Cons, family = binomial)
summary(PC.log)
exp(coef(PC.log))
exp(confint(PC.log))









######################
## Demographics ######
######################

#Age
describe(Pow$Age)

#Gender
Pow$Gender <- factor(Pow$GenFem, levels = c("0", "1", "2"), labels= c("Male", "Female", "Other"))

GenTab <- table(Pow$Gender)
prop.table(GenTab)

#Institution

Pow$Institution <- factor(Pow$InstID, levels = c("1", "2", "3", "4", "5"), labels = c("R1", "R2", "MA", "SLAC", "Other"))

InTab <- table(Pow$Institution)
prop.table(InTab)

#Rank

Pow$Rank <- factor(Pow$RankID, levels = c("1", "2", "3", "4", "5"), labels = c("Assistant", "Associate", "Full", "Emeritus", "Other"))

RTab <- table(Pow$Rank)
prop.table(RTab)



#### Graph it ####
# Fig 2

Condition = c("Findings", "Sample Size", "Consequences")
prob = c(.217, .610, .966)
lcl = c(.121, .474, .880)
ucl = c(.342, .735, .996)
dat <- data.frame(Condition, prob, lcl, ucl)

summary(dat)

dat$Cond <- factor(dat$Condition, levels = c("Findings", "Sample Size", "Consequences"))

quartz()
ggplot(data=dat, aes(x= Cond, y = prob, ymin = lcl, ymax=ucl))  + geom_pointrange(col="black", size=.666) + labs(x= "\nCondition", y = "Proportion Choosing \nLarger Sample Size Candidate \n   ", title="Hiring Preferences  \n") + theme_bw() + scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) + geom_hline(yintercept=.5, lty=3, col="darkgrey") + geom_hline(yintercept=0, lty=1, col="darkgrey") + geom_hline(yintercept=1, lty=1, col="darkgrey") + theme(plot.title = element_text(size = rel(1.2), face="bold"))




