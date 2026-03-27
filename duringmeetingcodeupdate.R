# Part 2 R Data #

# SET WD BY CLICKING "MORE" THEN "SET AS WORKIND DIRECTORY" #

setwd("C:/Users/liana/OneDrive/Desktop/Part_2_R_Files/RdataforP2")

# LOAD TIDYVERSE #

library(tidyverse)

# LOAD DATA #

originaldata2 <- read.csv("TotalBioticsData_Part_2csv.csv")

# MAKE DATA FRAME FOR GERMINATION DAY #

germdaydf <- originaldata2 %>% 
  select(Plant..,Treatment,Germinant.., Germination.Day)

# FIND AVERAGE GERM DAYS FOR EACH TREATMENT GROUP #

meangermday <- germdaydf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Germ_Day = mean(Germination.Day, na.rm = T))

# LET'S LOOK AT THE MEAN GERM DAY DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

germdaydfnona <- na.omit(germdaydf)

ggplot(germdaydfnona, aes(x = Treatment, y = Germination.Day, col = Treatment)) +
  geom_boxplot() 
  

# STATISTICAL TESTS: START WITH AN ANOVA #

germdayaov <- aov(data = germdaydfnona, Germination.Day ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(germdayaov)

# p = 3.31e-07 # 

# NOW WE NEED TO DO A POST-HOC TEST. LET'S DO BONFERRONI #
# BONFERRONI DOES A BUNCH OF T-TESTS BETWEEN VARIABLES TO SEE
# WHICH ONE IS SIGNIFICANT #

pairwise.t.test(germdaydfnona$Germination.Day,
                germdaydfnona$Treatment,
                p.adjust.method = "bonferroni")

# BEEFERONI :P LOOK AT RESULTS #

# ...

# NOW LET'S DO A TUKEYS #

TukeyHSD(germdayaov)

# ... 

###############################################################

# NOW LET'S DO ROOT WEIGHT #

# MAKE DATA FRAME FOR ROOT WEIGHT #

rootwtdf <- originaldata2 %>% 
  select(Plant..,Treatment,Germinant.., Root.Weight..g.)

# remove row 559 maybe outlier?

rootwtdf <- rootwtdf[-c(559,407,445),]


# FIND AVERAGE ROOT WEIGHTS FOR EACH TREATMENT GROUP #

meanrootwt <- rootwtdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Root_Wt = mean(Root.Weight..g., na.rm = T))

# LET'S LOOK AT THE MEAN ROOT WT DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

rootwtdfnona <- na.omit(rootwtdf)

ggplot(rootwtdfnona, aes(x = Treatment, y = Root.Weight..g.)) +
  geom_boxplot()

# STATISTICAL TESTS: START WITH AN ANOVA #

rootwtaov <- aov(data = rootwtdfnona, Root.Weight..g. ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(rootwtaov)

# rn it's sig with the 3 outliers removed. see line 73

# do a posty

TukeyHSD(rootwtaov)

##################################################################

# NOW LET'S DO ROOT LENGTH # 

# MAKE DATA FRAME FOR ROOT LENGTH #

rootlengthdf <- originaldata2 %>% 
  select(Plant..,Treatment,Germinant.., Root.Length..mm.)

# FIND AVERAGE ROOT LENGTHS FOR EACH TREATMENT GROUP #

meanrootlength <- rootlengthdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Root_Length = mean(Root.Length..mm., na.rm = T))

# LET'S LOOK AT THE MEAN ROOT LENGTH DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

rootlengthdfnona <- na.omit(rootlengthdf)

ggplot(rootlengthdfnona, aes(x = Treatment, y = Root.Length..mm.)) +
  geom_boxplot()

# STATISTICAL TESTS: START WITH AN ANOVA #

rootlengthaov <- aov(data = rootlengthdfnona, Root.Length..mm. ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(rootlengthaov)

# p= 0.0121 #

#posty 

TukeyHSD(rootlengthaov)
############################################################################

# NOW LET'S DO TOTAL WEIGHT #

# MAKE DATA FRAME FOR TOTAL WEIGHT #

totalweightdf <- originaldata2 %>% 
  select(Plant..,Treatment,Germinant.., Total.Weight..g.)

# FIND AVERAGE TOTAL WEIGHTS FOR EACH TREATMENT GROUP #

totalweightdf <- totalweightdf[-c(559,407,445),]


meantotalweight <- totalweightdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Total_Weight = mean(Total.Weight..g., na.rm = T))

# LET'S LOOK AT THE MEAN TOTAL WEIGHT DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

totalweightdfnona <- na.omit(totalweightdf)

ggplot(totalweightdfnona, aes(x = Treatment, y = Total.Weight..g.)) +
  geom_boxplot()

# STATISTICAL TESTS: START WITH AN ANOVA #

totalweightaov <- aov(data = totalweightdfnona, Total.Weight..g. ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(totalweightaov)

# p= 0.00273 #

TukeyHSD(totalweightaov)

###############################################################################

# LASTLY, LET'S DO ABOVEGROUD WEIGHT #

# MAKE DATA FRAME FOR AGW #

abovegroundwtdf <- originaldata2 %>% 
  select(Plant..,Treatment,Germinant.., AGB.Weight..g.)

# FIND AVERAGE ABOVEGROUND BIOMASS WEIGHTS FOR EACH TREATMENT GROUP #

abovegroundwtdf <- abovegroundwtdf[-c(559,407,445),]


meanabovegroundwt <- abovegroundwtdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Aboveground_Weight = mean(AGB.Weight..g., na.rm = T))

# LET'S LOOK AT THE MEAN AGB WEIGHT DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

abovegroundwtdfnona <- na.omit(abovegroundwtdf)

ggplot(abovegroundwtdfnona, aes(x = Treatment, y = AGB.Weight..g.)) +
  geom_boxplot() +
  theme_light()

# STATISTICAL TESTS: START WITH AN ANOVA #

abovegroundwtaov <- aov(data = abovegroundwtdfnona, AGB.Weight..g. ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(abovegroundwtaov)

# p = 

TukeyHSD(abovegroundwtaov)

##################################################################

# END CODE #