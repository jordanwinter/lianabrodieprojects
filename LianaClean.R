# Starting Clean 

#####libs####

library(tidyverse)
library(emmeans)
library(multcompView)
library(NatParksPalettes)

#load df

lianagerm <- read.csv("GerminationRate_Datasheetcsv.csv")
lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")


#################Germination Day Means##########################################

# MAKE DATA FRAME FOR GERMINATION DAY #

germdaydf <- lianabiotics %>% 
  select(Plant..,Treatment,Germinant.., Germination.Day)

# FIND AVERAGE GERM DAYS FOR EACH TREATMENT GROUP #

meangermday <- germdaydf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Germ_Day = mean(Germination.Day, na.rm = T))

# LET'S LOOK AT THE MEAN GERM DAY DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

germdaydfnona <- na.omit(germdaydf)

ggplot(germdaydfnona, aes(x = Treatment, y = Germination.Day, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) 


# STATISTICAL TESTS: START WITH AN ANOVA #

germdayaov <- aov(data = germdaydfnona, Germination.Day ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(germdayaov)

# p = 3.31e-07 # 

# NOW LET'S DO A TUKEYS #

TukeyHSD(germdayaov)

# add tukey sig indicators

# Tukey results
tuk <- TukeyHSD(germdayaov)

# Convert Tukey p-values to compact letter display
letters <- multcompLetters4(germdayaov, tuk)

# Extract letters for the Treatment factor
cld <- as.data.frame.list(letters$Treatment)
cld$Treatment <- rownames(cld)

# Get y-positions above each box
ypos <- germdaydfnona %>%
  group_by(Treatment) %>%
  summarise(y = max(Germination.Day, na.rm = TRUE) + 1, .groups = "drop")

plot_df <- left_join(cld, ypos, by = "Treatment")

ggplot(germdaydfnona, aes(x = Treatment, y = Germination.Day, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  geom_text(
    data = plot_df,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


######################Root Weight Normal#########################################

# NOW LET'S DO ROOT WEIGHT #

lianagerm <- read.csv("GerminationRate_Datasheetcsv.csv")
lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")


# MAKE DATA FRAME FOR ROOT WEIGHT #

rootwtdf <- lianabiotics %>% 
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

# add tukey sig indicators

# Tukey results
tuk <- TukeyHSD(rootwtaov)

# Convert Tukey p-values to compact letter display
letters <- multcompLetters4(rootwtaov, tuk)

# Extract letters for the Treatment factor
cld <- as.data.frame.list(letters$Treatment)
cld$Treatment <- rownames(cld)

# Get y-positions above each box
ypos <- rootwtdfnona %>%
  group_by(Treatment) %>%
  summarise(y = max(Root.Weight..g. + 0.0005, na.rm = TRUE), .groups = "drop")

plot_df <- left_join(cld, ypos, by = "Treatment")

ggplot(rootwtdfnona, aes(x = Treatment, y = Root.Weight..g., fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  coord_cartesian(ylim = c(0, 0.0075)) +
  geom_text(
    data = plot_df,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes =T,
    size = 5
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) 


######################Root Length Normal############################################

# NOW LET'S DO ROOT LENGTH # 

lianagerm <- read.csv("GerminationRate_Datasheetcsv.csv")
lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")

# MAKE DATA FRAME FOR ROOT LENGTH #

rootlengthdf <- lianabiotics %>% 
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

# Tukey results
tuk <- TukeyHSD(rootlengthaov)

# Convert Tukey p-values to compact letter display
letters <- multcompLetters4(rootlengthaov, tuk)

# Extract letters for the Treatment factor
cld <- as.data.frame.list(letters$Treatment)
cld$Treatment <- rownames(cld)

# Get y-positions above each box
ypos <- rootlengthdfnona %>%
  group_by(Treatment) %>%
  summarise(y = max(Root.Length..mm. + 10, na.rm = TRUE), .groups = "drop")

plot_df <- left_join(cld, ypos, by = "Treatment")

ggplot(rootlengthdfnona, aes(x = Treatment, y = Root.Length..mm., fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  geom_text(
    data = plot_df,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes =T,
    size = 5
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) 

#####################Total Weight Normal###########################################

# NOW LET'S DO TOTAL WEIGHT #

lianagerm <- read.csv("GerminationRate_Datasheetcsv.csv")
lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")

# MAKE DATA FRAME FOR TOTAL WEIGHT #

totalweightdf <- lianabiotics %>% 
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

# Tukey results
tuk <- TukeyHSD(totalweightaov)

# Convert Tukey p-values to compact letter display
letters <- multcompLetters4(totalweightaov, tuk)

# Extract letters for the Treatment factor
cld <- as.data.frame.list(letters$Treatment)
cld$Treatment <- rownames(cld)

# Get y-positions above each box
ypos <- totalweightdfnona %>%
  group_by(Treatment) %>%
  summarise(y = max(Total.Weight..g. + .005, na.rm = TRUE), .groups = "drop")

plot_df <- left_join(cld, ypos, by = "Treatment")

ggplot(totalweightdfnona, aes(x = Treatment, y = Total.Weight..g., fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  geom_text(
    data = plot_df,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes =T,
    size = 5
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) 

##########################Above Ground Biomass Normal#####################################################

# LASTLY, LET'S DO ABOVEGROUD WEIGHT #

lianagerm <- read.csv("GerminationRate_Datasheetcsv.csv")
lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")

# MAKE DATA FRAME FOR AGW #

abovegroundwtdf <- lianabiotics %>% 
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

# p = 1.75e-08

TukeyHSD(abovegroundwtaov)

# Tukey results
tuk <- TukeyHSD(abovegroundwtaov)

# Convert Tukey p-values to compact letter display
letters <- multcompLetters4(abovegroundwtaov, tuk)

# Extract letters for the Treatment factor
cld <- as.data.frame.list(letters$Treatment)
cld$Treatment <- rownames(cld)

# Get y-positions above each box
ypos <- abovegroundwtdfnona %>%
  group_by(Treatment) %>%
  summarise(y = max(AGB.Weight..g. + 0.0025, na.rm = TRUE), .groups = "drop")

plot_df <- left_join(cld, ypos, by = "Treatment")

ggplot(abovegroundwtdfnona, aes(x = Treatment, y = AGB.Weight..g., fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  geom_text(
    data = plot_df,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes =T,
    size = 5
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) 

######################Start of Ratio Work###############################################

# ok so lets make a df that has all biotics as a ratio of alive time?

lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")

##########Root Weight Ratio########

#Root wt

rootwtdf <- lianabiotics %>% 
  select(Plant..,Treatment,Germinant.., Root.Weight..g., Germination.Day) %>% 
  mutate(daysofgrowth = 39 - Germination.Day) %>% 
  mutate(rtwtratio = Root.Weight..g./daysofgrowth)

# remove row 559 maybe outlier?

rootwtdf <- rootwtdf[-c(559,407,445),]


# FIND AVERAGE ROOT WEIGHTS FOR EACH TREATMENT GROUP #

meanrootwt <- rootwtdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Root_Wt = mean(rtwtratio, na.rm = T))

# LET'S LOOK AT THE MEAN ROOT WT DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

rootwtdfnona <- na.omit(rootwtdf) %>% 
  filter(if_all(where(is.numeric), is.finite))

ggplot(rootwtdfnona, aes(x = Treatment, y = rtwtratio)) +
  geom_boxplot()

# STATISTICAL TESTS: START WITH AN ANOVA #

rootwtaov <- aov(data = rootwtdfnona, rtwtratio ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(rootwtaov)

# not significant if we do ratio for roowt

ggplot(rootwtdfnona, aes(x = Treatment, y = rtwtratio, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


#####################Root Length Ratio################################

#Root length

rootlengthdf <- lianabiotics %>% 
  select(Plant..,Treatment,Germinant.., Root.Length..mm., Germination.Day) %>% 
  mutate(daysofgrowth = 39 - Germination.Day) %>% 
  mutate(rtlnghthratio = Root.Length..mm./daysofgrowth)


# FIND AVERAGE ROOT LENGTHS FOR EACH TREATMENT GROUP #

meanrootlength <- rootlengthdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Root_Length = mean(rtlnghthratio, na.rm = T))

# LET'S LOOK AT THE MEAN ROOT LENGTH DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

rootlengthdfnona <- na.omit(rootlengthdf) %>% 
  filter(if_all(where(is.numeric), is.finite))

ggplot(rootlengthdfnona, aes(x = Treatment, y = rtlnghthratio)) +
  geom_boxplot()

# STATISTICAL TESTS: START WITH AN ANOVA #

rootlengthaov <- aov(data = rootlengthdfnona, rtlnghthratio ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(rootlengthaov)

# not significant as a ratio

ggplot(rootlengthdfnona, aes(x = Treatment, y = rtlnghthratio, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

##############Total Weight Ratio#############

# Total Wt

lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")

# MAKE DATA FRAME FOR TOTAL WEIGHT #

totalweightdf <- lianabiotics %>% 
  select(Plant..,Treatment,Germinant.., Total.Weight..g., Germination.Day) %>% 
  mutate(daysofgrowth = 39 - Germination.Day) %>% 
  mutate(totalweightratio = Total.Weight..g./daysofgrowth)


# FIND AVERAGE TOTAL WEIGHTS FOR EACH TREATMENT GROUP #

totalweightdf <- totalweightdf[-c(559,407,445),]


meantotalweight <- totalweightdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Total_Weight = mean(totalweightratio, na.rm = T))

# LET'S LOOK AT THE MEAN TOTAL WEIGHT DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

totalweightdfnona <- na.omit(totalweightdf) %>% 
  filter(if_all(where(is.numeric), is.finite))

ggplot(totalweightdfnona, aes(x = Treatment, y = totalweightratio)) +
  geom_boxplot()

# STATISTICAL TESTS: START WITH AN ANOVA #

totalweightaov <- aov(data = totalweightdfnona, totalweightratio ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(totalweightaov)

# is significant as a ratio! 0.0302

TukeyHSD(totalweightaov)

# Tukey results
tuk <- TukeyHSD(totalweightaov)

# Convert Tukey p-values to compact letter display
letters <- multcompLetters4(totalweightaov, tuk)

# Extract letters for the Treatment factor
cld <- as.data.frame.list(letters$Treatment)
cld$Treatment <- rownames(cld)

# Get y-positions above each box
ypos <- totalweightdfnona %>%
  group_by(Treatment) %>%
  summarise(y = max(totalweightratio + .0005, na.rm = TRUE), .groups = "drop")

plot_df <- left_join(cld, ypos, by = "Treatment")

ggplot(totalweightdfnona, aes(x = Treatment, y = totalweightratio, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  geom_text(
    data = plot_df,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes =T,
    size = 5
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) 

############Above Ground Biomass Ratio#############################

lianabiotics <- read.csv("TotalBioticsData_Part_2csv.csv")

# MAKE DATA FRAME FOR AGW #

abovegroundwtdf <- lianabiotics %>% 
  select(Plant..,Treatment,Germinant.., AGB.Weight..g., Germination.Day) %>% 
  mutate(daysofgrowth = 39 - Germination.Day) %>% 
  mutate(AGBratio = AGB.Weight..g./daysofgrowth)

# FIND AVERAGE ABOVEGROUND BIOMASS WEIGHTS FOR EACH TREATMENT GROUP #

abovegroundwtdf <- abovegroundwtdf[-c(559,407,445),]


meanabovegroundwt <- abovegroundwtdf %>% 
  group_by(Treatment) %>% 
  summarise(Mean_Aboveground_Weight = mean(AGBratio, na.rm = T))

# LET'S LOOK AT THE MEAN AGB WEIGHT DATA ON A PLOT #
# MAKE SURE THERE ARE NO NA VALUES #

abovegroundwtdfnona <- na.omit(abovegroundwtdf) %>% 
  filter(if_all(where(is.numeric), is.finite))

ggplot(abovegroundwtdfnona, aes(x = Treatment, y = AGBratio)) +
  geom_boxplot() +
  theme_light()

# STATISTICAL TESTS: START WITH AN ANOVA #

abovegroundwtaov <- aov(data = abovegroundwtdfnona, AGBratio ~ Treatment)

# LOOK AT IT IN A MORE DIGESTIBLE WAY #

summary(abovegroundwtaov)

# It's 0.0513 which is so annoying - trending towards significance as a ratio

ggplot(abovegroundwtdfnona, aes(x = Treatment, y = AGBratio, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_natparks_d("Yellowstone") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )




#############Germination Probability Stats and Graphs############ 

library(tidyverse)
library(emmeans)
library(ggplot2)
library(multcomp)
library(multcompView)
library(NatParksPalettes)

dat <- read.csv("GerminationRate_Datasheetcsv.csv")
dat <- dat[-c(2401,2402),]

dat$Treatment <- as.factor(dat$Treatment)
dat$Plant.. <- as.factor(dat$Plant..)

counts <- dat %>%
  group_by(Plant..) %>%
  summarise(nonzero_germ = sum(Germination.Day != 0),
            .groups = "drop") %>%
  mutate(total = 20)

joined <- left_join(dat, counts, by = "Plant..")

df_ <- joined %>%
  group_by(Plant.., Treatment) %>%
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data = df_,
  family = binomial
)

emm <- emmeans(m_glm, ~ Treatment)
tukey_letters <- cld(emm, adjust = "tukey", Letters = letters)

newdat <- data.frame(Treatment = levels(df_$Treatment))

preds <- predict(m_glm, newdata = newdat, type = "link", se.fit = TRUE)

ilogit <- function(x) exp(x) / (1 + exp(x))

newdat$pred_link <- preds$fit
newdat$se_link <- preds$se.fit
newdat$lwr_link <- newdat$pred_link - 1.96 * newdat$se_link
newdat$upr_link <- newdat$pred_link + 1.96 * newdat$se_link
newdat$pred_prob <- ilogit(newdat$pred_link)
newdat$lwr_prob <- ilogit(newdat$lwr_link)
newdat$upr_prob <- ilogit(newdat$upr_link)

newdat <- dplyr::left_join(
  newdat,
  dplyr::select(as.data.frame(tukey_letters), Treatment, .group),
  by = "Treatment"
)
ggplot(newdat, aes(x = Treatment, y = pred_prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = .22), size = 5) +
  labs(x = "Treatment",
       y = "Predicted germination probability",
       title = "Fixed-effect treatment estimates") +
  scale_fill_natparks_d("Yellowstone") +
  theme_bw()

n_treat <- length(unique(newdat$Treatment))

yellowstone_cols <- NatParksPalettes::natparks.pals("Yellowstone", n_treat)

yellowstone_cols
length(yellowstone_cols)

names(yellowstone_cols) <- levels(newdat$Treatment)

ggplot(newdat, aes(x = Treatment, y = pred_prob, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = .22), size = 5, color = "black") +
  scale_color_manual(values = yellowstone_cols, drop = FALSE) +
  labs(
    x = "Treatment",
    y = "Predicted germination probability",
    title = "Fixed-effect treatment estimates"
  ) +
  theme_bw()

# I don't know why I need to make both graphs for the theme to work but I don't care rn.