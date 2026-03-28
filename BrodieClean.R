### Brodie Germ Exp - Start Clean ###

##### Libs #####
library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)
library(NatParksPalettes)

#############Germination Probability for ALL Species Combined############
dat <- read.csv("Allelepathy and Germination (R_Setup).csv")

dat$Treatment <- as.factor(dat$Treatment)
dat$Dish. <- as.factor(dat$Dish.)
dat$Species <- as.factor(dat$Species)
dat$DayGerm <- as.numeric(dat$DayGerm)

counts <- dat %>% 
  group_by(Species,Treatment,Dish.) %>% 
  summarise(nonzero_germ = sum(DayGerm !=0),
            .groups = "drop") %>% 
  mutate(total = 4)

joined <- left_join(dat, counts, by = c("Species", "Treatment", "Dish."))

df_ <- joined %>%
  group_by(Species, Treatment) %>%
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data = df_,
  family = binomial
)

emm <- emmeans(m_glm, ~ Treatment)
#

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

yellowstone_cols <- NatParksPalettes::natparks.pals("Charmonix", n_treat)

yellowstone_cols
length(yellowstone_cols)

names(yellowstone_cols) <- levels(newdat$Treatment)

ggplot(newdat, aes(x = Treatment, y = pred_prob, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = .8), size = 5, color = "black") +
  scale_color_manual(values = yellowstone_cols, drop = FALSE) +
  labs(
    x = "Treatment",
    y = "Predicted germination probability for all species",
    title = "Fixed-effect treatment estimates"
  ) +
  theme_bw()

#############Germination Probability for RB ##########
dat <- read.csv("Allelepathy and Germination (R_Setup).csv")

dat <- dat %>% 
  filter(Species == "RB")
  
dat$Treatment <- as.factor(dat$Treatment)
dat$Dish. <- as.factor(dat$Dish.)
dat$Species <- as.factor(dat$Species)
dat$DayGerm <- as.numeric(dat$DayGerm)

counts <- dat %>% 
  group_by(Species,Treatment,Dish.) %>% 
  summarise(nonzero_germ = sum(DayGerm !=0),
            .groups = "drop") %>% 
  mutate(total = 4)

joined <- left_join(dat, counts, by = c("Species", "Treatment", "Dish."))

df_ <- joined %>%
  group_by(Species, Treatment) %>%
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data = df_,
  family = binomial
)

emm <- emmeans(m_glm, ~ Treatment)
#

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

yellowstone_cols <- NatParksPalettes::natparks.pals("Charmonix", n_treat)

yellowstone_cols
length(yellowstone_cols)

names(yellowstone_cols) <- levels(newdat$Treatment)

ggplot(newdat, aes(x = Treatment, y = pred_prob, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = .8), size = 5, color = "black") +
  scale_color_manual(values = yellowstone_cols, drop = FALSE) +
  labs(
    x = "Treatment",
    y = "Predicted germination probability for Roundhead Bushclover",
    title = "Fixed-effect treatment estimates"
  ) +
  theme_bw()

################# Germination Probability for BS############

dat <- read.csv("Allelepathy and Germination (R_Setup).csv")

dat <- dat %>% 
  filter(Species == "BS")

dat$Treatment <- as.factor(dat$Treatment)
dat$Dish. <- as.factor(dat$Dish.)
dat$Species <- as.factor(dat$Species)
dat$DayGerm <- as.numeric(dat$DayGerm)

counts <- dat %>% 
  group_by(Species,Treatment,Dish.) %>% 
  summarise(nonzero_germ = sum(DayGerm !=0),
            .groups = "drop") %>% 
  mutate(total = 4)

joined <- left_join(dat, counts, by = c("Species", "Treatment", "Dish."))

df_ <- joined %>%
  group_by(Species, Treatment) %>%
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data = df_,
  family = binomial
)

emm <- emmeans(m_glm, ~ Treatment)
#

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

yellowstone_cols <- NatParksPalettes::natparks.pals("Charmonix", n_treat)

yellowstone_cols
length(yellowstone_cols)

names(yellowstone_cols) <- levels(newdat$Treatment)

ggplot(newdat, aes(x = Treatment, y = pred_prob, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = .6), size = 5, color = "black") +
  scale_color_manual(values = yellowstone_cols, drop = FALSE) +
  labs(
    x = "Treatment",
    y = "Predicted germination probability for Broom Sedge",
    title = "Fixed-effect treatment estimates"
  ) +
  theme_bw()

############### Germination Probability for PB##########
dat <- read.csv("Allelepathy and Germination (R_Setup).csv")

dat <- dat %>% 
  filter(Species == "PB")

dat$Treatment <- as.factor(dat$Treatment)
dat$Dish. <- as.factor(dat$Dish.)
dat$Species <- as.factor(dat$Species)
dat$DayGerm <- as.numeric(dat$DayGerm)

counts <- dat %>% 
  group_by(Species,Treatment,Dish.) %>% 
  summarise(nonzero_germ = sum(DayGerm !=0),
            .groups = "drop") %>% 
  mutate(total = 4)

joined <- left_join(dat, counts, by = c("Species", "Treatment", "Dish."))

df_ <- joined %>%
  group_by(Species, Treatment) %>%
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data = df_,
  family = binomial
)

emm <- emmeans(m_glm, ~ Treatment)
#

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

yellowstone_cols <- NatParksPalettes::natparks.pals("Charmonix", n_treat)

yellowstone_cols
length(yellowstone_cols)

names(yellowstone_cols) <- levels(newdat$Treatment)

ggplot(newdat, aes(x = Treatment, y = pred_prob, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = .8), size = 5, color = "black") +
  scale_color_manual(values = yellowstone_cols, drop = FALSE) +
  labs(
    x = "Treatment",
    y = "Predicted germination probability for Prairie Blazing Star",
    title = "Fixed-effect treatment estimates"
  ) +
  theme_bw()

####### Germination Probability for LB############
dat <- read.csv("Allelepathy and Germination (R_Setup).csv")

dat <- dat %>% 
  filter(Species == "LB")

dat$Treatment <- as.factor(dat$Treatment)
dat$Dish. <- as.factor(dat$Dish.)
dat$Species <- as.factor(dat$Species)
dat$DayGerm <- as.numeric(dat$DayGerm)

counts <- dat %>% 
  group_by(Species,Treatment,Dish.) %>% 
  summarise(nonzero_germ = sum(DayGerm !=0),
            .groups = "drop") %>% 
  mutate(total = 4)

joined <- left_join(dat, counts, by = c("Species", "Treatment", "Dish."))

df_ <- joined %>%
  group_by(Species, Treatment) %>%
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data = df_,
  family = binomial
)

emm <- emmeans(m_glm, ~ Treatment)
#

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

yellowstone_cols <- NatParksPalettes::natparks.pals("Charmonix", n_treat)

yellowstone_cols
length(yellowstone_cols)

names(yellowstone_cols) <- levels(newdat$Treatment)

ggplot(newdat, aes(x = Treatment, y = pred_prob, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  geom_text(aes(label = .group, y = 1), size = 5, color = "black") +
  scale_color_manual(values = yellowstone_cols, drop = FALSE) +
  labs(
    x = "Treatment",
    y = "Predicted germination probability for Little Blue Stem",
    title = "Fixed-effect treatment estimates"
  ) +
  theme_bw()
