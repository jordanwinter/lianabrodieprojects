library(tidyverse)
library(survival)


df <- read.csv("GerminationRate_Datasheetcsv.csv")

df <- df[-c(2401,2402),]
max_day <- 32

df$event <- ifelse(df$Germination.Day > 0,1,0)

df$time <- ifelse(df$Germination.Day > 0, df$Germination.Day, max_day)



surv_obj <- Surv(time = df$time, event = df$event)

fit <- survfit(surv_obj ~ Treatment, data = df)

plot(fit, col = 1:length(levels(df$Treatment)))


survdiff(surv_obj ~ Treatment, data = df)  # log-rank test

###############################################################


# read data
dat <- read.csv("GerminationRate_Datasheetcsv.csv")

dat <- dat[-c(2401,2402),]
library(dplyr)

counts <- dat %>%
  group_by(`Plant..`) %>%
  summarise(nonzero_germ = sum(`Germination.Day` != 0),
            .groups = "drop")

counts <- counts %>% 
  mutate(total = 20)


joined <- left_join(dat, counts, by = "Plant..")

df_ <- joined %>%
  group_by(Plant.., Treatment) %>% 
  mutate(BlockID = cur_group_id()) %>%
  ungroup()

m_glm <- glm(
  cbind(nonzero_germ, total - nonzero_germ) ~ Treatment,
  data   = df_,
  family = binomial
)

summary(m_glm)
anova(m_glm, test = "Chisq")  # test treatment effect - it's significant 

library(emmeans)

emm <- emmeans(m_glm, ~ Treatment)

pairs(emm, adjust = "tukey")

library(ggplot2)

df_$Treatment <- as.factor(df_$Treatment) 

# For a GLM with only treatment as predictor:
newdat <- data.frame(
  Treatment = levels(df_$Treatment)
)

newdat$pred_link <- predict(m_glm, newdata = newdat, se.fit = TRUE)$fit
newdat$se_link   <- predict(m_glm, newdata = newdat, se.fit = TRUE)$se.fit

# 95% CI on link (logit) scale
newdat$lwr_link <- newdat$pred_link - 1.96 * newdat$se_link
newdat$upr_link <- newdat$pred_link + 1.96 * newdat$se_link

# back‑transform to probability
ilogit <- function(x) exp(x) / (1 + exp(x))

newdat$pred_prob <- ilogit(newdat$pred_link)
newdat$lwr_prob  <- ilogit(newdat$lwr_link)
newdat$upr_prob  <- ilogit(newdat$upr_link)

ggplot(newdat, aes(x = Treatment, y = pred_prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr_prob, ymax = upr_prob), width = 0.2) +
  labs(x = "Treatment",
       y = "Predicted germination probability",
       title = "Fixed-effect treatment estimates") +
  theme_bw()


