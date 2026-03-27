# Brodie Germ Exp

dat <- read.csv("Allelepathy and Germination (R_Setup).csv")

seeds_per_dish <- 4

# Create a germination proportion and a cbind(success, failure) column
dat$germinated <- dat$DayGerm
dat$total_seeds <- seeds_per_dish
dat$not_germinated <- dat$total_seeds - dat$germinated
  
counts <- dat %>%
  group_by(Species,Treatment,Dish.) %>%
  summarise(nonzero_germ = sum(DayGerm != 0),
            .groups = "drop")
  
  
dat$Dish. <- as.numeric(dat$Dish.)
counts$Dish. <- as.numeric(counts$Dish.)

joined <- left_join(dat, counts, by = c("Species","Treatment","Dish."))

dat$Species <- as.factor(dat$Species)
dat$Treatment <- as.factor(dat$Treatment)


# Fit GLM separately for each species
library(dplyr)

glm_by_species <- joined %>%
  group_by(Species) %>%
  do(
    model = glm(
      cbind(nonzero_germ, seeds_per_dish - nonzero_germ) ~ Treatment,
      family = binomial,
      data = joined
    )
  )

m_glm <- glm(
  cbind(nonzero_germ, seeds_per_dish - nonzero_germ) ~ Treatment,
  data   = joined,
  family = binomial
)

summary(m_glm)
anova(m_glm, test = "Chisq")

# See model summaries
summary(glm_by_species$model[[1]])  # first species
summary(glm_by_species$model[[2]])  # second species, etc.
summary(glm_by_species$model[[3]])
summary(glm_by_species$model[[4]])
# If instead you want a single model with Species and Treatment plus interaction:
glm_all <- glm(
  cbind(germinated, not_germinated) ~ Species * Treatment,
  family = binomial,
  data = dat
)
summary(glm_all)
