#a
cases <- c(12, 14, 33, 50, 67, 74, 123, 141, 165, 204, 253, 246, 240, 246, 232)
year  <- 1:15 + 1980

plot(year, cases, pch=19, xlab="Year", ylab="New AIDS cases",
     main="Figure 1. New AIDS cases by year")


#b
fit1 <- glm(cases ~ year, family=poisson)

plot(year, cases, pch=19, xlab="Year", ylab="New AIDS cases",
     main="Figure 2. Poisson fit: linear time trend")
lines(year, fitted(fit1), lwd=2)


#c
par(mfrow=c(2,2))
plot(fit1, which=1:4)
par(mfrow=c(1,1))


#d
fit2 <- glm(cases ~ year + I(year^2), family=poisson)

plot(year, cases, pch=19, xlab="Year", ylab="New AIDS cases",
     main="Figure 3. Poisson fit: quadratic time trend")
ord <- order(year)
lines(year[ord], fitted(fit2)[ord], lwd=2)


#e
anova(fit1, fit2, test="Chisq")


#2
#a
stand01 <- function(x) (x - min(x)) / (max(x) - min(x))

cereal2 <- data.frame(
  Shelf  = factor(cereal$Shelf),   # treat as categorical for multinomial
  sugar  = stand01(cereal$sugar_g  / cereal$size_g),
  fat    = stand01(cereal$fat_g    / cereal$size_g),
  sodium = stand01(cereal$sodium_mg/ cereal$size_g)
)


#b
par(mfrow=c(1,3))

boxplot(sugar ~ Shelf, data=cereal2, ylab="Sugar (scaled)", xlab="Shelf",
        main="Figure 4. Sugar by shelf", pars=list(outpch=NA))
stripchart(sugar ~ Shelf, data=cereal2, method="jitter", vertical=TRUE,
           pch=1, add=TRUE)

boxplot(fat ~ Shelf, data=cereal2, ylab="Fat (scaled)", xlab="Shelf",
        main="Figure 5. Fat by shelf", pars=list(outpch=NA))
stripchart(fat ~ Shelf, data=cereal2, method="jitter", vertical=TRUE,
           pch=1, add=TRUE)

boxplot(sodium ~ Shelf, data=cereal2, ylab="Sodium (scaled)", xlab="Shelf",
        main="Figure 6. Sodium by shelf", pars=list(outpch=NA))
stripchart(sodium ~ Shelf, data=cereal2, method="jitter", vertical=TRUE,
           pch=1, add=TRUE)

par(mfrow=c(1,1))


#d
library(nnet)

fit_full <- multinom(Shelf ~ sugar + fat + sodium, data=cereal2, trace=FALSE)

# LOVO: remove one predictor at a time and compare via AIC (quick) or CV log-loss (better)
fit_no_sugar  <- multinom(Shelf ~ fat + sodium, data=cereal2, trace=FALSE)
fit_no_fat    <- multinom(Shelf ~ sugar + sodium, data=cereal2, trace=FALSE)
fit_no_sodium <- multinom(Shelf ~ sugar + fat, data=cereal2, trace=FALSE)

AIC(fit_full, fit_no_sugar, fit_no_fat, fit_no_sodium)


# -----------------------------
# LOOCV log-loss for multinomial models
# -----------------------------

logloss <- function(probs, true_class) {
  eps <- 1e-15
  probs <- pmax(pmin(probs, 1 - eps), eps)
  -log(probs[true_class])
}

cv_multinom <- function(formula, data) {
  
  y <- data$Shelf
  n <- nrow(data)
  ll <- numeric(n)
  
  for(i in 1:n){
    
    train <- data[-i, ]
    test  <- data[i, , drop = FALSE]
    
    fit <- multinom(formula, data = train, trace = FALSE)
    
    # predicted probabilities
    p <- predict(fit, newdata = test, type = "probs")
    p <- as.matrix(p)
    
    # create probability vector for ALL shelf levels
    full_p <- rep(0, length(levels(y)))
    names(full_p) <- levels(y)
    
    # fill probabilities returned by model
    full_p[colnames(p)] <- p
    
    # compute log-loss
    ll[i] <- logloss(full_p, as.character(test$Shelf))
  }
  
  mean(ll)
}

# -----------------------------
# Run LOVO comparisons
# -----------------------------

cv_full      <- cv_multinom(Shelf ~ sugar + fat + sodium, cereal2)
cv_no_sugar  <- cv_multinom(Shelf ~ fat + sodium, cereal2)
cv_no_fat    <- cv_multinom(Shelf ~ sugar + sodium, cereal2)
cv_no_sodium <- cv_multinom(Shelf ~ sugar + fat, cereal2)

c(
  Full_Model = cv_full,
  No_Sugar   = cv_no_sugar,
  No_Fat     = cv_no_fat,
  No_Sodium  = cv_no_sodium
)


#e
summary(fit_full)

# Effect plot: predicted probs vs sugar, holding others at typical values
library(pdp)
pfun <- function(object, newdata){
  probs <- predict(object, newdata=newdata, type="probs")
  colMeans(probs)
}
pd <- partial(fit_full, pred.var="sugar", pred.fun=pfun, plot=FALSE,
              train=cereal2[,c("sugar","fat","sodium","Shelf")])

lattice::xyplot(yhat ~ sugar | yhat.id, data=pd, type="l",
                ylab="Predicted probability",
                main="Figure 7. Effect of sugar on shelf probabilities")


#f
# store mins/maxs from the per-gram variables BEFORE scaling
sugar_pg  <- cereal$sugar_g  / cereal$size_g
fat_pg    <- cereal$fat_g    / cereal$size_g
sodium_pg <- cereal$sodium_mg/ cereal$size_g

mins <- c(sugar=min(sugar_pg), fat=min(fat_pg), sodium=min(sodium_pg))
maxs <- c(sugar=max(sugar_pg), fat=max(fat_pg), sodium=max(sodium_pg))

scale01 <- function(x, vname) (x - mins[vname])/(maxs[vname]-mins[vname])

apple <- data.frame(
  sugar  = scale01(12/28, "sugar"),
  fat    = scale01(0.5/28, "fat"),
  sodium = scale01(130/28, "sodium")
)

predict(fit_full, newdata=apple, type="probs")


#3
#b
fit_full <- glm(cbind(O.ring, Number - O.ring) ~ Temp + Pressure,
                data=chall, family=binomial)
summary(fit_full)


#c
fit_temp <- glm(cbind(O.ring, Number - O.ring) ~ Temp,
                data=chall, family=binomial)

anova(fit_temp, fit_full, test="Chisq")


#e
predict(fit_temp, newdata=data.frame(Temp=31), type="response")


#f
fit_temp_probit <- glm(cbind(O.ring, Number - O.ring) ~ Temp,
                       data=chall, family=binomial(link="probit"))

c(
  logit  = predict(fit_temp, newdata=data.frame(Temp=31), type="response"),
  probit = predict(fit_temp_probit, newdata=data.frame(Temp=31), type="response")
)


#4

install.packages("broom")
#a
# =============================
# Question 4 (ISLR2::Bikeshare)
# Complete code to answer ALL parts
# =============================

rm(list = ls())

# ---- Packages ----
library(ISLR2)
library(MASS)      # glm.nb
library(ggplot2)   # plots
library(broom)     # tidy model outputs
library(dplyr)     # data wrangling

# ---- Load data ----
bike <- ISLR2::Bikeshare

# ---- Ensure factor encoding for categorical predictors ----
bike$workingday <- factor(bike$workingday)
bike$weathersit <- factor(bike$weathersit)
bike$mnth       <- factor(bike$mnth)
bike$hr         <- factor(bike$hr)

# (Optional) check structure
str(bike)

# =============================
# 1) Fit Poisson model + check overdispersion
# =============================
fit_pois <- glm(bikers ~ workingday + temp + weathersit + mnth + hr,
                data = bike, family = poisson)

disp <- sum(residuals(fit_pois, type = "pearson")^2) / df.residual(fit_pois)
disp

# If disp >> 1, Poisson is overdispersed -> Negative Binomial is appropriate.

# (Optional) quick residual checks
# par(mfrow=c(2,2)); plot(fit_pois)

# =============================
# 2) Fit Negative Binomial model (recommended)
# =============================
fit_nb <- glm.nb(bikers ~ workingday + temp + weathersit + mnth + hr,
                 data = bike, control = glm.control(maxit = 200))

summary(fit_nb)

# =============================
# 3) Well-formatted regression table (coefficients + SE)
#    (Use this table in your report; DON'T paste raw R output)
# =============================
tab_nb <- broom::tidy(fit_nb) %>%
  mutate(
    IRR = exp(estimate),
    IRR_low = exp(estimate - 1.96 * std.error),
    IRR_high = exp(estimate + 1.96 * std.error)
  )

tab_nb

# If you want a cleaner table:
tab_nb_select <- tab_nb %>%
  select(term, estimate, std.error, statistic, p.value, IRR, IRR_low, IRR_high)

tab_nb_select

# =============================
# 4) Compare Poisson vs NB (AIC) + show why NB is better
# =============================
AIC(fit_pois, fit_nb)

# =============================
# 5) Variable importance (practical approach)
#    Use ΔAIC from dropping each predictor group
# =============================

# Drop-one style comparisons (nested models)
fit_nb_no_workingday <- update(fit_nb, . ~ . - workingday)
fit_nb_no_temp       <- update(fit_nb, . ~ . - temp)
fit_nb_no_weather     <- update(fit_nb, . ~ . - weathersit)
fit_nb_no_mnth        <- update(fit_nb, . ~ . - mnth)
fit_nb_no_hr          <- update(fit_nb, . ~ . - hr)

AIC_full <- AIC(fit_nb)

imp_aic <- data.frame(
  Model = c("Full", "No workingday", "No temp", "No weathersit", "No mnth", "No hr"),
  AIC   = c(AIC_full,
            AIC(fit_nb_no_workingday),
            AIC(fit_nb_no_temp),
            AIC(fit_nb_no_weather),
            AIC(fit_nb_no_mnth),
            AIC(fit_nb_no_hr))
) %>%
  mutate(DeltaAIC = AIC - min(AIC))

imp_aic

# Bigger AIC increase when dropped => more important predictor group.

# =============================
# 6) Effect plots for TWO most important predictors
#    (Recommended candidates: hr and temp)
#    These are stakeholder-friendly.
# =============================

# Helper: build prediction grid with typical reference values
# Reference levels: first level of each factor by default.
ref_work <- levels(bike$workingday)[1]
ref_weat <- levels(bike$weathersit)[1]
ref_mnth <- levels(bike$mnth)[1]
ref_hr   <- levels(bike$hr)[1]

# ---- Effect plot: TEMP ----
temp_grid <- data.frame(
  temp = seq(min(bike$temp, na.rm=TRUE), max(bike$temp, na.rm=TRUE), length.out = 200),
  workingday = factor(ref_work, levels=levels(bike$workingday)),
  weathersit = factor(ref_weat, levels=levels(bike$weathersit)),
  mnth = factor(ref_mnth, levels=levels(bike$mnth)),
  hr = factor(ref_hr, levels=levels(bike$hr))
)

temp_grid$pred <- predict(fit_nb, newdata = temp_grid, type = "response")

p_temp <- ggplot(temp_grid, aes(x = temp, y = pred)) +
  geom_line() +
  labs(
    title = "Effect of temperature on expected bike rentals",
    x = "Temperature (temp)",
    y = "Predicted mean rentals (bikers)"
  )

p_temp

# ---- Effect plot: HOUR (hr) ----
hr_grid <- data.frame(
  hr = factor(levels(bike$hr), levels = levels(bike$hr)),
  temp = median(bike$temp, na.rm=TRUE),
  workingday = factor(ref_work, levels=levels(bike$workingday)),
  weathersit = factor(ref_weat, levels=levels(bike$weathersit)),
  mnth = factor(ref_mnth, levels=levels(bike$mnth))
)

hr_grid$pred <- predict(fit_nb, newdata = hr_grid, type = "response")

p_hr <- ggplot(hr_grid, aes(x = hr, y = pred, group = 1)) +
  geom_line() +
  labs(
    title = "Effect of hour of day on expected bike rentals",
    x = "Hour of day (hr)",
    y = "Predicted mean rentals (bikers)"
  )

p_hr

# =============================
# 7) Conditions where rentals are highest
#    Use model-based predictions for a few scenarios
# =============================

# Build a scenario function
scenario_pred <- function(workingday_level, weathersit_level, mnth_level, hr_level, temp_value) {
  nd <- data.frame(
    workingday = factor(workingday_level, levels = levels(bike$workingday)),
    weathersit = factor(weathersit_level, levels = levels(bike$weathersit)),
    mnth = factor(mnth_level, levels = levels(bike$mnth)),
    hr = factor(hr_level, levels = levels(bike$hr)),
    temp = temp_value
  )
  predict(fit_nb, newdata = nd, type = "response")
}

# Example scenarios (edit levels after you check levels(bike$...))
levels(bike$workingday)
levels(bike$weathersit)
levels(bike$mnth)
levels(bike$hr)

# Choose "warm temp" near upper quartile
warm_temp <- quantile(bike$temp, 0.75, na.rm = TRUE)

# Example predictions (you may need to adjust factor labels based on your levels)
pred1 <- scenario_pred(levels(bike$workingday)[1], levels(bike$weathersit)[1],
                       levels(bike$mnth)[7], levels(bike$hr)[18], warm_temp)

pred2 <- scenario_pred(levels(bike$workingday)[2], levels(bike$weathersit)[1],
                       levels(bike$mnth)[7], levels(bike$hr)[9], warm_temp)

pred1; pred2

# =============================
# 8) Optional: Check for zero-inflation evidence
#    (Simple diagnostic: proportion of zeros vs model-implied zeros)
# =============================

prop_zeros <- mean(bike$bikers == 0, na.rm = TRUE)
prop_zeros

# Expected zeros under NB depends on mu and theta; compute average implied P(Y=0)
mu_hat <- predict(fit_nb, type="response")
theta  <- fit_nb$theta
p0_nb  <- mean((theta / (theta + mu_hat))^theta)
p0_nb

# If observed zeros >> implied zeros, that suggests possible zero inflation.

# =============================
# 9) Save plots (for your report)
# =============================
ggsave("Q4_effect_temp.png", p_temp, width = 7, height = 4, dpi = 300)
ggsave("Q4_effect_hr.png",   p_hr,   width = 7, height = 4, dpi = 300)

# End of Question 4 code


