# Author: Gordon Burtch
# Subject: Simulating FE estimations with 0 outcome panels.

library(splitstackshape)
library(VGAM)
library(plm)
library(pglm)
library(texreg)

set.seed(1001)

# Function to take output of pglm and be able to use it with texreg
extract.pglm <- function (model, include.nobs = TRUE, include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$estimate)
  coefficients <- s$estimate[, 1]
  standard.errors <- s$estimate[, 2]
  significance <- s$estimate[, 4]
  loglik.value <- s$loglik
  n <- nrow(model$model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, loglik.value)
    gof.names <- c(gof.names, "Log-Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                     se = standard.errors, pvalues = significance, gof.names = gof.names, 
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)
}

# Let's build our panel dataset, such that we have a fixed subject ID and static feature (weight).
# We won't use the latter here; it's just a second feature that makes it easier to expand the dataframe.
weight <- rnorm(500, mean = 180, sd=30)
id <- factor(seq(1:500))
df <- data.frame(cbind(id,weight))
df <- expandRows(df, count=4,count.is.col = FALSE)
df$period <- rep(seq.int(nrow(df)/4),4)
df$period <- (df$period-1)%%4
df <- df[order(df$id,df$period),]
df$Treat <- round(runif(2000,0,1),0)

# Let's simulate the outcome. I'm using a zero-inflated Poisson distribution, to naturally obtain several 0 outcome panels.
df$Y <- rzipois(n = nrow(df), lambda = exp(0.1 + 0.6*df$Treat), p = 0.7)

# Let's ID panels that have variation in them (some non 0's)
ids_no0s <- df %>% group_by(id) %>% summarize(min_Y = min(Y), max_Y = max(Y)) %>% filter(min_Y!=0|max_Y!=0) %>% pull(id)
df_no0s <- df[df$id %in% ids_no0s,]

# Now we will estimate a conditional FE Poisson with all data, a dummy FE Poisson with all data, and a dummy FE removing the 0's.
poisson_all <- extract.pglm(pglm(data=df, Y~Treat, index=c("id","period"), model = "within", family="poisson"))
poisson_all_d <- glm(data=df, Y~Treat + factor(id) + factor(period), family="poisson")
poisson_all_d_no0s <- glm(data = df_no0s, Y~Treat + factor(id) + factor(period), family="poisson")
screenreg(list(poisson_all, poisson_all_d, poisson_all_d_no0s), omit.coef="factor", custom.model.names = c("Cond FE","Dummy FE","Dummy FE (Remove 0's)"))

