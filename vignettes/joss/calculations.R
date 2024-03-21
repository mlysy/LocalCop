## ---- libs

library(LocalCop)   # local likelihood estimation
library(VineCopula) # simulate copula data

## ---- dgm1-setup

set.seed(2024)

# simulation setting
family <- 3                    # Clayton Copula
n_obs <- 300                   # number of observations
eta_fun <- function(x) {       # calibration function
  sin(5*pi*x) + cos(8*pi*x^2)
}

## ---- dgm2-setup

set.seed(2024)

# simulation setting
family <- 3                    # Clayton Copula
n_obs <- 1000                  # number of observations
eta_fun <- function(x) {       # calibration function
  sin(5*pi*x) + cos(8*pi*x^2)
}


## ---- dgm

# simulate covariate values
x <- sort(runif(n_obs))

# simulate response data
eta_true <- eta_fun(x)                     # calibration parameter eta(x)
par_true <- BiCopEta2Par(family = family,  # copula parameter theta(x)
                         eta = eta_true)
udata <- VineCopula::BiCopSim(n_obs, family = family, par = par_true)

## ---- dgm-plot

# plot Kendall tau
tibble(
  x = x,
  tau = VineCopula::BiCopPar2Tau(family, par = par_true)
) %>%
  ggplot(aes(x = x, y = tau)) +
  geom_line() +
  ylim(c(0, 1)) +
  xlab(expression(x)) + ylab(expression(tau(x)))

## ---- select-precalc

# model selection and tuning
bandset <- c(.02, .05, .1, .2) # set of bandwidth parameters
famset <- c(1, 2, 3, 4, 5)     # set of copula families
kernel <- KernGaus             # kernel function
degree <- 1                    # degree of local polynomial
n_loo <- 100                   # number of LOO-CV observations
                               # (can be much smaller than n_obs)

## ---- select-calc

# calculate cv for each combination of family and bandwidth
cvselect <- CondiCopSelect(u1= udata[,1], u2 = udata[,2],
                           x = x, xind = n_loo,
                           kernel = kernel, degree = degree,
                           family = famset, band = bandset)

## ---- select-save

saveRDS(cvselect, file = "cvselect.rds")

## ---- select-load

cvselect <- readRDS("cvselect.rds")

## ---- select-plot

# plot cv results
fam_names <- c("Gaussian", "Student", "Clayton", "Gumbel", "Frank")
as_tibble(cvselect$cv) %>%
  mutate(
    family = factor(family, levels = c(1,2,3,4,5),
                    labels = fam_names),
    Bandwidth = factor(band),
  ) %>%
  ggplot(aes(x = family, y  = cv, fill = Bandwidth)) +
  geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("red", "green4", "blue", "orange")) +
  ## scale_fill_brewer(palette="Blues", direction=-1) +
  xlab("") + ylab("CV Likelihood") +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

## ---- locfit

# extract the selected family and bandwidth from cvselect
cv_res <- cvselect$cv
i_opt <- which.max(cv_res$cv)
fam_opt <- cv_res[i_opt,]$family
band_opt <- cv_res[i_opt,]$band

# calculate eta(x) on a grid of values
x0 <- seq(0, 1, by = 0.01)
copfit <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                         x = x, x0 = x0,
                         kernel = kernel, degree = degree,
                         family = fam_opt, band = band_opt)
# convert eta to Kendall tau
tau_loc <- BiCopEta2Tau(copfit$eta, family= fam_opt)

## ---- gamfit

# fit with gamCopula
gam_fit <- gamCopula::gamBiCopSelect(
  udata = udata[,1:2],
  smooth.covs = data.frame(x = x)
)
gam_res <- gam_fit$res
gam_pred <- gamCopula::gamBiCopPredict(
  object = gam_res,
  newdata = data.frame(x = x0),
  target = "tau"
)
tau_gam <- as.numeric(gam_pred$tau)

## ---- condfit

# fit with CondCopulas
cond_select <- CondCopulas::CKT.hCV.l1out(
  observedX1 = udata[,1], observedX2 = udata[,2],
  observedZ = x,
  # ZToEstimate = x0,
  range_h = bandset,
)

cond_par <- CondCopulas::estimateParCondCopula(
  observedX1 = udata[,1], observedX2 = udata[,2],
  observedX3 = x, newX3 = x0,
  family = fam_opt,
  h = cond_select$hCV,
  method = "mle"
)
tau_cond <- VineCopula::BiCopPar2Tau(
  family = fam_opt,
  par = cond_par
)

## ---- plotfit

plotfit <- function(dat) {
  dat %>%
    pivot_longer(cols = !x:True, names_to = "Estimator", values_to = "Est") %>%
    pivot_longer(cols = c("True", "Est"), names_to = "Method", values_to = "tau") %>%
    mutate(
      Method = ifelse(Method == "Est", Estimator, Method),
      Method = factor(Method, levels = c("True", "LocalCop", "gamCopula", "CondCopulas")),
      Estimator = factor(Estimator, levels = c("LocalCop", "gamCopula", "CondCopulas"))
    ) %>%
    ggplot(aes(x = x, y = tau, group = Method)) +
    geom_line(aes(color = Method), linewidth = 1) +
    scale_color_manual(
      name = expression("Kendall "*tau),
      values = c("black", "red", "green4", "blue")
    ) +
    facet_wrap(. ~ Estimator) +
    xlab(expression(x)) + ylab(expression(tau(x))) +
    theme_minimal() +
    theme(
      legend.position = "right",
      ## legend.title = element_blank(),
      strip.text = element_blank()
    )
}

## ---- scratch

xi <- .26
wgt <- KernWeight(x = x, x0 = xi, band = band_opt,
                  kernel = KernEpa, band_type = "constant")
adfun <- CondiCopLocFun(
  u1 = udata[,1],
  u2 = udata[,2],
  wgt = wgt,
  family = fam_opt,
  x = x,
  x0 = xi,
  degree = 0,
  eta = c(0, 0),
  nu = 1
)


optim_fun(adfun)

CondiCopLocFit(
  u1 = udata[,1],
  u2 = udata[,2],
  family = fam_opt,
  x = x,
  x0 = xi,
  degree = 0,
  band = band_opt,
  optim_fun = optim_fun
)

optim(par = c(.00711, -9.2171), fn = adfun$fn, method = "BFGS")


