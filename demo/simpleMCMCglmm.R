# sample models I used as test cases for developing the code
set.seed(1234)
dat <- mtcars[sample(1:32, 1000, replace = TRUE), ]
dat <- within(dat, {
  qsec <- scale(qsec)
  hp <- scale(hp)
  mpg <- scale(mpg)
  disp <- scale(disp)
})
dat$ID <- factor(rep(letters, length.out = 1000))
dat$cyl <- factor(dat$cyl)

# estimate model
set.seed(10)
m.simple <- MCMCglmm(cyl ~ qsec, family = "ordinal",
  data = dat, prior = list(
  B = list(mu = c(0, 0), V = diag(2) * 1e10),
  R = list(V = 1, fix = 1)),
  nitt = 55000, thin = 100, burnin = 5000, verbose=FALSE)

# rescale back to standard normal
# based on 1 for standard normal plus fixed residual variance
m.simple$Sol <- m.simple$Sol/sqrt(1 + 1)
m.simple$CP <- m.simple$CP/sqrt(1 + 1)
m.simple$VCV <- m.simple$VCV / 2

# estimate model in classical way
malt <- MASS::polr(cyl ~ qsec, data = dat, method = "probit", Hess=TRUE)


# compare Bayesian to classical estimates
summary(m.simple)
summary(malt)
