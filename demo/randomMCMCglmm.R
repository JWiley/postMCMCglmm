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

# set seed and estimate model
set.seed(10)
m <- MCMCglmm(cyl ~ hp, random = ~ ID, family = "ordinal",
  data = dat, prior = list(
  B = list(mu = c(0, 0), V = diag(2) * 1e10),
  R = list(V = 1, fix = 1),
  G = list(G1 = list(V = 1, nu = .002))), pr=TRUE,
  nitt = 55000, thin = 100, burnin = 5000, verbose=FALSE)

# rescale back to standard normal
# based on 1 for standard normal plus fixed residual variance
m.random <- m
m.random$Sol <- m.random$Sol/sqrt(1 + 1)
m.random$CP <- m.random$CP/sqrt(1 + 1)
m.random$VCV <- m.random$VCV / 2

# view a summary
summary(m.random)
