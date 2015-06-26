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
m2 <- MCMCglmm(cyl ~ trait + qsec + drat,
  random = ~ ID, rcov = ~ us(trait):units,
  family = "categorical", data = dat,
  prior = list(
    B = list(mu = c(0, 0, 0, 0), V = diag(4) * c(1e4, 1e4, rep(1e2, 2))),
    R = list(V = diag(2), nu = .002, fix = 1),
    G = list(G1 = list(V = 1, nu = 1))
  ), pr=TRUE, pl=TRUE,
  nitt = 55000, thin = 20, burnin = 5000, verbose=TRUE)

# set seed and estimate model
set.seed(10)
m3 <- MCMCglmm(cyl ~ trait + qsec + drat,
  random = ~ ID, rcov = ~ us(trait):units,
  family = "categorical", data = dat,
  prior = list(
    B = list(mu = c(0, 0, 0, 0), V = diag(4) * c(1e4, 1e4, rep(1e2, 2))),
    R = list(V = diag(2), nu = .002, fix = 1),
    G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V=100))
  ), pr=TRUE, pl=TRUE,
  nitt = 55000, thin = 20, burnin = 5000, verbose=TRUE)

