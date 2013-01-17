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
m.full <- MCMCglmm(cyl ~ qsec + mpg + drat, random = ~ ID, family = "ordinal",
  data = dat, prior = list(
  B = list(mu = c(0, 0, 0, 0), V = diag(4) * 1e2),
  R = list(V = 1, fix = 1),
  G = list(G1 = list(V = 1, nu = .002))), pr=TRUE,
  nitt = 55000, thin = 20, burnin = 5000, verbose=FALSE)

# rescale back to standard normal
# based on 1 for standard normal plus fixed residual variance
m.full$Sol <- m.full$Sol/sqrt(1 + 1)
m.full$CP <- m.full$CP/sqrt(1 + 1)
m.full$VCV <- m.full$VCV / 2


## model summary
summary(m.full)

## predictor means
colMeans(dat[, c("qsec", "mpg", "drat")])
## predictor ranges
lapply(dat[, c("qsec", "mpg", "drat")], range)

## new data for prediction
myX <- as.matrix(data.frame("(Intercept)" = 1, qsec = 0, mpg = 0, drat = seq(from = 2.76, to = 4.93, by = .1), check.names=FALSE))

## get the predicted values (fixed effects only)
pred1 <- predict2(m.full, X = myX, Z = NULL, use = "all", type = "response", varepsilon = 1)

## summarize them
(spred1 <- summary(pred1))

## can combine predicted probs + HPD intervals with prediction data
## to be able to plot
preddat <- as.data.frame(cbind(do.call(rbind, rep(list(myX), 3)), do.call(rbind, spred1)))
## need an indicator for which level of the outcome
preddat$outcome <- factor(rep(1:3, each = nrow(myX)))

## look at the data
head(preddat)

## plot it
require(ggplot2)
ggplot(preddat, aes(x = drat, y = M, colour = outcome)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = outcome), alpha = .25) +
  geom_line(size=2)


## now get average marginal recycled probabilities
## note that this can be a bit slow
tmpres <- lapply(c("qsec", "mpg", "drat"), function(i) {
  tmp <- recycler(m, index = i, twiddle = .01)
  tmp <- summary(tmp, marginalize = TRUE)
  tmp <- do.call("rbind", tmp)
  (tmp2 <- sprintf("%0.2f [%0.2f, %0.2f]", tmp[, 1], tmp[, 2], tmp[, 3]))
})

## setup final "output" table
finaltable <- matrix(NA, nrow = 3, ncol = 4, dimnames = list(
  NULL, c("Variable", "L1", "L2", "L3")))

finaltable[1, ] <- c("Var1", tmpres[[1]])
finaltable[2, ] <- c("Var2", tmpres[[2]])
finaltable[3, ] <- c("Var3", tmpres[[3]])
finaltable[,1] <- c("qsec", "mpg", "disp")
finaltable[] <- gsub("0\\.", ".", finaltable)

## these are the average changes in probability
## of membership in each group, for a one unit change
## of the predictor, on average across the sample
## with 95% HPD intervals
finaltable

## write to clipboard on Windows (e.g., for collaborators who want in Excel)
#write.table(finaltable, file = "clipboard", na = "", sep = "\t", row.names=FALSE, col.names=TRUE)
