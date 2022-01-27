## -------------------------------------------------------------------------------
data(stamps1, package = "multimode")
str(stamps1)
Thickness <- stamps1$thickness


## ----fig.keep="none"------------------------------------------------------------
dens <- densityMclust(Thickness)


## -------------------------------------------------------------------------------
summary(dens, parameters = TRUE)


## ----stamp1, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Histogram of thickness for the Hidalgo 1872 \\code{stamps1}\\ dataset, with the GMM density estimate superimposed."----
br <- seq(min(Thickness), max(Thickness), length.out = 21)
plot(dens, what = "density", data = Thickness, breaks = br)


## -------------------------------------------------------------------------------
with(dens$parameters, 
     data.frame(mean = mean,
                sd = sqrt(variance$sigmasq),
                CoefVar = sqrt(variance$sigmasq)/mean*100))


## -------------------------------------------------------------------------------
x <- c(0.07, 0.08, 0.1, 0.12)
predict(dens, newdata = x, what = "dens")
predict(dens, newdata = x, logarithm = TRUE)
predict(dens, newdata = x, what = "cdens")


## -------------------------------------------------------------------------------
predict(dens, newdata = x, what = "z")


## ----stamp2, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Histograms of thickness by overprinted year for the Hidalgo 1872 \\code{stamps1}\\ dataset, with mixture component-density estimates superimposed."----
Year <- stamps1$year
table(Year)
h1 <- hist(Thickness[Year == "1872"], breaks = br, plot = FALSE)
h1$density <- h1$density*prop.table(table(Year))[1]
h2 <- hist(Thickness[Year == "1873-1874"], breaks = br, plot = FALSE)
h2$density <- h2$density*prop.table(table(Year))[2]
x <- seq(min(Thickness)-diff(range(Thickness))/10, 
         max(Thickness)+diff(range(Thickness))/10, length = 200)
cdens <- predict(dens, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro))
col <- adjustcolor(mclust.options("classPlotColors")[1:2], alpha = 0.3)
ylim <- range(h1$density, h2$density, cdens)
plot(h1, xlab = "Thickness", freq = FALSE, main = "", border = FALSE, 
     col = col[1], xlim = range(x), ylim = ylim)
plot(h2, add = TRUE, freq = FALSE, border = FALSE, col = col[2])
matplot(x, cdens, type = "l", lwd = 1, lty = 1, col = 1, add = TRUE)
box()
legend("topright", legend = levels(Year), col = col, pch = 15, inset = 0.02,
       title = "Overprinted years:", title.adj = 0.2)


## -------------------------------------------------------------------------------
data(acidity, package = "BNPdensity")


## -------------------------------------------------------------------------------
summary(mclustBIC(acidity), k = 5)


## ----fig.keep="none"------------------------------------------------------------
dens_E2 <- densityMclust(acidity, G = 2, modelNames = "E")
summary(dens_E2, parameters = TRUE)


## ----acidity_LRT, cache=TRUE----------------------------------------------------
mclustBootstrapLRT(acidity, modelName = "V")


## ----fig.keep="none"------------------------------------------------------------
dens_V3 <- densityMclust(acidity, G = 3, modelNames = "V")
summary(dens_V3, parameters = TRUE)


## ----acidity_dens, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Density estimates provided by (a) the model (\\code{E},2) supported by BIC, and (b) the model (\\code{V},3) supported by the LRT, for the \\code{acidity}\\ data."----
plot(dens_E2, what = "density",
     ylim = c(0, max(dens_E2$density, dens_V3$density)))
rug(acidity)
plot(dens_V3, what = "density", 
     ylim = c(0, max(dens_E2$density, dens_V3$density)))
rug(acidity)


## ----echo=FALSE, eval=FALSE-----------------------------------------------------
## x0 <- seq(-2,8,length.out=1000)
## hist(acidity, breaks = 15, probability = TRUE,
##      main = NULL, border = "white", col = "lightgrey")
## box()
## lines(x0, predict(dens_E2, x0), col = "red3")
## lines(x0, predict(dens_V3, x0), col = "dodgerblue2")


## ----acidity_E2_diagn, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Density estimate diagnostics for the (\\code{E},2) model estimated for the \\code{acidity}\\ data: (a) estimated CDF and empirical distribution function, (b) Q-Q plot of sample quantiles vs quantiles from density estimation."----
plot(dens_E2, what = "diagnostic", type = "cdf")
plot(dens_E2, what = "diagnostic", type = "qq")


## ----acidity_V3_diagn, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Density estimate diagnostics for the (\\code{V},3) model estimated for the \\code{acidity}\\ data: (a) estimated CDF and empirical distribution function, (b) Q-Q plot of sample quantiles vs quantiles from density estimation."----
plot(dens_V3, what = "diagnostic", type = "cdf")
plot(dens_V3, what = "diagnostic", type = "qq")


## ----eval=FALSE, echo=FALSE-----------------------------------------------------
## ks.test(acidity, function(x) cdfMclust(dens_E2, data = x)$y)
## ks.test(acidity, function(x) cdfMclust(dens_V3, data = x)$y)


## ----fig.keep="none"------------------------------------------------------------
data(faithful, package = "datasets")
plot(faithful)
dens <- densityMclust(faithful)
summary(dens, parameters = TRUE)


## ----echo=FALSE, out.width="0.49\\textwidth"------------------------------------
plot(faithful)


## ----echo=FALSE, out.width="0.49\\textwidth"------------------------------------
plot(dens, what = "dens")


## ----echo=FALSE, out.width="0.49\\textwidth"------------------------------------
plot(dens, what = "density", type = "image")


## ----echo=FALSE, out.width="0.49\\textwidth"------------------------------------
plot(dens, what = "density", type = "persp")


## ----fig.keep="none"------------------------------------------------------------
plot(dens, what = "density")


## ----fig.keep="none"------------------------------------------------------------
plot(dens, what = "density", type = "image")
plot(dens, what = "density", type = "persp")


## ----eval=FALSE, echo=FALSE-----------------------------------------------------
## plot(dens, what = "density", data = faithful, grid = 200, points.cex = 0.5,
##      drawlabels = FALSE)
## plot(dens, what = "density", type = "hdr", prob = c(0.25, 0.5, 0.75))
## plot(dens, what = "density", type = "image", col = "steelblue", grid = 200)
## plot(dens, what = "density", type = "persp", theta = -25, phi = 20)


## ----faithful_hdr, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Highest density regions from the density estimated on the \\code{faithful}\\ data at probability levels 0.25, 0.5, and 0.75."----
plot(dens, what = "density", type = "hdr")


## -------------------------------------------------------------------------------
data(aircraft, package = "sm")
X <- log(subset(aircraft, subset = (Period == 3), select = 3:8))
PCA <- princomp(X, cor = TRUE)
summary(PCA, loadings = TRUE, cutoff = 0)
Z <- PCA$scores[,1:3]
colnames(Z) <- paste0("PC", 1:ncol(Z))


## -------------------------------------------------------------------------------
BIC <- mclustBIC(Z)
summary(BIC, k = 5)


## ----aircraft, fig.width=8, fig.height=7, out.width="\\textwidth", fig.cap="Scatterplot matrix of selected principal components for the \\code{aircraft}\\ data with bivariate highest density regions at probability levels 0.25, 0.5, and 0.75."----
densAircraft <- densityMclust(Z, G = 5, modelNames = "VVE", plot = FALSE)
plot(densAircraft, what = "density", type = "hdr", 
     data = Z, points.cex = 0.5)


## -------------------------------------------------------------------------------
library("mclustAddons")


## ----fig.keep="none"------------------------------------------------------------
data("suicide", package = "mclustAddons")
dens <- densityMclust(suicide)
rug(suicide)            # add data points at the bottom of the graph
abline(v = 0, lty = 3)  # draw a vertical line at the natural boundary


## ----fig.keep="none"------------------------------------------------------------
bdens <- densityMclustBounded(suicide, lbound = 0)
summary(bdens, parameters = TRUE)
plot(bdens, what = "density")
rug(suicide)            # add data points at the bottom of the graph
abline(v = 0, lty = 3)  # draw a vertical line at the natural boundary


## ----suicide_dens, echo=FALSE, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Density estimates for the \\code{suicide}\\ data. Panel (a) shows the default estimate which ignores the lower boundary of the variable, while panel (b) shows the density estimated accounting for the natural boundary of the variable."----
plot(dens, what = "density")
rug(suicide)
abline(v  =0, lty = 3)
#
plot(bdens, what = "density")
rug(suicide)
abline(v  =0, lty = 3)


## ----racial_dens, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Density estimate for the \\code{racial}\\ data obtained by taking into account the natural boundary of the proportion."----
data("racial", package = "mclustAddons")
bdens <- densityMclustBounded(racial$PropWhite, lbound = 0, ubound = 1)
plot(bdens, what = "density", 
     lwd = 2, col = "dodgerblue2",
     data = racial$PropWhite, breaks = 15,
     xlab = "Proportion of white student enrolled in schools")
rug(racial$PropWhite)        # add data points at the bottom of the graph
abline(v = c(0,1), lty = 3)  # draw a vertical line at the natural boundary


## ----hdr_mixdens, fig.pos="hb", fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Density of the univariate two-component Gaussian mixture $f(x) = 0.7~\\phi(x | \\mu = 0, \\sigma = 1) + 0.3~\\phi(x | \\mu = 4, \\sigma = 1)$."----
f <- function(x) 
  0.7*dnorm(x, mean = 0, sd = 1) + 0.3*dnorm(x, mean = 4, sd = 1)
curve(f, from = -4, to = 8)


## ----echo=-1--------------------------------------------------------------------
set.seed(20190525)
par <- list(pro = c(0.7,0.3), mean = c(0,4), 
            variance = mclustVariance("E", G = 2))
par$variance$sigmasq <- c(1,1)
x <- sim(modelName = "E", parameters = par, n = 1e4)[,2]


## ----hdr_levels, echo=-1, fig.width=6, fig.height=5, out.width="\\textwidth", fig.cap="Highest density regions at specified probability levels from the density of a univariate two-component Gaussian mixture."----
par(mfrow=c(2,2))
prob <- c(0.25, 0.5, 0.75, 0.95)
(hdr <- hdrlevels(f(x), prob))
for (j in seq(prob))
{
  curve(f, from = -4, to = 8)
  mtext(side = 3, paste0(prob[j]*100, "% HDR"), adj = 0)
  abline(h = hdr[j], lty = 2)
  rug(x, col = "lightgrey")
  rug(x[f(x) >= hdr[j]])
}


## ----fig.keep="none"------------------------------------------------------------
dens <- densityMclust(x, plot = FALSE)
hdrlevels(predict(dens, x), prob)

