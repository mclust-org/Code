## -------------------------------------------------------------------------------
data(diabetes, package = "rrcov")
X <- diabetes[,1:5]
Class <- diabetes$group
table(Class)


## ----diabetes_clpairs, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Pairwise scatterplots for the \\code{diabetes}\\ data with points marked according to the true classification."----
clp <- clPairs(X, Class, lower.panel = NULL)
clPairsLegend(0.1, 0.3, class = clp$class, col = clp$col, pch = clp$pch)


## -------------------------------------------------------------------------------
mod <- Mclust(X, G = 3, modelNames = "VVV")


## -------------------------------------------------------------------------------
summary(mod)


## -------------------------------------------------------------------------------
summary(mod, parameters = TRUE)


## -------------------------------------------------------------------------------
table(Class, Cluster = mod$classification)
adjustedRandIndex(Class, mod$classification)


## ----diabetes_classif, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Scatterplot matrix of variables in the \\code{diabetes}\\ data with points marked according to the \\pkg{mclust} classification, and ellipses corresponding to projections of the \\pkg{mclust} cluster covariances."----
plot(mod, what = "classification")


## ----fig.keep="none"------------------------------------------------------------
plot(mod, what = "classification", fillEllipses = TRUE)


## ----diabetes_classif_fill, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Scatterplot of a pair of variables in the \\code{diabetes}\\ data with points marked according to the \\pkg{mclust} classification, and filled ellipses corresponding to \\pkg{mclust} cluster covariances."----
plot(mod, what = "classification", dimens = c(3,4), fillEllipses = TRUE)


## ----diabetes_uncert, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Scatterplot of a pair of variables in the \\code{diabetes}\\ data with points marked according to the clustering, and point size reflecting the corresponding uncertainty of the MAP classification for the \\pkg{mclust} model."----
plot(mod, dimens = c(3,4), what = "uncertainty")


## ----thyroid_clpairs, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Pairwise scatterplots showing the classification for the \\code{thyroid}\\ gland data."----
data(thyroid, package = "mclust")
X <- data.matrix(thyroid[,2:6])
Class <- thyroid$Diagnosis

clp <- clPairs(X, Class, lower.panel = NULL,
               symbols = c(0,1,2), 
               colors = c("gray50", "black", "red3")) 
clPairsLegend(0.1, 0.3, title = "Thyroid diagnosis:", class = clp$class, 
              col = clp$col, pch = clp$pch)


## -------------------------------------------------------------------------------
mod <- Mclust(X)


## -------------------------------------------------------------------------------
mod$BIC


## -------------------------------------------------------------------------------
summary(mod$BIC, k = 5)


## ----thyroid_bic, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="BIC traces for the GMMs estimated for the \\code{thyroid}\\ data."----
plot(mod, what = "BIC", 
     legendArgs = list("bottomright", ncol = 5))


## -------------------------------------------------------------------------------
summary(mod, parameters = TRUE)


## ----thyroid_classification, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Scatterplot matrix for the \\code{thyroid}\\ data with points marked according to the GMM (\\code{VVI},3) clustering, and ellipses corresponding to projections of the estimated cluster covariances."----
plot(mod, what = "classification")


## -------------------------------------------------------------------------------
table(Class, Cluster = mod$classification)
adjustedRandIndex(Class, mod$classification)


## ----thyroid_postprobs, fig.width=8, fig.height=5, out.width="\\textwidth", fig.cap="Estimated posterior conditional probabilities of class membership for each of the three clusters determined by the MAP classification of observations in the \\code{thyroid}\\ data. The panels correspond to the different clusters."----
z  <- mod$z               # posterior conditional probabilities
cl <- mod$classification  # MAP clustering
G  <- mod$G               # number of clusters
sclass <- 10 # class separation
sedge <- 3   # edge spacing
L <- nrow(z) + G*(sclass+2*sedge)
plot(1:L, runif(L), ylim = c(0,1), type = "n",axes = FALSE, 
     ylab = "Posterior conditional probabilities", xlab = "")
axis(2)
col <- mclust.options("classPlotColors")
l <- sclass
for (k in 1:G)
{
 i <- which(cl == k)
 ord <- i[order(z[i,k], decreasing = TRUE)]
 for (j in 1:G)
    points((l+sedge)+1:length(i), z[ord,j], pch=as.character(j), col=col[j])
 rect(l, 0, l+2*sedge+length(i), 1, 
      border = col[k], col = col[k], lwd = 2, density = 0)
 l <- l + 2*sedge + length(i) + sclass
}


## ----echo=FALSE, eval=FALSE-----------------------------------------------------
## plotCondProbs <- function(z, col = mclust.options("classPlotColors"),
##                           nlevels = 11, ...)
## {
##   z <- as.matrix(z)
##   z <- sweep(z, MARGIN = 1, FUN = "/", STATS = rowSums(z))
##   cl <- map(z)
##   G <- ncol(z)
##   n <- nrow(z)
##   ord <- NULL
##   for (k in seq(G))
##   {
##     i <- which(cl == k)
##     ord <- c(ord, i[order(z[i,k], decreasing = TRUE)])
##   }
##   z <- z[ord,]
##   col <- col[seq(unique(cl))]
##   cnames <- colnames(z)
##   if (is.null(cnames))
##     cnames <- paste("Cluster", seq(1, G))
##   rnames <- rownames(z)
##   if (is.null(rnames))
##     rnames <- seq(n)
##   z1 <- seq(n)
##   z2 <- seq(G)
##   levels <- seq(0, 1, length = nlevels)
##   step <- diff(z1)[1]
##   yrange <- range(z1) + c(-step/2, step/2)
##   step <- diff(z2)[1]
##   xrange <- range(z2) + c(-step/2, step/2)
##   image <- as.data.frame(matrix(NA, nrow = n, G))
##   for (k in seq(G))
##   {
##     image[,k] <- cut(z[,k], breaks = levels, include.lowest = TRUE,
##                      labels = sapply(levels[-1], function(a) adjustcolor(col[k], alpha.f = a)))
##   }
##   image <- as.matrix(image)
## 
##   oldpar <- par(no.readonly = TRUE)
##   on.exit(par(oldpar))
##   par(mar = c(1,5,1,1), oma = c(0,0,1,0), xpd = NA)
##   plot(z, type = "n", xaxt = "n", yaxt = "n", ann = FALSE,
##        xlim = xrange, ylim = yrange, xaxs = "i", yaxs = "i")
##   rasterImage(image, xrange[1], yrange[1], xrange[2], yrange[2],
##               interpolate = FALSE)
##   axis(2, at = seq(n, 1), labels = rnames, cex.axis = 0.8, las = 2)
##   rect(par("usr")[1], par("usr")[4],
##        par("usr")[2], par("usr")[4]+2.5*max(strheight(cnames, units = "user")),
##        col = "lightgrey")
##   mtext(cnames, side = 3, line = 0.25, at = seq(G))
##   par(xpd = FALSE)
##   abline(h = n - rev(cumsum(table(cl)))[-1]+0.5, lty = 2)
## }
## # Almost works, but there is a problem if saved as pdf (smooth colors along the rows) but ok if saved as png
## # TODO: alternatively this function could be used (and added to mclust)
## plotCondProbs(mod$z)


## -------------------------------------------------------------------------------
data(wine, package = "gclus")
Class <- factor(wine$Class, levels = 1:3,
                labels = c("Barolo", "Grignolino", "Barbera"))
X <- data.matrix(wine[,-1])


## ----wine_clust, cache=TRUE-----------------------------------------------------
mod <- Mclust(X)
summary(mod$BIC, k = 3)


## ----wine_bic, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="BIC plot for models fitted to the \\code{wine}\\ data."----
plot(mod, what = "BIC", 
     ylim = range(mod$BIC[,-(1:2)], na.rm = TRUE),
     legendArgs = list(x = "bottomleft"))


## -------------------------------------------------------------------------------
summary(mod)
table(Class, mod$classification)
adjustedRandIndex(Class, mod$classification)


## ----wine_heatmap, fig.width=8, fig.height=5, out.width="\\textwidth", fig.cap="Heatmap of normalized cluster means for the clustering model fitted to the \\code{wine}\\ data."----
norm01 <- function(x) (x - min(x))/(max(x) - min(x))
M <- apply(t(mod$parameters$mean), 2, norm01)
heatmap(M, Rowv = NA, scale = "none", margins = c(8,2),
        labRow = paste("Cluster", 1:mod$G), cexRow = 1.2)


## ----wine_scatterplot, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Coordinate projection plot of selected features showing the clusters for the \\code{wine}\\ data."----
plot(mod, what = "classification", dimens = c(1,2,7))


## ----faithful_data, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Scatterplot of the Old Faithful data available in the \\code{faithful}\\ dataset."----
data(faithful, package = "datasets")
plot(faithful)


## ----faithful_BIC, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="BIC traces for the GMMs estimated on the \\code{faithful}\\ dataset."----
BIC <- mclustBIC(faithful)
BIC
plot(BIC)


## -------------------------------------------------------------------------------
mod1 <- Mclust(faithful, x = BIC)


## -------------------------------------------------------------------------------
ICL <- mclustICL(faithful)
ICL
mod2 <- Mclust(faithful, G = 2, modelNames = "VVE")


## ----faithful_ICL, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="ICL traces for the GMMs estimated for the \\code{faithful}\\ dataset."----
plot(ICL)


## ----faithful_GMMs, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap="", fig.cap="Scatterplots for the \\code{faithful}\\ dataset with points and ellipses corresponding to the classification from the best estimated GMMs selected by BIC (a) and ICL (b)."----
plot(mod1, what = "classification", fillEllipses = TRUE)
plot(mod2, what = "classification", fillEllipses = TRUE)


## ----faithful_lrt, cache=TRUE---------------------------------------------------
LRT <- mclustBootstrapLRT(faithful, modelName = "VVV")
LRT


## ----faithful_lrt_boot, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap="", fig.cap="Histograms of LRTS bootstrap distributions for testing the number of mixture components in the \\code{faithful}\\ data assuming the \\code{VVV}\\ model. The dotted vertical lines refer to the sample values of LRTS."----
plot(LRT, G = 1)
plot(LRT, G = 2)


## ----fig.keep="none"------------------------------------------------------------
data(hemophilia, package = "rrcov")
X <- hemophilia[,1:2]
Class <- as.factor(hemophilia$gr) 
clp <- clPairs(X, Class, symbols = c(16,0), colors = "black")
clPairsLegend(0.8, 0.2, class = clp$class, col = clp$col, pch = clp$pch)


## -------------------------------------------------------------------------------
mod <- Mclust(X, G = 2, modelName = "VVV")
summary(mod, parameters = TRUE)


## ----fig.keep="none"------------------------------------------------------------
plot(mod, what = "classification", fillEllipses = TRUE)


## ----hemophilia, echo=FALSE, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap="", fig.cap="Scatterplots for the \\code{hemophilia}\\ data displaying the true class membership (a) and the classification obtained from the fit of a GMM (b)."----
clp <- clPairs(X, Class, symbols = c(16,0), colors = "black")
clPairsLegend(0.8, 0.2, class = clp$class, col = clp$col, pch = clp$pch)
plot(mod, what = "classification", fillEllipses = TRUE)


## ----mclust_hemophilia_boot, cache=TRUE-----------------------------------------
boot <- MclustBootstrap(mod, nboot = 999, type = "bs")


## -------------------------------------------------------------------------------
summary(boot, what = "se")


## -------------------------------------------------------------------------------
summary(boot, what = "ci", conf.level = 0.9)


## ----hemophilia_boot1, fig.width=10, fig.height=4, out.width="\\textwidth", fig.cap="Bootstrap distribution for the mixture proportions of the GMM fitted to the \\code{hemophilia}\\ data. The vertical dotted lines indicate the MLEs, the bottom line indicates the percentile confidence intervals, and the square shows the center of the bootstrap distribution."----
par(mfrow = c(1,2))
plot(boot, what = "pro")


## ----hemophilia_boot2, fig.width=10, fig.height=8, out.width="\\textwidth", fig.cap="Bootstrap distribution for the mixture component means of the GMM fitted to the \\code{hemophilia}\\ data. The vertical dotted lines indicate the MLEs, the bottom line indicates the percentile confidence intervals, and the square shows the center of the bootstrap distribution."----
par(mfcol = c(2,2))
plot(boot, what = "mean")


## ----mclust_hemophilia_wlboot, cache=TRUE---------------------------------------
wlboot <- MclustBootstrap(mod, nboot = 999, type = "wlbs")
summary(wlboot, what = "se")


## ----hemophilia_boot3, echo=-1, fig.width=6, fig.height=5, out.width="0.8\\textwidth", fig.cap="Bootstrap percentile intervals for the means of the GMM fitted to the \\code{hemophilia}\\ data. Solid lines refer to the nonparametric bootstrap, dashed lines to the weighted likelihood bootstrap."----
par(mfrow = c(1,2), mar = c(4,4,1,1))
boot.ci <- summary(boot, what = "ci")
wlboot.ci <- summary(wlboot, what = "ci")
for (j in 1:mod$G)
{ 
  plot(1:mod$G, mod$parameters$mean[j,], col = 1:mod$G, pch = 15,
       ylab = colnames(X)[j], xlab = "Mixture component",
       ylim = range(boot.ci$mean, wlboot.ci$mean), 
       xlim = c(.5,mod$G+.5), xaxt = "n")
  points(1:mod$G+0.2, mod$parameters$mean[j,], col = 1:mod$G, pch = 15)
  axis(side = 1, at = 1:mod$G)
  with(boot.ci, 
       errorBars(1:G, mean[1,j,], mean[2,j,], col = 1:G))
  with(wlboot.ci, 
       errorBars(1:G+0.2, mean[1,j,], mean[2,j,], col = 1:G, lty = 2))
}


## ----precip_dotchart, fig.width=5, fig.height=7, out.width="0.8\\textwidth", fig.cap="Dot chart of average annual rainfall (in inches) for 70 US cities."----
data(precip, package = "datasets")
dotchart(sort(precip), cex = 0.6, pch = 19,
         xlab = "Average annual rainfall (in inches)")


## ----precip_bic, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="BIC traces for GMMs fitted to the \\code{precip}\\ data."----
mod <- Mclust(precip)
summary(mod, parameters = TRUE)
plot(mod, what = "BIC", legendArgs = list(x = "bottomleft"))


## ----precip_classif_uncert, echo=-1, fig.width=10, fig.height=4, out.width="\\textwidth", fig.subcap="", fig.cap="Classification (left) and uncertainty (right) plots for the \\code{precip}\\ data."----
par(mfrow=c(1,2))
plot(mod, what = "classification")
plot(mod, what = "uncertainty")


## ----precip_dotchart_cluster, fig.width=5, fig.height=7, out.width="0.8\\textwidth", fig.cap="Dot chart of average annual rainfall (in inches) for 70 US cities grouped by the estimated clustering partitions."----
x <- data.frame(precip, clusters = mod$classification)
rownames(x) <- make.unique(names(precip)) # correct duplicated names 
x <- x[order(x$precip),]
dotchart(x$precip, labels = rownames(x),
         groups = factor(x$clusters, labels = c("Cluster 1", "Cluster 2")),
         cex = 0.6, pch = 19, 
         color = mclust.options("classPlotColors")[x$clusters],
         xlab = "Average annual rainfall (in inches)")


## -------------------------------------------------------------------------------
data(EuroUnemployment, package = "mclust")
summary(EuroUnemployment)


## -------------------------------------------------------------------------------
HC_EII <- hc(EuroUnemployment, modelName = "EII", use = "VARS")
HC_VVV <- hc(EuroUnemployment, modelName = "VVV", use = "VARS")


## ----HCmerge, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap="", fig.cap="Dendrograms with height corresponding to the number of groups for model-based hierarchical clustering of the \\code{EuroUnemployment}\\ data with the \\code{EII}\\ (a) and \\code{VVV}\\ (b) models."----
plot(HC_EII, what="merge", labels = TRUE, hang = 0.02)
plot(HC_VVV, what="merge", labels = TRUE, hang = 0.02)


## ----HCloglik, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap="", fig.cap="Dendrograms with height corresponding to the classification log-likelihood for model-based hierarchical clustering of the \\code{EuroUNemployment}\\ data with the \\code{EII}\\ (a) and \\code{VVV}\\ (b) models. The log-likelihood is undefined at the other levels of the trees due to the presence of singletons."----
plot(HC_EII, what="loglik")
plot(HC_VVV, what="loglik")


## -------------------------------------------------------------------------------
data(HRstars, package = "GDAdata")
set.seed(0)
initial <- kmeans(HRstars[,-1], centers = 100)$cluster
HC_VVV <- hc(HRstars[,-1], modelName = "VVV", 
             partition = initial, use = "VARS")


## -------------------------------------------------------------------------------
data(flea, package = "tourr")
X <- data.matrix(flea[,1:6])
Class <- factor(flea$species, 
                labels = c("Concinna","Heikertingeri","Heptapotamica")) 
table(Class)


## ----flea_clpairs, fig.width=7, fig.height=7, out.width="\\textwidth", fig.cap="Scatterplot matrix for the \\code{flea}\\ dataset with points marked according to the true classes."----
col <- mclust.options("classPlotColors")
clp <- clPairs(X, Class, lower.panel = NULL, gap = 0,
               symbols = c(16,15,17), 
               colors = adjustcolor(col, alpha.f = 0.5))
clPairsLegend(x = 0.1, y = 0.3, class = clp$class, 
              col = col, pch = clp$pch,
              title = "Flea beatle species")


## -------------------------------------------------------------------------------
# set the default for the current session
mclust.options("hcUse" = "VARS")
mod1 <- Mclust(X)
# or specify the initialization method only for this model
# mod1 <- Mclust(X, initialization = list(hcPairs = hc(X, use = "VARS")))
summary(mod1)
table(Class, mod1$classification)
adjustedRandIndex(Class, mod1$classification)


## -------------------------------------------------------------------------------
mod2 <- Mclust(X[,6:1])
summary(mod2)
table(Class, mod2$classification)
adjustedRandIndex(Class, mod2$classification)


## -------------------------------------------------------------------------------
mod3 <- Mclust(X, initialization = list(hcPairs = hcRandomPairs(X, seed = 1)))
summary(mod3)
table(Class, mod3$classification)
adjustedRandIndex(Class, mod3$classification)


## ----flea_rndstarts, echo=-1, cache=TRUE----------------------------------------
set.seed(20190603)
BIC <- NULL
for (i in 1:50)
{
  # get BIC table from initial random start
  BIC0 <- mclustBIC(X, verbose = FALSE,
                    initialization = list(hcPairs = hcRandomPairs(X)))
  # update BIC table by merging best BIC values for each
  # G and modelNames
  BIC  <- mclustBICupdate(BIC, BIC0)
}
summary(BIC, k = 5)


## ----flea_rndstarts_model, cache=TRUE-------------------------------------------
mod4 <- Mclust(X, x = BIC)
summary(mod4)
table(Class, mod4$classification)
adjustedRandIndex(Class, mod4$classification)


## -------------------------------------------------------------------------------
mclust.options("hcUse" = "SVD")  # restore the default
mod5 <- Mclust(X)
# or specify only for this model fit
# mod5 <- Mclust(X, initialization = list(hcPairs = hc(X, use = "SVD")))
summary(mod5)
table(Class, mod5$classification)
adjustedRandIndex(Class, mod5$classification)  


## -------------------------------------------------------------------------------
data(iris, package = "datasets")
str(iris)
ms <- mstep(iris[,1:4], modelName = "VVV",
            z = unmap(iris$Species))
str(ms,1)
es <- estep(iris[,1:4], modelName = "VVV", 
            parameters = ms$parameters)
str(es,1)

