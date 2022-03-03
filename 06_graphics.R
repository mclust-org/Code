## ---------------------------------------------------------------------------
data("Snapper", package = "FSAdata")
x <- Snapper[,1]
mod <- Mclust(x, G = 4, modelNames = "V")
summary(mod, parameters = TRUE)


## ----fishery2, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Classification (a) and uncertainty (b) plots created with \\code{mclust1Dplot()}\\ for the GMM fit to the fishery data."----
mclust1Dplot(x, what = "classification", 
             parameters = mod$parameters, z = mod$z,
             xlab = "Fish length")

mclust1Dplot(x, what = "uncertainty",
             parameters = mod$parameters, z = mod$z,
             xlab = "Fish length")


## ----faithful2Dplots1, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Classification (a) and uncertainty (b) plots created with \\code{mclust2Dplot()}\\ for the model fitted with \\code{Mclust()}\\ to the \\code{faithful}\\ dataset."----
mod <- Mclust(faithful)

mclust2Dplot(data = faithful, what = "classification", 
             parameters = mod$parameters, z = mod$z)

mclust2Dplot(data = faithful, what = "uncertainty",
             parameters = mod$parameters, z = mod$z)


## ----eval=FALSE-------------------------------------------------------------
## surfacePlot(data = faithful, parameters = mod$parameters,
##             what = "density", type = "contour",
##             transformation = "log")
## 
## surfacePlot(data = faithful, parameters = mod$parameters,
##             what = "density", type = "image")
## 
## surfacePlot(data = faithful, parameters = mod$parameters,
##             what = "density", type = "persp")
## 
## surfacePlot(data = faithful, parameters = mod$parameters,
##             what = "uncertainty", type = "image",
##             transformation = "sqrt")


## ----eval=FALSE, echo=FALSE-------------------------------------------------
## # Example fo tuning the col.palette arg
## surfacePlot(data = faithful, parameters = mod$parameters,
##              what = "density", type = "image",
##              col.palette = function(...) hcl.colors(..., "Geyser"))


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "density", type = "contour",
            transformation = "log")


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "density", type = "image")


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "density", type = "persp")


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "uncertainty", type = "image",
            transformation = "sqrt")


## ---------------------------------------------------------------------------
data("iris", package = "datasets")
mod <- Mclust(iris[,1:4], G = 3)
summary(mod)


## ----fig.keep="none"--------------------------------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "classification",
          parameters = mod$parameters, z = mod$z)

coordProj(data = iris[,1:4], dimens = c(2,4), what = "uncertainty",
          parameters = mod$parameters, z = mod$z)

coordProj(data = iris[,1:4], dimens = c(2,4), what = "error",
          parameters = mod$parameters, z = mod$z, truth = iris$Species)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "classification",
          parameters = mod$parameters, z = mod$z)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "uncertainty",
          parameters = mod$parameters, z = mod$z)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "error",
          parameters = mod$parameters, z = mod$z, truth = iris$Species)


## ----echo=-1----------------------------------------------------------------
set.seed(1)
nrow(iris)
Q <- randomOrthogonalMatrix(nrow(iris), 2)
dim(Q)
QTQ <- crossprod(Q)                # equivalently t(Q) %*% Q
zapsmall(QTQ)                      # 2 x 2 identity matrix


## ----fig.keep="none"--------------------------------------------------------
randProj(data = iris[,1:4], seeds = c(1,13,79,201), what = "classification",
         parameters = mod$parameters, z = mod$z)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
randProj(data = iris[,1:4], seed = 1, what = "classification",
         parameters = mod$parameters, z = mod$z)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
randProj(data = iris[,1:4], seed = 13, what = "classification",
         parameters = mod$parameters, z = mod$z)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
randProj(data = iris[,1:4], seed = 79, what = "classification",
         parameters = mod$parameters, z = mod$z)


## ----echo=FALSE, out.width="0.48\\textwidth"--------------------------------
randProj(data = iris[,1:4], seed = 201, what = "classification",
         parameters = mod$parameters, z = mod$z)


## ----iris_crimcoords, fig.width=6, fig.height=5, out.width="0.6\\textwidth", fig.cap="Discriminant coordinates or crimcoords projection for the clustering for the 3-group model-based clustering of the \\code{iris}\\ dataset."----
crimcoords(iris[,1:4], mod$classification)


## ---------------------------------------------------------------------------
data("thyroid", package = "mclust")
CRIMCOORDS <- crimcoords(thyroid[,-1], thyroid$Diagnosis, plot = FALSE)
zapsmall(eval <- CRIMCOORDS$evalues)
eval[1:2]/sum(eval)
CRIMCOORDS$basis


## ----thyroid_crimcoords, fig.width=6, fig.height=5, out.width="0.6\\textwidth", fig.cap="Discriminant coordinates or crimcoords projection for the diagnosis classes of the \\code{thyroid}\\ dataset. Group centroids are represented with the $+$ symbol."----
plot(CRIMCOORDS)
points(CRIMCOORDS$means %*% CRIMCOORDS$basis, pch = 3, cex = 1.5)
legend("topright", legend = levels(thyroid$Diagnosis), inset = 0.02,
       col = mclust.options("classPlotColors")[1:3],
       pch = mclust.options("classPlotSymbols")[1:3])


## ---------------------------------------------------------------------------
data("wine", package = "gclus")
Class <- factor(wine$Class, levels = 1:3,
                labels = c("Barolo", "Grignolino", "Barbera"))
X <- data.matrix(wine[,-1])
mod <- Mclust(X, G = 3, modelNames = "VVE")
table(Class, Cluster = mod$classification)


## ----graphs_wine_dr---------------------------------------------------------
drmod <- MclustDR(mod, lambda = 1)
summary(drmod)


## ----eval=FALSE-------------------------------------------------------------
## plot(drmod, what = "contour")


## ----eval=FALSE-------------------------------------------------------------
## plot(drmod, what = "boundaries")


## ----eval=FALSE-------------------------------------------------------------
## miscl <- classError(Class, mod$classification)$misclassified
## points(drmod$dir[miscl,], pch = 1, cex = 2)


## ----echo=FALSE, out.width="0.49\\textwidth"--------------------------------
plot(drmod, what = "contour")


## ----echo=FALSE, out.width="0.49\\textwidth"--------------------------------
plot(drmod, what = "boundaries")


## ---------------------------------------------------------------------------
mod <- Mclust(iris[, 1:4], G = 3)
drmod <- MclustDR(mod, lambda = .5)
summary(drmod)


## ----irisdr3, fig.width=6, fig.height=5, out.width="0.6\\textwidth", fig.cap="Plot of generalized eigenvalues from \\code{MclustDR()}\\ and the corresponding contributions from means and variances for the 3-group model-based clustering of the \\code{iris}\\ dataset."----
plot(drmod, what = "evalues")


## ----irisdr4, fig.width=8, fig.height=7, out.width="\\textwidth", fig.cap="The \\code{iris}\\ data projected onto the principal eigenvectors from \\code{MclustDR()}. The colors and symbols correspond to the 3-group model-based clustering."----
plot(drmod, what = "pairs")


## ----graphs_iris_drsubsel, cache=TRUE---------------------------------------
sdrmod <- MclustDRsubsel(drmod, verbose = TRUE)
summary(sdrmod)


## ---------------------------------------------------------------------------
zapsmall(cor(drmod$dir, sdrmod$dir)^2)


## ----eval=FALSE-------------------------------------------------------------
## plot(sdrmod, what = "contour", nlevels = 7)
## plot(sdrmod, what = "classification")
## plot(sdrmod, what = "boundaries")


## ----echo=FALSE, out.width="0.49\\textwidth"--------------------------------
plot(sdrmod, what = "contour", nlevels = 7)


## ----echo=FALSE, out.width="0.49\\textwidth"--------------------------------
plot(sdrmod, what = "boundaries")


## ----echo=FALSE, out.width="0.49\\textwidth"--------------------------------
plot(sdrmod, what = "classification")


## ----echo=FALSE, out.width="0.49\\textwidth"--------------------------------
plot(sdrmod, what = "density")


## ----eval=FALSE-------------------------------------------------------------
## plot(sdrmod, what = "density")


## ---------------------------------------------------------------------------
data("banknote", package = "mclust")
mod <- MclustDA(data = banknote[, -1], class = banknote$Status)
summary(mod)


## ---------------------------------------------------------------------------
drmod <- MclustDR(mod, lambda = .5)
summary(drmod)


## ----banknote1, fig.width=6, fig.height=5, out.width="0.6\\textwidth", fig.cap="Plot of generalized eigenvalues from \\code{MclustDR()}\\ applied to a classification model for the \\code{banknote}\\ data."----
plot(drmod, what = "evalues")


## ----banknote2, fig.width=8, fig.height=7, out.width="\\textwidth", fig.cap="Pairs plot of points projected onto the directions estimated with \\code{MclustDR()}\\ for the \\code{banknote}\\ dataset."----
plot(drmod, what = "pairs", lower.panel = NULL)
clPairsLegend(0.1, 0.4, class = levels(drmod$classification), 
              col = mclust.options("classPlotColors")[1:2],
              pch = mclust.options("classPlotSymbols")[1:2],
              title = "Swiss banknote data")


## ----graphs_banknote_drsubsel, cache=TRUE-----------------------------------
sdrmod <- MclustDRsubsel(drmod, verbose = TRUE)
summary(sdrmod)


## ---------------------------------------------------------------------------
zapsmall(cor(drmod$dir, sdrmod$dir)^2)


## ----banknote3, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap=rep("",2), fig.cap="Contour plot of mixture densities for each class (a) and MAP classification regions (b) drawn on the projection subspace estimated with \\code{MclustDRsubsel()}\\ for the \\code{banknote}\\ dataset."----
plot(sdrmod, what = "contour", nlevels = 15)
plot(sdrmod, what = "classification")


## ----ggplot_faithful, fig.width=6, fig.height=5, out.width="0.7\\textwidth", fig.cap="Scatterplot of the \\code{faithful}\\ data with points marked according to the GMM clusters identified by \\code{Mclust}."----
mod <- Mclust(faithful)
DF <- data.frame(mod$data, cluster = factor(mod$classification))
library("ggplot2")
ggplot(DF, aes(x = eruptions, y = waiting, 
               colour = cluster, shape = cluster)) +
  geom_point()


## ----ggplot_faithful_bic, fig.cap="BIC traces for the GMMs estimated for the \\code{faithful}\\ data.", fig.width=7, fig.height=5, out.width="\\textwidth"----
library("tidyr")
DF <- data.frame(mod$BIC[], G = 1:nrow(mod$BIC))
DF <- pivot_longer(DF, cols = 1:14, names_to = "Model", values_to = "BIC")
DF$Model <- factor(DF$Model, levels = mclust.options("emModelNames"))

ggplot(DF, aes(x = G, y = BIC, colour = Model, shape = Model)) +
  geom_point() + 
  geom_line() +
  scale_shape_manual(values = mclust.options("bicPlotSymbols")) +
  scale_color_manual(values = mclust.options("bicPlotColors")) +
  scale_x_continuous(breaks = unique(DF$G)) +
  xlab("Number of mixture components") +
  guides(shape = guide_legend(ncol=2))


## ---------------------------------------------------------------------------
mod <- Mclust(iris[, 1:4], G = 3)


## ----ggplot_iris_latent_profile, fig.cap="Latent profiles plot for the (\\code{VEV},3) model estimated for the \\code{iris}\\ data.", fig.width=6, fig.height=5, out.width="0.9\\textwidth"----
means <- data.frame(Profile = 1:mod$G, t(mod$parameters$mean))
means <- pivot_longer(means, cols = -1, 
                      names_to = "Variable",
                      values_to = "Mean")
means$Profile  <- factor(means$Profile)
means$Variable <- factor(means$Variable, 
                         levels = rownames(mod$parameters$mean))
means

ggplot(means, aes(Variable, Mean, group = Profile, 
                  shape = Profile, color = Profile)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = NULL, y = "Latent profiles means") +
  scale_color_manual(values = mclust.options("classPlotColors")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom")


## ----ggplot_iris_chull, fig.width=7, fig.height=5, out.width="0.9\\textwidth", fig.cap="Scatterplot of the first two GMMDR directions with added convex hulls for the \\code{iris}\\ data classes."----
damod <- MclustDA(iris[, 1:4], iris$Species)
drmod <- MclustDR(damod)
DF1 <- data.frame(drmod$dir[, 1:2], cl = damod$class)
DF2 <- do.call("rbind", by(DF1, DF1[, 3], 
                           function(x) x[chull(x), ]))
ggplot() + 
  geom_point(data = DF1, aes(x = Dir1, y = Dir2, color = cl, shape = cl)) + 
  geom_polygon(data = DF2, aes(x = Dir1, y = Dir2, fill = cl), alpha = 0.4) 


## ----ggplot_faithful_density, fig.cap="Plot of histogram and density estimated by \\code{densityMclust()}\\ for the waiting time of the \\code{faithful}\\ data.", fig.width=6, fig.height=5, out.width="0.7\\textwidth"----
mod <- densityMclust(faithful$waiting)
x <- extendrange(faithful$waiting, f = 0.1)
x <- seq(x[1], x[2], length.out = 101)
pred <- data.frame(x, density = predict(mod, newdata = x))
ggplot(faithful, aes(waiting)) +
  geom_histogram(aes(y = stat(density)), bins = 15, 
                 fill = "slategray3", colour = "grey92") +
  geom_line(data = pred, aes(x, density))


## ----echo=FALSE, eval=FALSE-------------------------------------------------
## # sarebbe carino fare un grafico del tipo in
## # https://ggplot2.tidyverse.org/reference/geom_density.html
## # https://stackoverflow.com/questions/40131443/kernel-density-estimation-in-ggplot-with-geom-density
## mod <- densityMclust(faithful$waiting)
## x <- extendrange(faithful$waiting, f = 0.1)
## x <- seq(x[1], x[2], length.out = 101)
## pred <- predict(as.Mclust(mod), newdata = x)
## str(pred)
## 
## matplot(x, pred$z)
## 
## DF <- data.frame(waiting = x, comp = pred$classification, z = pred$z)
## head(DF)
## DF <- tidyr::pivot_longer(DF, cols = 3:4, values_to = "z")
## DF
## 
## ggplot(DF, aes(x = waiting, y = stat(count), fill = factor(comp))) +
## #  geom_point()
##   geom_density(position = "fill")
## 
## A final example involves ``faceting'', an efficient technique for presenting information in panels conditioning on one or more variables.
## Consider the bootstrap procedure discussed in Section~\ref{sec:resamp}, where the \code{MclustBootstrap()} function is applied to the two-component \code{VVV} model for the hemophilia data from the \pkg{rrcov} package
## \citep{Rpkg:rrcov}. Using the information returned and stored in \code{boot}, a \pkg{ggplot2} plot of the bootstrap distribution for the mixing proportions can be obtained as follows:
## 

## ----ggplot_hemophilia_boot1, echo=-1, fig.cap="Bootstrap distribution for the mixture proportions of (\\code{VVV},2) model fitted to the \\code{hemophilia}\\ data.", fig.width=9, fig.height=5, out.width="0.9\\textwidth"----
load_cache("mclust_hemophilia_boot", path = "./Cache/")
DF <- data.frame(mixcomp = rep(1:boot$G, each = boot$nboot), 
                 pro = as.vector(boot$pro))
ggplot(DF, aes(x = pro)) + 
  geom_histogram(aes(y = stat(density)), bins = 15, 
                 fill = "slategray3", colour = "grey92") +
  facet_grid(~ mixcomp) +
  xlab("Mixing proportions") +
  ylab("Density of bootstrap distribution") 


## ----ggplot_hemophilia_boot2, echo=-(1:2), fig.cap="Bootstrap distribution for the mixture component means of (\\code{VVV},2) model fitted to the \\code{hemophilia}\\ data.", fig.width=8, fig.height=7, out.width="\\textwidth"----
DF0 <- data.frame(rbind(cbind(mixcomp = 1, boot$mean[, , 1]), 
                        cbind(mixcomp = 2, boot$mean[, , 2])))
DF0 <- tidyr::pivot_longer(DF0, cols = 2:3, names_to = "variable", values_to = "mean")
DF <- rbind(
  data.frame("mixcomp"  = 1,
             "variable" = rep(colnames(boot$mean[, , 1]), 
                              each = dim(boot$mean)[1]),
             "mean"     = as.vector(boot$mean[, , 1])),
  data.frame("mixcomp"  = 2,
             "variable" = rep(colnames(boot$mean[, , 2]), 
                              each = dim(boot$mean)[1]),
             "mean"     = as.vector(boot$mean[, , 2])))
ggplot(DF, aes(x = mean)) +
   geom_histogram(aes(y = stat(density)), bins = 15,
                  fill = "slategray3", colour = "grey92") +
   facet_grid(mixcomp ~ variable, scales = "free_x") +
   xlab("Means of mixture") +
   ylab("Density of bootstrap distribution")


## ---------------------------------------------------------------------------
mclust.options("bicPlotColors")
mclust.options("classPlotColors")


## ----echo=FALSE, eval=FALSE-------------------------------------------------
## # tentativo di mostrare un grafico con i colori, ma non mi convince
## bicPlotColors <- mclust.options("bicPlotColors")
## classPlotColors <- mclust.options("classPlotColors")
## 
## WongPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## 
## 
## palette(WongPalette)
## plot(1:8, pch = 15, cex = 2, col = 1:8)
## 
## palette.colors(8, palette = "Okabe-Ito")
## 
## plot(1:14, pch = 15, cex = 2, col = 1:8)
## 
## bicPlotColorsWong <- c(WongPalette, WongPalette[1:6])
## classPlotColorsWong <- WongPalette[-1]
## 
## par(mfrow=c(2, 2), mar = c(1, 1, 3, 0))
## n = 14
## image(1:n, 1, as.matrix(1:n), col = bicPlotColors[1:n],
##       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
## title("bicPlotColors", font = 1)
## image(1:n, 1, as.matrix(1:n), col = bicPlotColorsWong,
##       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
## title("bicPlotColorsWong", font = 1)
## n = 7
## image(1:n, 1, as.matrix(1:n), col = classPlotColors[1:n],
##       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
## title("classPlotColors", font = 1)
## image(1:n, 1, as.matrix(1:n), col = classPlotColorsWong,
##       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
## title("classPlotColorsWong", font = 1)


## ---------------------------------------------------------------------------
palette.colors(palette = "Okabe-Ito")


## ----echo=FALSE, eval=FALSE-------------------------------------------------
## # OLD
## # Wong color-blind-friendly palette
## WongPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
##                  "#0072B2", "#D55E00", "#CC79A7")
## 
## bicPlotColorsWong <- c(bicPlotColors[1:2], WongPalette[c(2:8, 2:6)])
## classPlotColorsWong <- WongPalette[-1]

## ---------------------------------------------------------------------------
# get and save default palettes
bicPlotColors <- mclust.options("bicPlotColors")
classPlotColors <- mclust.options("classPlotColors")

# set Okabe-Ito palette for use in mclust
bicPlotColors_Okabe_Ito <-
   palette.colors(palette = "Okabe-Ito")[c(9, 1, 2:8, 2:6, 9, 1)]
names(bicPlotColors_Okabe_Ito) <- names(bicPlotColors)

classPlotColorsWong <- palette.colors(palette = "Okabe-Ito")[-1]

mclust.options("bicPlotColors" = bicPlotColors_Okabe_Ito)
mclust.options("classPlotColors" = classPlotColorsWong)


## ----cbPalette, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = rep("",2), fig.cap="\\mclust plots with a color-blind-friendly palette."----
mod <- Mclust(iris[,3:4])
plot(mod, what = "BIC")
plot(mod, what = "classification")


## ---------------------------------------------------------------------------
mclust.options("bicPlotColors" = bicPlotColors)
mclust.options("classPlotColors" = classPlotColors)

