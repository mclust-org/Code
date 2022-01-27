## ----chevron1, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Simulated minefield data (a), marked to distinguish the known noise (b)."----
data(chevron, package = "mclust")
summary(chevron)
noise <- with(chevron, class == "noise")
X <- chevron[,2:3]
plot(X, cex = 0.5)
plot(X, cex = 0.5, col = ifelse(noise, "grey", "black"), 
     pch = ifelse(noise, 3, 1))


## ----fig.keep="none"------------------------------------------------------------
library(prabclus)
nnc <- NNclean(X, k = 5)
table(nnc$z)
clPairs(X, nnc$z, colors = c("darkgrey", "black"), symbols = c(3, 1))


## ----fig.keep="none"------------------------------------------------------------
modNoise <- Mclust(X, initialization = list(noise = (nnc$z == 0)))
summary(modNoise$BIC)


## ----fig.keep="none"------------------------------------------------------------
summary(modNoise, parameters = TRUE)
addmargins(table(chevron$class, modNoise$classification), 2)
plot(modNoise, what = "classification")


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
clPairs(X, nnc$z, colors = c("darkgrey", "black"), symbols = c(3, 1))


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(modNoise, what = "classification")


## -------------------------------------------------------------------------------
hypvol(X)


## -------------------------------------------------------------------------------
library(cluster)
ehull <- ellipsoidhull(as.matrix(X))
volume(ehull)
modNoise.ehull <- Mclust(X, Vinv = 1/volume(ehull),
                         initialization = list(noise = (nnc$z == 0)))
summary(modNoise.ehull)
tab <- table(chevron$class, modNoise.ehull$classification)
addmargins(tab, 2)


## -------------------------------------------------------------------------------
library(geometry)
chull <- convhulln(X, options = "FA")
chull$vol


## ----cache=TRUE-----------------------------------------------------------------
modNoise.chull <- Mclust(X, Vinv = 1/chull$vol,
                         initialization = list(noise = (nnc$z == 0)))
summary(modNoise.chull)
tab <- table(chevron$class, modNoise.chull$classification)
addmargins(tab, 2)


## -------------------------------------------------------------------------------
library(covRobust)
nnve <- cov.nnve(X, k = 5)
table(nnve$classification)


## -------------------------------------------------------------------------------
modNoise.nnve <- Mclust(X, initialization = 
                             list(noise = (nnve$classification == 0)))
summary(modNoise.nnve$BIC)
summary(modNoise.nnve)
addmargins(table(chevron$class, modNoise.nnve$classification), 2)


## ----echo=-1, fig.keep="none"---------------------------------------------------
set.seed(123)
(Sigma <- array(c(2,0.5,0.5,0.5, 1,0,0,0.1, 2,-0.5,-0.5,0.5), 
                dim = c(2,2,3)))
var.decomp <- sigma2decomp(Sigma)
str(var.decomp)
par <- list(pro = c(1/3, 1/3, 1/3), 
            mean = cbind(c(0,3), c(3,0), c(-3,0)),
            variance = var.decomp)
data <- sim(par$variance$modelName, parameters = par, n = 200)
noise <- matrix(runif(100, -10, 10), nrow = 50, ncol = 2)
X <- rbind(data[,2:3], noise)
cluster <- c(data[,1], rep(0, 50))
clPairs(X, ifelse(cluster == 0, 0, 1), 
        colors = "black", symbols = c(16,1), cex = c(0.5,1))


## -------------------------------------------------------------------------------
nnc <- NNclean(X, k = 5)
modNoise <- Mclust(X, 
  initialization = list(noise = (nnc$z == 0)))
summary(modNoise$BIC)
summary(modNoise, parameters = TRUE)


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
clPairs(X, ifelse(cluster == 0, 0, 1), 
        colors = "black", symbols = c(16,1), cex = c(0.5,1))


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(modNoise, what = "classification")


## ----fig.keep="none"------------------------------------------------------------
plot(modNoise, what = "classification")
table(cluster, Classification = modNoise$classification)


## ----galaxies-------------------------------------------------------------------
data(galaxies, package = "MASS") 
# this fix a typographic error in the data, see help(galaxies, package = "MASS")
galaxies[78] <- 26960 
galaxies <- galaxies / 1000


## ----echo=-1, fig.keep="none"---------------------------------------------------
set.seed(20190410)
BIC <- NULL
for(i in 1:50)
{
  BIC0 <- mclustBIC(galaxies, verbose = FALSE,
                    initialization = list(hcPairs = hcRandomPairs(galaxies)))
  BIC  <- mclustBICupdate(BIC, BIC0)
}
summary(BIC, k = 5)
plot(BIC)

mod <- densityMclust(galaxies, x = BIC)
summary(mod, parameters = TRUE)

plot(mod, what = "density", data = galaxies, breaks = 11)
rug(galaxies)


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(BIC)


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(mod, what = "density", data = galaxies, breaks = 11)
rug(galaxies)


## ----echo=-1, fig.keep="none"---------------------------------------------------
set.seed(20190411)
BICp <- NULL
for(i in 1:50)
{
  BIC0p <- mclustBIC(galaxies, verbose = FALSE, 
                     prior = priorControl(),
            initialization = list(hcPairs = hcRandomPairs(galaxies)))
  BICp  <- mclustBICupdate(BICp, BIC0p)
}
summary(BICp, k = 5)
plot(BICp)

modp <- densityMclust(galaxies, x = BICp)
summary(modp, parameters = TRUE)

plot(modp, what = "density", data = galaxies, breaks = 11)
rug(galaxies)


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(BICp)


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(modp, what = "density", data = galaxies, breaks = 11)
rug(galaxies)


## -------------------------------------------------------------------------------
defaultPrior(galaxies, G = 3, modelName = "V")


## ----olive oil------------------------------------------------------------------
data(olive, package = "pgmm")
# recode of labels for Region and Area
olive$Region <- factor(olive$Region, levels = 1:3, 
                       labels = c("South", "Sardinia", "North"))
olive$Area <- factor(olive$Area, levels = 1:9, 
                     labels = c("North Apulia", "Calabria", "South Apulia", 
                                "Sicily", "Inland Sardinia", "Costal Sardinia", 
                                "East Liguria", "West Liguria", "Umbria"))
with(olive, table(Area, Region))


## ----oliveoil01, fig.keep="none"------------------------------------------------
X <- scale(olive[,3:10])

BIC <- mclustBIC(X, G = 1:15)
summary(BIC)
plot(BIC, legendArgs = list(x = "bottomright", ncol = 5))

BICp <- mclustBIC(X, G = 1:15, prior = priorControl())
summary(BICp)
plot(BICp, legendArgs = list(x = "bottomright", ncol = 5))


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(BIC, legendArgs = list(x = "bottomright", ncol = 5))


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(BICp, legendArgs = list(x = "bottomright", ncol = 5))


## -------------------------------------------------------------------------------
mod <- Mclust(X, x = BIC)
summary(mod)

table(Region = olive$Region, cluster = mod$classification)
adjustedRandIndex(olive$Region, mod$classification)

table(Area = olive$Area, cluster = mod$classification)
adjustedRandIndex(olive$Area, mod$classification)


## -------------------------------------------------------------------------------
modp <- Mclust(X, x = BICp)
summary(modp)

table(Region = olive$Region, cluster = modp$classification)
adjustedRandIndex(olive$Region, modp$classification)

table(Area = olive$Area, cluster = modp$classification)
adjustedRandIndex(olive$Area, modp$classification)


## ----ex4-1_fig1, fig.width=6, fig.height=5, out.width="0.6\\textwidth", fig.cap="The \\code{ex4.1}\\ simulated dataset."----
data(Baudry_etal_2010_JCGS_examples, package = "mclust")
plot(ex4.1)


## ----cache=TRUE-----------------------------------------------------------------
mod_ex4.1 <- Mclust(ex4.1)
summary(mod_ex4.1)


## -------------------------------------------------------------------------------
CLUSTCOMBI <- clustCombi(mod_ex4.1)
summary(CLUSTCOMBI)


## ----ex4-1_clustcombi, fig.width=6, fig.height=9, out.width="\\textwidth", fig.cap="Hierarchy of mixture component combinations for the \\code{ex4.1}\\ dataset."----
par(mfrow = c(3,2), mar = c(4,4,3,1))
plot(CLUSTCOMBI, ex4.1, what = "classification")


## ----ex4-1_clustcombitree, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap=rep("",2), fig.cap="Dendrogram tree plots of the merging structure of mixture component combinations for the simulated \\code{ex4.1}\\ dataset."----
plot(CLUSTCOMBI, what = "tree")
plot(CLUSTCOMBI, what = "tree", type = "rectangle", yaxis = "step")


## ----ex4-1_clustcombioptim, fig.width=6, fig.height=5, out.width="0.6\\textwidth", fig.cap="Entropy plot for selecting the final clustering by merging mixture components for the simulated \\code{ex4.1}\\ dataset."----
optimClust <- clustCombiOptim(CLUSTCOMBI, reg = 2, plot = TRUE)
str(optimClust)


## ----fig.keep="none"------------------------------------------------------------
data(faithful, package = "datasets")
mod <- Mclust(faithful)
clPairs(faithful, mod$classification)
plot(as.densityMclust(mod), what = "density", add = TRUE)


## ----fig.keep="none", cache=TRUE------------------------------------------------
GMMHD <- gmmhd(mod)
summary(GMMHD)


## ----fig.keep="none"------------------------------------------------------------
plot(GMMHD, what = "mode")
plot(GMMHD, what = "cores")
plot(GMMHD, what = "clusters")


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
clPairs(faithful, mod$classification)
plot(as.densityMclust(mod), what = "density", add = TRUE)


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(GMMHD, what = "mode")


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(GMMHD, what = "cores")


## ----echo=FALSE, out.width="0.48\\textwidth"------------------------------------
plot(GMMHD, what = "clusters")


## ----cache=TRUE-----------------------------------------------------------------
data(yeast, package="MixSAL")
dup <- duplicated(yeast)
sum(dup) # count replicated observations
mod <- Mclust(yeast[!dup,-1])
summary(mod)
adjustedRandIndex(yeast$Site[!dup], mod$classification)


## ----cache=TRUE-----------------------------------------------------------------
GMMHD <- gmmhd(mod)
summary(GMMHD)


## -------------------------------------------------------------------------------
table(yeast$Site[!dup], cluster = GMMHD$cluster)
adjustedRandIndex(yeast$Site[!dup], GMMHD$cluster)


## ----gmmhd_yeast, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Plot of yeast data projected along the first two GMMDR directions. Points are marked by cluster cores in panel (a) and by final clustering in panel (b): \\textcolor{dodgerblue2}{$\\bullet$} = ME3, \\textcolor{red3}{$\\blacktriangle$} = CYT, \\textcolor{gray50}{\\footnotesize$\\square$} = unlabeled."----
plot(GMMHD, what = "cores", 
     col = c("grey50", "red3", "dodgerblue2"), 
     pch = c(0, 17, 16))
surfacePlot(GMMHD$x, parameters = GMMHD$MclustDR$parameters,
            what = "density", type = "contour", nlevels = 15,
            col = "black", drawlabels = FALSE, add = TRUE)
plot(GMMHD, what = "clusters", 
     col = c("red3", "dodgerblue2"), pch = c(17, 16))


## ----cache=TRUE-----------------------------------------------------------------
data(Baudry_etal_2010_JCGS_examples, package = "mclust")
mod <- Mclust(ex4.1)
GMMHD <- gmmhd(mod)
summary(GMMHD)
optimClust <- clustCombiOptim(clustCombi(mod), reg = 2)
table(GMMHD = GMMHD$cluster, CLUSTCOMBI = optimClust$cluster)


## -------------------------------------------------------------------------------
mod <- Mclust(faithful)
sim0 <- sim(modelName = mod$modelName, 
            parameters = mod$parameters,
            n = nrow(faithful), seed = 0)
sim1 <- sim(modelName = mod$modelName, 
            parameters = mod$parameters,
            n = nrow(faithful), seed = 1)


## ----fig.keep="none"------------------------------------------------------------
xlim <- range(c(faithful[,1],sim0[,2],sim1[,2]))
ylim <- range(c(faithful[,2],sim0[,3],sim1[,3]))

mclust2Dplot(data = sim0[,-1], parameters = mod$parameters,
             classification = sim0[,1], xlim = xlim, ylim = ylim)

mclust2Dplot(data = sim1[,-1], parameters = mod$parameters,
             classification = sim1[,1], xlim = xlim, ylim = ylim)


## -------------------------------------------------------------------------------
head(sim0)


## -------------------------------------------------------------------------------
par <- list(
  pro = c(0.5, 0.3, 0.2),
  mean = matrix(c(0,0,3,3,-4,1), nrow = 2, ncol = 3),
  variance = sigma2decomp(matrix(c(1,0.6,0.6,1.5), nrow = 2, ncol = 2), 
                          G = 3))
str(par)


## -------------------------------------------------------------------------------
sim <- sim(modelName = "EEI", parameters = par, n = 10000, seed = 123)
cluster <- sim[,1]
x <- sim[,2:3]


## ----eval=FALSE-----------------------------------------------------------------
## mclust.options(subset = Inf)  # no subsetting
## system.time(mod1 <- Mclust(x))
## ##    user  system elapsed
## ## 172.981   0.936 173.866
## summary(mod1$BIC)
## ## Best BIC values:
## ##              EEI,3         EEE,3        EVI,3
## ## BIC      -76523.89 -76531.836803 -76541.59835
## ## BIC diff      0.00     -7.948006    -17.70955
## 
## mclust.options(subset = 2000)  # reset to default setting
## system.time(mod2 <- Mclust(x))
## ##  user  system elapsed
## ## 14.247   0.196  14.395
## summary(mod2$BIC)
## ## Best BIC values:
## ##              EEI,3         EEE,3        EVI,3
## ## BIC      -76524.09 -76531.845624 -76541.89691
## ## BIC diff      0.00     -7.754063    -17.80535
## 
## s <- sample(1:nrow(x), size = 1000)  # one-time subsetting
## system.time(mod3 <- Mclust(x, initialization = list(subset = s)))
## ##   user  system elapsed
## ## 12.091   0.146  12.197
## summary(mod3$BIC)
## ## Best BIC values:
## ##              EEI,3         EEE,3        VEI,3
## ## BIC      -76524.17 -76532.043460 -76541.78983
## ## BIC diff      0.00     -7.874412    -17.62078


## ----eval=FALSE-----------------------------------------------------------------
## table(mod1$classification, mod2$classification)
## ##       1    2    3
## ##  1 5090    0    0
## ##  2    0    0 3007
## ##  3   10 1893    0


## -------------------------------------------------------------------------------
data(stlouis, package = "mix")
x <- data.frame(stlouis[,-(1:3)], row.names = NULL)
table(complete.cases(x))
apply(x, 2, function(x) prop.table(table(complete.cases(x))))


## ----stlouis_missing, fig.width=6, fig.height=6, out.width="0.7\\textwidth", fig.cap="Image plot of the missing values pattern for the \\code{stlouis}\\ dataset."----
library(ggplot2)
df <- data.frame(obs = rep(1:nrow(x), times = ncol(x)),
                 var = rep(colnames(x), each = nrow(x)),
                 missing = as.vector(is.na(x)))
ggplot(data = df, aes(x = var, y = obs)) +
  geom_tile(aes(fill = missing)) +
  scale_fill_manual(values = c("lightgrey", "black")) +
  labs(x = "Variables", y = "Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(margin = margin(b = 10))) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank())


## -------------------------------------------------------------------------------
ximp <- imputeData(x, seed = 123)


## ----stlouis_imputePairs, fig.width=7, fig.height=6, out.width="0.9\\textwidth", fig.cap="Pairs plot showing an instance of applying \\code{imputeData()}\\ to the continuous variables in the \\code{stlouis}\\ dataset. The open black circles correspond to non-missing data, while the filled red circles correspond to imputed missing values."----
imputePairs(x, ximp)


## -------------------------------------------------------------------------------
x <- as.matrix(iris[,1:4])
mod <- Mclust(x, G = 3, modelNames = "EEE")
table(iris$Species, mod$classification)
adjustedRandIndex(iris$Species, mod$classification)
mod$parameters$mean  # component means


## ----echo=-1--------------------------------------------------------------------
set.seed(20171211)
M <- sample(c(TRUE,FALSE), size = prod(dim(x)), 
            replace = TRUE, prob = c(0.1,0.9))
x[M] <- NA
table(cmpObs <- complete.cases(x))


## ----echo=-1--------------------------------------------------------------------
set.seed(20171212)
nImp  <- 100
muImp <- array(NA, c(ncol(x), 3, nImp))
clImp <- array(NA, c(nrow(x), nImp))
for(i in 1:nImp)
{ 
  xImp <- imputeData(x, verbose = FALSE)
  modImp <- Mclust(xImp, G = 3, modelNames = "EEE", verbose = FALSE)
  if (i == 1) clImp[,i] <- modImp$classification
  mcl <- matchCluster(clImp[,1], modImp$classification)
  clImp[,i]  <- mcl$cluster
  muImp[,,i] <- modImp$parameters$mean[,mcl$ord]
}

# majority rule
cl <- apply(clImp, 1, function(x) majorityVote(x)$majority)
table(iris$Species, cl)
adjustedRandIndex(iris$Species, cl)

# pooled estimate of cluster means
apply(muImp, 1:2, mean)

