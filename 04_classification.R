## -------------------------------------------------------------------------------
data(wdbc, package = "mclust")
X <- wdbc[,c("Texture_mean", "Area_extreme", "Smoothness_extreme")]
Class <- wdbc[,"Diagnosis"]


## -------------------------------------------------------------------------------
set.seed(123)
train <- sample(1:nrow(X), size = round(nrow(X)*2/3), replace = FALSE)
X_train <- X[train,]
Class_train <- Class[train]
tab <- table(Class_train)
cbind(Counts = tab, "%" = prop.table(tab)*100)
X_test <- X[-train,]
Class_test <- Class[-train]
tab <- table(Class_test)
cbind(Counts = tab, "%" = prop.table(tab)*100)


## ----wdbc1, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Pairwise scatterplot matrix of selected features for the breast cancer data with points distinguished by tumor diagnosis."----
clp <- clPairs(X_train, Class_train, lower.panel = NULL)
clPairsLegend(0.1, 0.3, col = clp$col, pch = clp$pch, 
              class = ifelse(clp$class == "B", "Benign", "Malign"),
              title = "Breast cancer diagnosis:")


## -------------------------------------------------------------------------------
mod1 <- MclustDA(X_train, Class_train, modelType = "EDDA")
summary(mod1)


## -------------------------------------------------------------------------------
summary(mod1, parameters = TRUE)


## -------------------------------------------------------------------------------
summary(mod1, newdata = X_test, newclass = Class_test)


## ----wdbc2, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Pairwise scatterplot matrix of selected features for the breast cancer training data with points distinguished by observed classes and ellipses representing the Gaussian distribution estimated for each class by EDDA."----
plot(mod1, what = "scatterplot")


## ----wdbc3, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Pairwise scatterplot matrix of selected features for the breast cancer training data with points distinguished by observed classes. Filled black points represent those cases misclassified by the fitted EDDA model."----
plot(mod1, what = "error")


## -------------------------------------------------------------------------------
mod2 <- MclustDA(X_train, Class_train)
summary(mod2, newdata = X_test, newclass = Class_test)


## ----wdbc4, fig.width=7, fig.height=6, out.width="\\textwidth", fig.cap="Scatterplots of selected features for the breast cancer training data with points distinguished by observed classes and ellipses representing the Gaussian distribution estimated for each class by MclustDA."----
plot(mod2, what = "scatterplot")


## ----fig.keep='none'------------------------------------------------------------
plot(mod2, what = "scatterplot", dimens = c(1,2))


## -------------------------------------------------------------------------------
mod3 <- MclustDA(X_train, Class_train, G = 2, modelNames = "EEE")
summary(mod3, newdata = X_test, newclass = Class_test)


## ----cvmclustda_ce, echo=-1, cache=TRUE-----------------------------------------
set.seed(20190520)
cv1 <- cvMclustDA(mod1)
str(cv1)
cv2 <- cvMclustDA(mod2)
str(cv2)


## -------------------------------------------------------------------------------
unlist(cv1[c("ce", "se.ce", "brier", "se.brier")])
unlist(cv2[c("ce", "se.ce", "brier", "se.brier")])


## ----cvmclustda, echo=-1, cache=TRUE--------------------------------------------
set.seed(2)
models <- mclust.options("emModelNames")
tab_CE <- tab_Brier <- matrix(NA, nrow = length(models)+1, ncol = 5)
rownames(tab_CE) <- rownames(tab_Brier) <- 
  c(paste0("EDDA[", models, "]"), "MCLUSTDA")
colnames(tab_CE) <- colnames(tab_Brier) <- 
  c("Train", "10-fold CV", "se(CV)", "lower", "upper")
for (i in seq(models))
{
  mod <- MclustDA(X, Class, modelType = "EDDA", 
                  modelNames = models[i], verbose = FALSE)
  pred <- predict(mod, X)
  cv <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  #
  tab_CE[i,1] <- classError(pred$classification, Class)$errorRate
  tab_CE[i,2] <- cv$ce
  tab_CE[i,3] <- cv$se.ce
  tab_CE[i,4] <- cv$ce - cv$se.ce
  tab_CE[i,5] <- cv$ce + cv$se.ce
  #
  tab_Brier[i,1] <- BrierScore(pred$z, Class)
  tab_Brier[i,2] <- cv$brier
  tab_Brier[i,3] <- cv$se.brier
  tab_Brier[i,4] <- cv$brier - cv$se.brier
  tab_Brier[i,5] <- cv$brier + cv$se.brier
}
i <- length(models)+1
mod <- MclustDA(X, Class, modelType = "MclustDA", verbose = FALSE)
pred <- predict(mod, X)
cv <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
#
tab_CE[i,1] <- classError(pred$classification, Class)$errorRate
tab_CE[i,2] <- cv$ce
tab_CE[i,3] <- cv$se.ce
tab_CE[i,4] <- cv$ce - cv$se.ce
tab_CE[i,5] <- cv$ce + cv$se.ce
#
tab_Brier[i,1] <- BrierScore(pred$z, Class)
tab_Brier[i,2] <- cv$brier
tab_Brier[i,3] <- cv$se.brier
tab_Brier[i,4] <- cv$brier - cv$se.brier
tab_Brier[i,5] <- cv$brier + cv$se.brier


## -------------------------------------------------------------------------------
tab_CE


## ----cvmclustda_fig1, fig.width=6, fig.height=4, out.width="\\textwidth", fig.cap="Training and cross-validated misclassification error rates of Gaussian mixture classification models for the breast cancer data."----
library(ggplot2)
df <- data.frame(rownames(tab_CE), tab_CE)
colnames(df) <- c("model", "train", "cv", "se", "lower", "upper")
df$model <- factor(df$model, levels = rev(df$model))
ggplot(df, aes(x = model, y = cv, ymin = lower, ymax = upper)) +
  geom_point(aes(shape = "s1", color = "c1")) + 
  geom_errorbar(width = 0.5, col = "dodgerblue3") + 
  geom_point(aes(y = train, shape = "s2", color = "c2")) +
  scale_y_continuous(breaks = seq(0,0.2,by=0.01), lim = c(0,NA)) +
  scale_color_manual(name = "", 
                     breaks = c("c1", "c2"),
                     values = c("dodgerblue3", "black"),
                     labels = c("CV", "Train")) +
  scale_shape_manual(name = "", 
                     breaks = c("s1", "s2"),
                     values = c(19, 0),
                     labels = c("CV", "Train")) +
  ylab("Classification error") + xlab("") + coord_flip() +
  theme(legend.position = "top")


## -------------------------------------------------------------------------------
tab_Brier


## ----cvmclustda_fig2, fig.width=6, fig.height=4, out.width="\\textwidth", fig.cap="Training and cross-validated Brier scores of Gaussian mixture classification models for the breast cancer data."----
df <- data.frame(rownames(tab_Brier), tab_Brier)
colnames(df) <- c("model", "train", "cv", "se", "lower", "upper")
df$model <- factor(df$model, levels = rev(df$model))
ggplot(df, aes(x = model, y = cv, ymin = lower, ymax = upper)) +
  geom_point(aes(shape = "s1", color = "c1")) + 
  geom_errorbar(width = 0.5, col = "dodgerblue3") + 
  geom_point(aes(y = train, shape = "s2", color = "c2")) +
  scale_y_continuous(breaks = seq(0,0.2,by=0.01), lim = c(0,NA)) +
  scale_color_manual(name = "", 
                     breaks = c("c1", "c2"),
                     values = c("dodgerblue3", "black"),
                     labels = c("CV", "Train")) +
  scale_shape_manual(name = "", 
                     breaks = c("s1", "s2"),
                     values = c(19, 0),
                     labels = c("CV", "Train")) +
  ylab("Brier score") + xlab("") + coord_flip() +
  theme(legend.position = "top")


## -------------------------------------------------------------------------------
# confusion matrix
(tab <- table(Predict = cv1$classification, Class = Class_train))
tab[2,2]/sum(tab[,2])  # sensitivity
tab[1,1]/sum(tab[,1])  # specificity


## -------------------------------------------------------------------------------
threshold <- seq(0, 1, by = 0.01)
sensitivity <- specificity <- rep(NA, length(threshold))
for(i in 1:length(threshold))
{
  pred <- factor(ifelse(cv1$z[,"M"] > threshold[i], "M", "B"),
                 levels = c("B", "M"))
  tab <- table(pred, Class_train)
  sensitivity[i] <- tab[2,2]/sum(tab[,2])
  specificity[i] <- tab[1,1]/sum(tab[,1])
}


## ----eval=FALSE-----------------------------------------------------------------
## plot(1-specificity, sensitivity, type = "l", lwd = 2)  # ROC curve
## abline(h = c(0, 1), v = c(0, 1), lty = 3)  # limits of [0,1]x[0,1] region
## abline(a = 0, b = 1, lty = 2)  # line of random classification


## ----echo=-1--------------------------------------------------------------------
# ModelMetrics::auc(ifelse(Class_train=="M", 1, 0), cv1$z[,2])
auc_approx <- function(tpr, fpr)
{
  x <- 1-fpr
  y <- tpr
  dx <- c(diff(x), 0)
  dy <- c(diff(y), 0)
  sum(y * dx) + sum(dy * dx)/2
}

auc_approx(tpr = sensitivity, fpr = 1-specificity)


## ----echo=FALSE-----------------------------------------------------------------
threshold <- seq(0, 1, by = 0.01)
sensitivity2 <- specificity2 <- rep(NA, length(threshold))
for(i in 1:length(threshold))
{
  pred <- factor(ifelse(cv2$z[,"M"] > threshold[i], "M", "B"),
                 levels = c("B", "M"))
  tab <- table(pred, Class_train)
  sensitivity2[i] <- tab[2,2]/sum(tab[,2])
  specificity2[i] <- tab[1,1]/sum(tab[,1])
}


## ----wdbc_roc, echo=FALSE, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap="", fig.cap="ROC curves from the cross-validated predictions of selected (a) EDDA and (b) MclustDA models for the breast cancer data."----
plot(1-specificity, sensitivity, type = "l", lwd = 2)
abline(h = c(0, 1), v = c(0, 1), lty = 3)
abline(a = 0, b = 1, lty = 2)
plot(1-specificity2, sensitivity2, type = "l", lwd = 2,
     xlab = "1-specificity", ylab = "sensitivity") 
abline(h = c(0, 1), v = c(0, 1), lty = 3)
abline(a = 0, b = 1, lty = 2)


## -------------------------------------------------------------------------------
J <- sensitivity + specificity - 1 
threshold[which.max(J)]     # optimal threshold
sensitivity[which.max(J)]   # sensitivity at optimal threshold
specificity[which.max(J)]   # specificity at optimal threshold


## ----bankruptcy1, fig.width=7, fig.height=6, out.width="0.7\\textwidth", fig.cap="Scatterplot of financial ratios with points distinguished by observed classes."----
data(bankruptcy, package = "MixGHD")
X <- bankruptcy[,-1]
Class <- factor(bankruptcy$Y, levels = c(1:0), 
                labels = c("solvent", "bankrupt"))
cl <- clPairs(X, Class)
legend("bottomright", legend = cl$class, 
       pch = cl$pch, col = cl $col, inset = 0.02)


## -------------------------------------------------------------------------------
mod <- MclustDA(X, Class, modelType = "EDDA")
summary(mod)


## ----bankruptcy2, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Scatterplots of financial ratios with points distinguished by observed classes. Panel (a) shows the ellipses implied by the estimated model. Panel (b) includes black points corresponding to the misclassified observations."----
plot(mod, what = "scatterplot")
plot(mod, what = "error")


## ----echo=FALSE, eval=FALSE-----------------------------------------------------
## pred <- predict(mod)
## table(Class, Predicted = pred$classification)
## 
## # equal costs of misclassification
## (C <- 1 - diag(2))
## pred <- predict(mod, prop = mod$prop*rowSums(C))
## (tab <- table(Class, Predicted = pred$classification))
## sum(tab * C) # Total cost of misclassification


## -------------------------------------------------------------------------------
(C <- matrix(c(0,1,10,0), nrow = 2, ncol = 2, byrow = TRUE))


## -------------------------------------------------------------------------------
rowSums(C)


## -------------------------------------------------------------------------------
pred <- predict(mod)
(tab <- table(Class, Predicted = pred$classification))
sum(tab * C)


## -------------------------------------------------------------------------------
pred <- predict(mod, prop = mod$prop*rowSums(C))
(tab <- table(Class, Predicted = pred$classification))
sum(tab * C)


## ----simimbal1, echo=-1, fig.width=6, fig.height=5, out.width="0.49\\textwidth", fig.subcap = c("", ""), fig.cap="Histograms for synthetic datasets with observations sampled from two different Gaussian distributions. In the training data set, the cases ($y=1$) and the controls ($y=0$) are sampled in about the same proportions (a), whereas the cases ($y=1$) account for 10\\% of the observations in the whole population (b)."----
set.seed(20190515)
# generate training data from a balanced case-control sample
n_train <- 1000
class_train <- factor(sample(0:1, size = n_train, prob = c(0.5, 0.5), 
                             replace = TRUE))
x_train <- ifelse(class_train == 1, rnorm(n_train, mean = 3, sd = 1), 
                                    rnorm(n_train, mean = 0, sd = 1))

hist(x_train[class_train == 0], breaks = 11, xlim = range(x_train), 
     main = "", xlab = "x", 
     col = adjustcolor("dodgerblue2", alpha.f = 0.5), border = "white")
hist(x_train[class_train == 1], breaks = 11, add = TRUE,
     col = adjustcolor("red3", alpha.f = 0.5), border = "white")
box()

# generate test data from mixture f(x) = 0.9 * N(0,1) + 0.1 * N(3,1)
n <- 10000
mixpro <- c(0.9, 0.1)
class_test <- factor(sample(0:1, size = n, prob = mixpro, 
                            replace = TRUE))
x_test <- ifelse(class_test == 1, rnorm(n, mean = 3, sd = 1), 
                                  rnorm(n, mean = 0, sd = 1))
hist(x_test[class_test == 0], breaks = 15, xlim = range(x_test), 
     main = "", xlab = "x", 
     col = adjustcolor("dodgerblue2", alpha.f = 0.5), border = "white")
hist(x_test[class_test == 1], breaks = 11, add = TRUE,
     col = adjustcolor("red3", alpha.f = 0.5), border = "white")
box()


## -------------------------------------------------------------------------------
mod <- MclustDA(x_train, class_train)
summary(mod, parameters = TRUE)


## -------------------------------------------------------------------------------
pred <- predict(mod, newdata = x_test)
classError(pred$classification, class_test)$error
BrierScore(pred$z, class_test)


## ----simimbal, cache=TRUE-------------------------------------------------------
priorProp <- seq(0.01, 0.99, by = 0.01)
CE <- BS <- rep(as.double(NA), length(priorProp))
for (i in seq(priorProp))
{
  pred <- predict(mod, newdata = x_test, 
                  prop = c(1-priorProp[i], priorProp[i]))
  CE[i] <- classError(pred$classification, class = class_test)$error
  BS[i] <- BrierScore(pred$z, class_test)
}


## ----simimbal2, fig.width=6, fig.height=5, out.width="0.7\\textwidth",fig.cap="Classification error and Brier score as functions of the prior probability for the minority class. The vertical segments show the biased sample proportion of cases in the training set (dashed line), and the sample proportion of cases in the test set (dotted line), which is usually unknown."----
matplot(priorProp, cbind(CE,BS), type = "l", lty = 1, lwd = 2, xaxt = "n",
        xlab = "Class prior probability", ylab = "", ylim = c(0,max(CE,BS)), 
        col = c("red3", "dodgerblue3"),
        panel.first = 
          { abline(h = seq(0,1,by=0.05), col = "grey", lty = 3)
            abline(v = seq(0,1,by=0.05), col = "grey", lty = 3) 
          })
axis(side = 1, at = seq(0, 1, by = 0.1))
abline(v = mod$prop[2],             # training proportions
       lty = 2, lwd = 2)            
abline(v = mean(class_test == 1),   # test proportions (usually unknown)
       lty = 3, lwd = 2)   
legend("topleft", legend = c("ClassError", "BrierScore"),
       col = c("red3", "dodgerblue3"), lty = 1, lwd = 2, inset = 0.02)


## -------------------------------------------------------------------------------
(priorProbs <- classPriorProbs(mod, x_test))


## -------------------------------------------------------------------------------
pred <- predict(mod, newdata = x_test, prop = priorProbs)
classError(pred$classification, class = class_test)$error
BrierScore(pred$z, class_test)


## -------------------------------------------------------------------------------
(prior_test <- prop.table(table(class_test)))
pred <- predict(mod, newdata = x_test, prop = prior_test)
classError(pred$classification, class = class_test)$error
BrierScore(pred$z, class_test)


## -------------------------------------------------------------------------------
data(wdbc, package="mclust")
x <- with(wdbc, 
    0.2322*Texture_mean + 0.01117*Area_extreme + 68.37*Smoothness_extreme)
Class <- wdbc[,"Diagnosis"]
mod <- MclustDA(x, Class, modelType = "MclustDA")
summary(mod)


## ----wdbcuniv1, fig.width=6, fig.height=5, out.width="0.8\\textwidth", fig.cap="Densities for the benign and malignant tumors estimated using the univariate feature extracted from the breast cancer dataset."----
(prop <- mod$prop)
col <- mclust.options("classPlotColors")
x0 <- seq(0, max(x)*1.1, length = 1000)
par1 <- mod$models[["B"]]$parameters
f1 <- dens(par1$variance$modelName, data = x0, parameters = par1)
par2 <- mod$models[["M"]]$parameters
f2 <- dens(par2$variance$modelName, data = x0, parameters = par2)
matplot(x0, cbind(prop[1]*f1, prop[2]*f2), type = "l", lty = 1, 
        col = col, ylab = "Class density", xlab = "x")
legend("topright", title = "Diagnosis:", legend = names(prop), 
       col = col, lty = 1, inset = 0.02)


## ----echo=-1--------------------------------------------------------------------
set.seed(20190520)
cv <- cvMclustDA(mod)  # by default: prop = mod$prop
unlist(cv[c("error", "se.error")])


## ----echo=-1--------------------------------------------------------------------
set.seed(20190520)
cv <- cvMclustDA(mod, prop = c(0.5, 0.5))
unlist(cv[c("error", "se.error")])


## -------------------------------------------------------------------------------
x0 <- seq(min(x), max(x), length = 1000)
pred <- predict(mod, newdata = x0) 
(threshold1 <- approx(pred$z[,2], x0, xout = 0.5)$y)
pred <- predict(mod, newdata = x0, prop = c(0.5, 0.5))
(threshold2 <- approx(pred$z[,2], x0, xout = 0.5)$y)


## ----wdbcuniv2, echo=-1, fig.width=8, fig.height=4, out.width="\\textwidth", fig.cap="Distribution of the training data conditional on the true classes and on the predicted classes for the univariate feature extracted from the breast cancer data."----
par(mfrow=c(1,2))
plot(mod, what = "scatterplot", main = TRUE)
abline(v = threshold1, lty = 2)
plot(mod, what = "classification", main = TRUE)
abline(v = threshold1, lty = 2)


## ----eval=FALSE, echo=FALSE-----------------------------------------------------
## z <- predict(mod, newdata = x0)$z
## matplot(x0, z, type = "n", yaxt = "n")
## axis(side = 2, at = seq(0, 1, by = 0.1))
## rect(xleft = par("usr")[1], xright = par("usr")[2],
##      ybottom = 0.3, ytop = 0.7,
##      col = adjustcolor("grey", alpha.f = 0.5), border = NA)
## abline(h = seq(0, 1, by = 0.1), lty = 3)
## matpoints(x0, z, type = "l", lty = 1, col = col)
## legend("topright", title = "Diagnosis", legend = names(prop),
##        col = col, lty = 1, inset = 0.02)
## abline(v = threshold, lty = 2)


## ----wdbcuniv3, echo=-1, cache=TRUE, fig.width=6, fig.height=5, out.width="0.8\\textwidth", fig.cap="The cross-validated  misclassification error rate as a function of the probability threshold for the univariate feature extracted from the breast cancer data."----
set.seed(1)
threshold <- seq(0.1, 0.9, by = 0.05)
ngrid <- length(threshold)
cv <- data.frame(threshold, error = numeric(ngrid))
cverr <- cvMclustDA(mod, verbose = FALSE)
for (i in seq(threshold))
{
  cv$error[i] <- classError(ifelse(cverr$z[,2] > threshold[i], "M", "B"),
                            Class)$errorRate
}  
min(cv$error)
threshold[which.min(cv$error)]

ggplot(cv, aes(x = threshold, y = error)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(0,1,by=0.1)) +
  ylab("CV misclassification error") + 
  xlab("Threshold probability of class (M)")


## ----wdbcuniv4, echo=-1, cache=TRUE, fig.width=6, fig.height=5, out.width="0.8\\textwidth", fig.cap="Plot of the CV misclassification error rate as a function of the class prior probability, with error bars shown at $\\pm$ one standard error of the CV procedure, for the univariate feature extracted from the breast cancer data."----
set.seed(1)
priorProb <- seq(0.1, 0.9, by = 0.05)
ngrid <- length(priorProb)
cv_error2 <- data.frame(priorProb, cv = numeric(ngrid), 
                        lower = numeric(ngrid), upper = numeric(ngrid))
for (i in seq(priorProb))
{
  cv <- cvMclustDA(mod, prop = c(1-priorProb[i], priorProb[i]), 
                   verbose = FALSE)
  cv_error2$cv[i]    <- cv$ce
  cv_error2$lower[i] <- cv$ce - cv$se.ce
  cv_error2$upper[i] <- cv$ce + cv$se.ce
}  
min(cv_error2$cv)
priorProb[which.min(cv_error2$cv)]

ggplot(cv_error2, aes(x = priorProb, y = cv, ymin = lower, ymax = upper)) +
  geom_point() + geom_linerange() +
  scale_x_continuous(breaks = seq(0,1,by=0.1)) +
  ylab("CV misclassification error") + 
  xlab("Class (M) prior probability")


## ----ssc_example, echo=FALSE, fig.width=5.5, fig.height=5, out.width="0.49\\textwidth", fig.subcap = rep("",2), fig.cap="Example of two-class simulated dataset with (a) all labeled data (points are marked according to the true classes) and (b) only partial knowledge of labeled data (unlabeled data are shown as grey squares)."----
n <- 200
pars <- list(pro = c(0.5, 0.5),
             mean = matrix(c(-1,1), nrow = 2, ncol = 2, byrow = TRUE),
             variance = mclustVariance("EII", d = 2, G = 2))
pars$variance$sigmasq <- 1
data <- sim("EII", parameters = pars, n = n, seed = 12)
class <- data[,1]
X <- data[,-1]
clPairs(X, class, symbols = c(1,2))

# Randomly remove labels
cl <- class; cl[sample(1:n, size = 195)] <- NA
# table(cl, useNA = "ifany")
clPairs(X, ifelse(is.na(cl), 0, class),
        symbols = c(0, 16, 17), colors = c("grey", 4, 2))


## ----ssc_example2, echo=FALSE, fig.width=5.5, fig.height=5, out.width="0.7\\textwidth", fig.cap="Classification boundaries for the two-class simulated dataset obtained (i) under the assumption of full knowledge of class labels (cyan solid line), (ii) using only the labeled data (black dashed line), and (iii) both labeled and unlabeled data (black dotted line)."----
mod      <- MclustDA(X, class, modelType = "EDDA") 
mod_EDDA <- MclustDA(X[!is.na(cl),], cl[!is.na(cl)], modelType = "EDDA", modelNames = "EII")
mod_SSC  <- MclustSSC(X, cl)

ngrid = 100
xgrid = seq(-4, 4, length.out = ngrid)
ygrid = seq(-5, 5, length.out = ngrid)
xygrid = expand.grid(xgrid, ygrid)

pred_mod  <- predict(mod, newdata = xygrid)
pred_EDDA <- predict(mod_EDDA, newdata = xygrid)
pred_SSC  <- predict(mod_SSC, newdata = xygrid)

col = mclust.options("classPlotColors")[class]
pch = class
pch[!is.na(cl)] = ifelse(cl[!is.na(cl)] == 1, 19, 17)

plot(X, pch = pch, col = col)
contour(xgrid, ygrid, matrix(pred_mod$z[,1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lwd = 2, col = "cyan2")
contour(xgrid, ygrid, matrix(pred_EDDA$z[,1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lty = 2, lwd = 2)
contour(xgrid, ygrid, matrix(pred_SSC$z[,1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lty = 3, lwd = 2)


## -------------------------------------------------------------------------------
data(olive, package = "pgmm")
X <- olive[,3:10]
class <- factor(olive$Region, levels = 1:3, 
                labels = c("South", "Sardinia", "North"))
table(class)


## -------------------------------------------------------------------------------
mod_EDDA_full <- MclustDA(X, class, modelType = "EDDA")
pred_EDDA_full <- predict(mod_EDDA_full, newdata = X)
classError(pred_EDDA_full$classification, class)$errorRate
BrierScore(pred_EDDA_full$z, class)


## ----echo=-1--------------------------------------------------------------------
set.seed(20200807)
pct_labeled_data <- 10
n <- nrow(X)
cl <- class
cl[sample(1:n, round(n*(1-pct_labeled_data/100)))] <- NA
table(cl, useNA = "ifany")


## ----olive_ssc, cache=TRUE, out.width="70%", fig.cap="BIC values for the semi-supervised classification models fitted to the Italian olive oils data using 10\\% of labeled data."----
mod_SSC <- MclustSSC(X, cl)
plot(mod_SSC, what = "BIC")
mod_SSC$BIC
pickBIC(mod_SSC$BIC, 5) - max(mod_SSC$BIC)  # BIC diff for the top-5 models


## -------------------------------------------------------------------------------
summary(mod_SSC)


## -------------------------------------------------------------------------------
pred_SSC <- predict(mod_SSC, newdata = X[is.na(cl),])
table(Predicted = pred_SSC$classification, Class = class[is.na(cl)])
classError(pred_SSC$classification, class[is.na(cl)])$errorRate
BrierScore(pred_SSC$z, class[is.na(cl)])


## ----olive_ssc_sim, echo=-1, cache=TRUE, fig.width=7, fig.height=5, out.width="80%", fig.cap="Brier score values for the EDDA classification model on labeled data and the semi-supervised classification (SSC) model as a function of the percentage of labeled data."----
set.seed(20201013)
pct_labeled_data <- c(5, seq(10, 90, by = 10), 95)
BS <- matrix(as.double(NA), nrow = length(pct_labeled_data), ncol = 2,
             dimnames = list(pct_labeled_data, c("EDDA", "SSC")))
for (i in seq(pct_labeled_data))
{
  cl <- class
  labeled <- sample(1:n, round(n*pct_labeled_data[i]/100))
  cl[-labeled] <- NA
  # Classification on labeled data
  mod_EDDA  <- MclustDA(X[labeled,], cl[labeled], 
                        modelType = "EDDA")
  # prediction for the unlabeled data
  pred_EDDA <- predict(mod_EDDA, newdata = X[-labeled,])
  BS[i,1]   <- BrierScore(pred_EDDA$z, class[-labeled])
  # Semi-supervised classification
  mod_SSC  <- MclustSSC(X, cl)
  # prediction for the unlabeled data
  pred_SSC <- predict(mod_SSC, newdata = X[-labeled,])
  BS[i,2]  <- BrierScore(pred_SSC$z, class[-labeled])
}
BS

matplot(pct_labeled_data, BS, type = "b", 
        lty = 1, pch = c(19, 15), col = c(2,4), xaxt = "n",
        xlab = "Percentage of labeled data", ylab = "Brier score")
axis(side = 1, at = pct_labeled_data)
abline(h = BrierScore(pred_EDDA_full$z, class), lty = 2)
legend("topright", pch = c(19, 15), col = c(2,4), lty = 1, 
       legend = c("EDDA", "SSC"), inset = 0.02)

