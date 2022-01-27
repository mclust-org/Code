## ----fishery, echo=FALSE, fig.width = 9, fig.height = 7, width="\\textwidth", fig.cap="Distribution of fish lengths for a sample of snappers with estimated Gaussian mixtures from 1 up to 4 number of mixture components. Dashed lines represent component densities, solid lines the mixture densities."----
data(Snapper, package = "FSAdata")
x <- Snapper[,1]

mod1 <- densityMclust(x, G = 1, plot = FALSE)
mod2 <- densityMclust(x, G = 2, plot = FALSE)
mod3 <- densityMclust(x, G = 3, plot = FALSE)
mod4 <- densityMclust(x, G = 4, plot = FALSE)

par(mfrow = c(2,2), mar = c(2,2,1,1), oma = c(2,2,0,0))
x0 <- extendrange(x, f = 0.1)
x0 <- seq(x0[1], x0[2], length = 1000)
#
cdens <- predict(mod1, newdata = x0, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*mod1$parameters$pro))
plot(mod1, x, what = "density", lwd = 2, breaks = 20, ylim = c(0,0.4))
matplot(x0, t(cdens), type = "l", lty = 2, col = 1, add = TRUE)
# title("G = 1", cex.main = 0.9)
#
cdens <- predict(mod2, newdata = x0, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*mod2$parameters$pro))
plot(mod2, x, what = "density", lwd = 2, breaks = 20, ylim = c(0,0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
# title("G = 2", cex.main = 0.9)
#
cdens <- predict(mod3, newdata = x0, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*mod3$parameters$pro))
plot(mod3, x, what = "density", lwd = 2, breaks = 20, ylim = c(0,0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
# title("G = 3", cex.main = 0.9)
#
cdens <- predict(mod4, newdata = x0, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*mod4$parameters$pro))
plot(mod4, x, what = "density", lwd = 2, breaks = 20, ylim = c(0,0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
# title("G = 4", cex.main = 0.9)
#
mtext("Fish length (in)", side = 1, line = 1, outer = TRUE)
mtext("Density", side = 2, line = 1, outer = TRUE)

