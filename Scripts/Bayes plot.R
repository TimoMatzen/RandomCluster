#########################################################################
##########           Code for Plotting the BF             ###############
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Load packages
library(plotrix)

# Load the data
load("bayes5nod.Rdata")
bayes <- bayes5nod

# Set p.grid, equal to sampling grid
p.grid <- seq(0,.9, 0.05)

# Calculate the line means
meanl1 <- apply(bayes[,,1], 2, mean)
meanl2 <- apply(bayes[,,2], 2, mean)
meanl3 <- apply(bayes[,,3],2, mean)

## Calculate Credible intervals

# Calculate SDs
sdl1 <- apply(log(bayes[,,1]), 2, sd)
# Divide by calc ME
upperl1 <- log(meanl1) + (1.96 * (sdl1/sqrt(10)))
lowerl1 <- log(meanl1) - (1.96 * (sdl1/sqrt(10)))

# Calculate SDs
sdl2 <- apply(log(bayes[,,2]), 2, sd)
# Divide by calc ME
upperl2 <- log(meanl2) + (1.96 * (sdl2/sqrt(10)))
lowerl2 <- log(meanl2) - (1.96 * (sdl2/sqrt(10)))

sdl3 <- apply(log(bayes[,,3]), 2, sd)
# Divide by calc ME
upperl3 <- log(meanl3) + (1.96 * (sdl3/sqrt(10)))
lowerl3 <- log(meanl3) - (1.96 * (sdl3/sqrt(10)))

# Set grey color
greycol <- rgb(red = 190, green = 190, 
               blue = 190, alpha = 170, maxColorValue = 255)

# Set plotting margins
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# Plot first line
plot( x = p.grid, y = log(meanl3), xlim = c(0,.9), ylim = c(-5,25 ), xlab = "", 
     ylab = "", cex.lab = 1.3, cex.axis = 1.3, las = 1, yaxt = "n", bty = "n", 
     type = "b", pch = 2, bg = "grey", axes = FALSE, lwd = 2, 
     main = expression(paste( ~logBF[1][0], ~"for ", 
    "RC sampled data with 5 nodes", sep = " ")), 
     cex.main = 2)
# Make axis
axis(2, at = seq(-5, 20, 5), labels = seq(-5, 20, 5))
axis(1, at = c(0,p.grid), labels = c(0, p.grid), las = 2)
# Add confidence bands
polygon(c(p.grid, rev(p.grid)), c(upperl3,
                                  rev(lowerl3)),
        col="#7C060799",border=NA)
polygon(c(p.grid, rev(p.grid)), c(upperl2,
                                 rev(lowerl2)),
        col="#F1F1F199",border=NA)
polygon(c(p.grid, rev(p.grid)), c(upperl1,
                                 rev(lowerl1)),
        col="#1F28A299",border=NA)
# Plot BF for lower n
lines(x = p.grid, y = log(meanl2), type = "b", lwd = 2, col = "darkgrey")
lines(x = p.grid, y = log(meanl1), type = "b", lwd = 2, col = greycol, pch = 3)

ablineclip(h = 0, lty = 2, x2 = 1e+05, y2 = 0)
# Add texts and legend
mtext(expression(logBF[1][0]), side = 2, line = 2.5, las = 0, cex = 1.7)
mtext("P value", side = 1, 
      line = 3.5, las = 1, cex = 1.3)
text(0.6, 18, expression("Max" ~logBF[1][0] ~"= 4.75"), cex = 1.3)
text(.71, -3, "Evidence H0", pos = 4, cex = 1.3)
text(.71, 3, "Evidence H1", pos = 4, cex = 1.3)
arrows(.7, -1, .7, -4, length = 0.25, angle = 30, code = 2, lwd = 2)
arrows(.7, 1, .7, 4, length = 0.25, angle = 30, code = 2, lwd = 2)
legend(x = 0.5, 16, legend =  c(expression(logBF[1][0] ~"(n = 60)"), 
                                expression(logBF[1][0] ~"(n = 30)"),
                                expression(logBF[1][0] ~"(n = 10)")), 
       lwd = 3, col = c("black", "darkgrey", greycol), 
       bty = "n", x.intersp = 0.5, 
       cex = 1.2, pch =  c(2,1,3), lty= rep(NA, 3))




