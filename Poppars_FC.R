# Estimation of population parameters (recruitment, survival, abundance)
#   using the Schwarz-Arnason model

library(jagsUI)
library(IPMbook)  # for helper functions
library(mcmcOutput)  # for plotting functions

# Load and prepare data
# =====================

sims <- read.csv("FC_3years.csv")
str(sims)
# 'data.frame':	127 obs. of  3 variables:
# $ y1: int  1 1 1 1 1 1 1 1 1 1 ...
# $ y2: int  1 0 0 0 1 0 1 0 0 0 ...
# $ y3: int  0 0 0 0 0 0 0 0 0 0 ...

y <- as.matrix(sims)
( nInd <- nrow(y) )   # 127
( nYears <- ncol(y) ) # 3

# Prepare for the JAGS run
# ========================
# Augment the data
nAug <- 50
yAug <- rbind(y, matrix(0, nAug, nYears))
jdata <- list(y = yAug, nInd = nInd, nYears = nYears, M = nInd+nAug,
    z = zKnown(yAug), w=rep(c(1, NA), c(nInd, nAug)))

wanted <- c("b", "phi", "p", "omega", "N", "B", "Nsuper")

# null model, b(.)phi(.)p(.)
# --------------------------

( out_b.phi.p. <- jags(jdata, NULL, wanted, "SA_b(.)phi(.)p(.).jags", DIC=FALSE,
    n.chains=3, n.iter=1e5, n.adapt=1000, n.thin=1, parallel=TRUE) )  # 1 mins for 1e5

#           mean     sd    2.5%     50%   97.5% overlap0 f  Rhat  n.eff
# b[1]     0.289  0.042   0.210   0.287   0.377    FALSE 1 1.001   2533
# b[2]     0.356  0.021   0.312   0.357   0.395    FALSE 1 1.001   2533
# b[3]     0.356  0.021   0.312   0.357   0.395    FALSE 1 1.001   2533
# phi      0.350  0.070   0.228   0.345   0.502    FALSE 1 1.000 300000
# p        0.758  0.088   0.608   0.749   0.941    FALSE 1 1.000  24162
# omega    0.888  0.074   0.733   0.897   0.995    FALSE 1 1.000  17770
# N[1]    45.150  6.028  35.000  45.000  58.000    FALSE 1 1.000  16574
# N[2]    72.925  8.067  58.000  73.000  89.000    FALSE 1 1.000  38893
# N[3]    81.460  9.774  64.000  81.000 101.000    FALSE 1 1.000  28610
# B[1]    45.150  6.028  35.000  45.000  58.000    FALSE 1 1.000  16574
# B[2]    54.728  5.652  45.000  55.000  66.000    FALSE 1 1.000  16150
# B[3]    58.109  5.219  49.000  58.000  69.000    FALSE 1 1.000  13128
# Nsuper 157.986 12.636 132.000 159.000 177.000    FALSE 1 1.000  16315


plot(out_b.phi.p.)  # a bit poor mixing, n.eff's big

# Check that augmentation is adequate
max(out_b.phi.p.$sims.list$Nsuper)  # 177
nInd+nAug                           # 177

#save(out_b.phi.p., file="out_b(.)phi(.)p(.).RData")
#load("out_b(.)phi(.)p(.).RData")
# b(t)phi(t)p(.) model
# --------------------

( out_bt_phit_p. <- jags(jdata, NULL, wanted, "SA_b(t)phi(t)p(.).jags", DIC=FALSE,
    n.chains=3, n.iter=1e6, n.adapt=1000, parallel=TRUE) )  # 1 mins for 1e5
    
#           mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# b[1]     0.288  0.044   0.207   0.286   0.378    FALSE 1 1.000 83361
# b[2]     0.334  0.048   0.242   0.333   0.428    FALSE 1 1.000  9834
# b[3]     0.379  0.046   0.288   0.378   0.470    FALSE 1 1.000 19675
# phi[1]   0.460  0.112   0.264   0.452   0.700    FALSE 1 1.000  7429
# phi[2]   0.291  0.080   0.157   0.284   0.468    FALSE 1 1.000 12982
# p        0.752  0.088   0.602   0.742   0.937    FALSE 1 1.001  1786
# omega    0.892  0.074   0.736   0.902   0.996    FALSE 1 1.001  1576
# N[1]    45.591  6.266  35.000  45.000  59.000    FALSE 1 1.001  2503
# N[2]    74.035  8.926  59.000  74.000  92.000    FALSE 1 1.001  2345
# N[3]    81.522  9.676  64.000  82.000 100.000    FALSE 1 1.001  1661
# B[1]    45.591  6.266  35.000  45.000  59.000    FALSE 1 1.001  2503
# B[2]    52.879  6.054  43.000  52.000  66.000    FALSE 1 1.000  9061
# B[3]    60.116  6.217  49.000  60.000  73.000    FALSE 1 1.001  1981
# Nsuper 158.586 12.558 133.000 160.000 177.000    FALSE 1 1.001  1491

plot(out_bt_phit_p.)
# Convergence ok, some n.eff's too small

( mco <- mcmcOutput(out_bt_phit_p.))
crosscorrPlot(mco, c("b", "phi"))

b <- mco$b
plot(b[,1], b[,2])  # positive correlation
plot(b[,2], b[,3])  # no correlation

#save(out_bt_phit_p., file="out_b(t)phi(t)p(.).RData")

# Fully time-dependent model, b(t)phi(t)p(t)
# -----------------------------------------=
# Not all parameters identifiable

( out_bt_phit_pt <- jags(jdata, NULL, wanted, "SA_b(t)phi(t)p(t).jags", DIC=FALSE,
    n.chains=3, n.iter=1e6, n.adapt=10000, n.thin=10, parallel=TRUE) )  # 26 mins for 1e6

#           mean     sd    2.5%     50%   97.5% overlap0 f Rhat  n.eff
# b[1]     0.319  0.099   0.180   0.298   0.562    FALSE 1    1  32476
# b[2]     0.329  0.098   0.128   0.329   0.518    FALSE 1    1  65585
# b[3]     0.353  0.073   0.220   0.348   0.505    FALSE 1    1   6392
# phi[1]   0.495  0.133   0.274   0.482   0.796    FALSE 1    1  60699
# phi[2]   0.281  0.092   0.141   0.267   0.498    FALSE 1    1  14896
# p[1]     0.698  0.185   0.341   0.710   0.986    FALSE 1    1  28285
# p[2]     0.707  0.130   0.462   0.708   0.942    FALSE 1    1  21389
# p[3]     0.787  0.131   0.528   0.796   0.989    FALSE 1    1   8496
# omega    0.916  0.061   0.776   0.928   0.997    FALSE 1    1 102505
# N[1]    52.149 16.551  34.000  47.000  95.000    FALSE 1    1  44337
# N[2]    79.692 14.898  58.000  77.000 114.000    FALSE 1    1  18890
# N[3]    78.977 13.850  61.000  76.000 111.000    FALSE 1    1   6530
# B[1]    52.149 16.551  34.000  47.000  95.000    FALSE 1    1  44337
# B[2]    53.563 15.412  22.000  53.000  85.000    FALSE 1    1 116164
# B[3]    57.335 11.250  40.000  55.000  83.000    FALSE 1    1   6538
# Nsuper 163.047 10.355 140.000 165.000 177.000    FALSE 1    1 103725

plot(out_bt_phit_pt)

plot(out_bt_phit_pt, c("phi[10]", "p[11]", "b[1]", "p[1]"))
plot(out_bt_phit_pt, c("b[1]", "b[2]", "b[11]", "phi[10]"))

max(out_bt_phit_pt$sims.list$Nsuper)  # 177
nInd+nAug                             # 177

( mco <- mcmcOutput(out_bt_phit_pt))
crosscorrPlot(mco, c("b", "phi", "p"), main="Cross-correlations for b(t)phi(t)p(t)")

#save(out_bt_phit_pt, file="out_b(t)phi(t)p(t).RData")




# This is the best model for Non-spatial CJS in MLE
# MLE phi = 0.27
# Bayes phi = 0.35

( out_b.phi.pt <- jags(jdata, NULL, wanted, "SA_b(.)phi(.)p(t).jags", DIC=FALSE,
                         n.chains=3, n.iter=1e5, n.adapt=10000, n.thin=10, parallel=TRUE) )  # 1 mins for 1e5

#           mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# b[1]     0.295  0.074   0.180   0.285   0.465    FALSE 1 1.001  2731
# b[2]     0.352  0.037   0.267   0.358   0.410    FALSE 1 1.001  2731
# b[3]     0.352  0.037   0.267   0.358   0.410    FALSE 1 1.001  2731
# phi      0.354  0.075   0.228   0.348   0.518    FALSE 1 1.000 12275
# p[1]     0.747  0.162   0.426   0.761   0.989    FALSE 1 1.000  3950
# p[2]     0.750  0.095   0.572   0.749   0.933    FALSE 1 1.000 30000
# p[3]     0.748  0.122   0.530   0.742   0.976    FALSE 1 1.000 30000
# omega    0.904  0.066   0.760   0.915   0.996    FALSE 1 1.000 30000
# N[1]    47.180 11.642  34.000  44.000  77.000    FALSE 1 1.000  5636
# N[2]    73.535  8.259  59.000  73.000  90.000    FALSE 1 1.000 16481
# N[3]    82.808 12.989  62.000  82.000 109.000    FALSE 1 1.000 30000
# B[1]    47.180 11.642  34.000  44.000  77.000    FALSE 1 1.000  5636
# B[2]    54.847  8.183  40.000  55.000  71.000    FALSE 1 1.000 12824
# B[3]    58.887  7.919  46.000  58.000  75.000    FALSE 1 1.000  7571
# Nsuper 160.914 11.197 137.000 162.000 177.000    FALSE 1 1.000 30000

plot(out_bt_phit_pt)
#save(out_b.phi.pt, file="out_b(.)phi(.)p(t).RData")


