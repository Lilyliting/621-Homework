library(quantmod)
library(Sim.DiffProc)

setwd("/Users/apple/Desktop/621final")


# 1 -----------------------------------------------------------------------


# XLF
# Financial Select Sector SPDR Fund

# Top 20 tickers:
# JPM	JP Morgan Chase & Co
# BRK-B	Berkshire Hathaway B
# WFC	Wells Fargo & Co
# BAC	Bank of America Corp
# C	Citigroup Inc
# GS	Goldman Sachs Group Inc
# USB	US Bancorp
# CB	Chubb Limited
# MS	Morgan Stanley
# AXP	American Express Co
# PNC	PNC Finl Services Group
# AIG	American Intl Group Inc
# MET	Metlife Inc
# BK	The Bank of New York Mellon Corp
# SCHW	Schwab Charles Corp
# BLK	BlackRock Inc
# PRU	Prudential Financial Inc
# CME	CME Group Inc A
# COF	Capital One Financial
# MMC	Marsh & McLennan Companies


# Download data
sbls <- c("JPM", "BRK-B", "WFC", "BAC", "C", "GS", "USB", "CB", "MS", "AXP", 
          "PNC", "AIG", "MET", "BK", "SCHW", "BLK", "PRU", "CME", "COF", "MMC")
tickers <- data.frame(matrix(ncol = 21, nrow = 1347))
ReturnMatrix <- data.frame(matrix(ncol = 20, nrow = 1347))

for (i in 1:20) {
    sbl <- sbls[i]
    temp <- getSymbols(Symbols = sbl, from = "2012-01-01", to = "2017-05-10", auto.assign = F)
    ReturnMatrix[, i] <- (Cl(temp)-Op(temp))/Op(temp)
    tickers[, (i+1)] <- Cl(temp)
}

tickers[, 1] <- rownames(as.data.frame(temp))
colnames(tickers) <- c("Date", sbls)
colnames(ReturnMatrix) <- sbls
head(tickers)
# All prices are above $5

# Principal component analysis
tickers.pca <- prcomp(ReturnMatrix, scale. = T)
plot(tickers.pca, type = "l")
pca.summary <- summary(tickers.pca)
pca.summary
# PC7

importance <- as.data.frame(pca.summary$importance)
weight <- importance[2, 1:7]
atic <- abs(tickers.pca$rotation)
stdize <- sweep(atic, 2, colSums(atic), "/")

# Weighted sum of PC1-PC7
weightsum <- atic[, 1:7] %*% matrix(as.numeric(weight), ncol=1)
odr <- order(weightsum, decreasing = T)
order.pca <- stdize[odr, ]
rownames(order.pca)[1:4]

# 2 -----------------------------------------------------------------------

n <- c(1,3,5,7)
fit.data <- ts(tickers[, (n+1)], deltat=1/255)

model.match <- as.data.frame(matrix(nrow = 4, ncol = 2))
model.match[, 1] <- colnames(fit.data)
best.n <- c()
for (i in 1:4) {
    # model A
    fx1 <- expression( theta[1]*x )
    gx1 <- expression( theta[2]*x )
    mod1 <- fitsde(data=fit.data[, i], drift=fx1, diffusion=gx1,
                   start = list(theta1=.1, theta2=.1), pmle="euler")
    
    # model B
    fx2 <- expression( theta[1]+theta[2]*x )
    gx2 <- expression( theta[3]*x^theta[4] )
    mod2 <- fitsde(data=fit.data[, i], drift=fx2, diffusion=gx2,
                   start = list(theta1=.1,theta2=.1,theta3=.1,theta4=.1), pmle="euler")
    
    # model C
    fx3 <- expression( theta[1]*x )
    gx3 <- expression( theta[2]+theta[3]*x^theta[4] )
    mod3 <- fitsde(data=fit.data[, i], drift=fx3, diffusion=gx3,
                   start = list(theta1=.1,theta2=.1,theta3=.1,theta4=.1), pmle="euler")
    
    # model D
    fx4 <- expression( theta[1]*x )
    gx4 <- expression( theta[2]*x^(3/2) )
    mod4 <- fitsde(data = fit.data[, i], drift=fx4, diffusion=gx4,
                   start = list(theta1=.1,theta2=.1), pmle="euler")
    
    # model E
    fx5 <- expression( theta[1]+theta[2]*x )
    gx5 <- expression( (theta[3]+theta[4]*log(x))*x )
    mod5 <- fitsde(data=fit.data[, i], drift=fx5, diffusion=gx5,
                   start = list(theta1=.1,theta2=.1,theta3=.1,theta4=.1), pmle="euler")
    
    # Computes AIC
    AIC <- c(AIC(mod1),AIC(mod2),AIC(mod3),AIC(mod4),AIC(mod5))
    k <- which.min(AIC)
    best.n[i] <- k
    best <- paste0("model", k)
    model.match[i, 2] <- best
}
model.match # All model 1

# Coefficients
coefs <- data.frame(matrix(nrow = 4, ncol = 3))
coefs[, 1] <- model.match[, 1]

for (i in 1:4) {
    # model A
    fx1 <- expression( theta[1]*x )
    gx1 <- expression( theta[2]*x )
    mod1 <- fitsde(data=fit.data[, i], drift=fx1, diffusion=gx1,
                   start = list(theta1=.1, theta2=.1), pmle="euler")
    coefs[i, 2:3] <- coef(mod1)
}

coefs

# 3 -----------------------------------------------------------------------
corrmatrix <- cor(fit.data)
corrmatrix

# 4 -----------------------------------------------------------------------

S0 <- fit.data[1, ]
monte_carlo_corr <- function(S0, Tm=1, dt=1/255, n=1000, corr=corrmatrix) {
    # S0: Initial stock prices
    # Tm: Time to maturity
    # dt: Time interval
    # n:  Quantity of paths
    # corr: Correlation matrix
    R <- chol(corrmatrix)
    L <- t(R)
    STj <- data.frame(matrix(ncol = 4, nrow = n))
    for (j in 1:n) {
        dZt <- matrix(rnorm(4*Tm/dt)*sqrt(dt), nrow = 4)
        dWt <- L %*% dZt
        
        St <- S0
        
        for (i in 1:(Tm/dt)) {
            drift <- St*coefs[, 2]*dt
            diffusion <- coefs[, 3]*dWt[, i]
            diffusion[1] <- diffusion[1]*St[1]
            diffusion[2] <- diffusion[2]*St[2]
            diffusion[3] <- diffusion[3]*St[3]
            diffusion[4] <- diffusion[4]*St[4]
            St <- St + drift + diffusion
        }
        STj[j, ] <- St
    }
    return(STj)
}
STj <- monte_carlo_corr(S0)
STj

statistics <- do.call(data.frame, 
                      list(mean = apply(STj, 2, mean),
                           sd = apply(STj, 2, sd),
                           skewness = apply(STj, 2, skewness),
                           kurtosis = apply(STj, 2, kurtosis)))
statistics

# 5 -----------------------------------------------------------------------
XLF <- getSymbols(Symbols = "XLF", from = "2012-01-01", to = "2017-05-10", auto.assign = F)
XLF.close <- ts(as.numeric(XLF$XLF.Close), deltat=1/255)

# model A - geometric Brownian motion
fx1 <- expression( theta[1]*x )
gx1 <- expression( theta[2]*x )
mod1 <- fitsde(data=XLF.close, drift=fx1, diffusion=gx1,
               start = list(theta1=.1, theta2=.1), pmle="euler")
coef <- coef(mod1)

S0 <- XLF.close[1]
monte_carlo <- function(S0, Tm=1, dt=1/255, n=1000) {
    ST <- data.frame(matrix(ncol = 1, nrow = n))
    for (j in 1:n) {
        dWt <- rnorm(Tm/dt)*sqrt(dt)
        St <- S0
        
        for (i in 1:(Tm/dt)) {
            drift <- St*coef[1]*dt
            diffusion <- St*coef[2]*dWt[i]
            St <- St + drift + diffusion
        }
        ST[j, 1] <- St
    }
    return(ST)
}
XLF.MC <- monte_carlo_corr(S0)
XLF.MC[, 1]

MCsimulations <- cbind(STj, XLF.MC[, 1])
colnames(MCsimulations) <- c("JPM", "WFC", "C", "USB", "XLF")
write.csv(MCsimulations, file = "MCsimulations.csv", row.names = F)

# 6 -----------------------------------------------------------------------
Multivariate.data <- as.data.frame(cbind(fit.data, XLF.close))
write.csv(Multivariate.data, "P2_Monte_Carlo.csv", row.names = F)
lm.fit <- lm(XLF.close~. + 0, data = Multivariate.data)
w <- lm.fit$coefficients
w

# 7 -----------------------------------------------------------------------

basket <- MCsimulations # Weighted ST
for (i in 1:4) {
    basket[, i] <- MCsimulations[, i]*w[i]
}

mean(pmax(rowSums(basket) - MCsimulations[, 5], 0))
mean(pmax(MCsimulations[, 5] - rowSums(basket), 0))
