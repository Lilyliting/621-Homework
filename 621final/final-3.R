
# Dealing with data -------------------------------------------------------

library(xlsx)
library(plotly)
library(quantmod)
SPX <- read.xlsx("SPX.xls", sheetName = "QuoteData-1", header=F) # read data
#SPX <- read.xlsx("AAPL.xls", sheetName = "QuoteData-1", header=F) # read data
SPX <- SPX[colSums(!is.na(SPX)) > 0] # omit NA columns
SPX <- SPX[rowSums(!is.na(SPX)) > 0, ] # omit NA rows
date0 <- as.numeric(as.character(SPX[1, 1])) 
price0 <-as.numeric(as.character(SPX[1, 2])) 
r0 <- as.numeric(as.character(SPX[1, 3]))/100 
div <- 0

colnames(SPX) <- c("Date", "Tm", "K", "Price")
SPX <- SPX[-c(1, 2), ]
for (i in 1:4) { # Convert factor to numeric
    SPX[, i] <- as.numeric(as.character(SPX[, i]))
}
SPX <- SPX[!duplicated(SPX[, 2:3], fromLast = T), ]

head(SPX)

# a -----------------------------------------------------------------------
Option_BSM <- function(S0, K, Tm, sigma, r=r0, div=0) {
    # Pricing by Black Scholes
    d1 <- (log(S0/K) + (r - div + sigma^2/2)*Tm)/sigma/sqrt(Tm)
    d2 <- d1 - sigma*sqrt(Tm)
    p <- S0*exp(-div*Tm)*pnorm(d1) - K*exp(-r*Tm)*pnorm(d2)
    return(p)
}

fsigma <- function(sigma, Tm, K, price) {
    # Epsilon. To calculate implied vol
    S0 <- price0
    price.by.bs <- Option_BSM(S0, K, Tm, sigma)
    ans <- price.by.bs - price
    return(ans)
}

imp_vol <- function(Tm, K, price) {
    # use bisection method to calculate implied vols
    a <- 0.0001
    b <- 10
    epsilon <- abs(a - b)
    while(epsilon > 1e-4) {
        mid <- (a + b)/2
        if (fsigma(mid, Tm, K, price)*fsigma(b, Tm, K, price) < 0 ) a <- mid
        else b <- mid
        epsilon <- abs(a - b)
    }
    return(a)
}

Implied_Vol <- SPX[, 1]
for(i in 1:nrow(SPX)) {
    price <- SPX[i, 4]
    Tm <- SPX[i, 2]
    K <- SPX[i, 3]
    Implied_Vol[i] <- imp_vol(Tm, K, price)
}
head(Implied_Vol)

newdf <- data.frame(SPX[, -1], Implied_Vol)
newdf <- subset(newdf, Implied_Vol>0.0001) #eliminate meaningless implied vol
head(newdf)


# To find 20 different strike prices and 4 different maturities

todr <- order(table(newdf$Tm), decreasing = T)[1:4]
T1 <- as.numeric(names(table(newdf$Tm)[todr]))

K1 <- newdf$K[newdf$Tm==T1[1]]
K2 <- newdf$K[newdf$Tm==T1[2]]
K3 <- newdf$K[newdf$Tm==T1[3]]
K4 <- newdf$K[newdf$Tm==T1[4]]
temp <- intersect(intersect(K1,K2),intersect(K3,K4))
l <- length(temp)
common_element=temp[(floor(l/2)-9):(floor(l/2)+10)]

dfselect <- subset(newdf, Tm %in% T1 & K %in% common_element)
a <- order(dfselect[, 1], dfselect[, 2])
dfselect <- dfselect[a, ]

if (nrow(dfselect) < 80) {print("No enough option data")}

plot_ly(x=dfselect$K, y=dfselect$Tm, z=dfselect$Implied_Vol)
# b -----------------------------------------------------------------------

library(akima)
library(rgl)
x <- dfselect$K 
y <- dfselect$Tm 
z <- dfselect$Implied_Vol

# linear interpolation
n_interpolation <- 200
spline_interpolated <- interp(x, y, z,
                              xo=seq(min(x), max(x), length = n_interpolation),
                              yo=seq(min(y), max(y), length = n_interpolation),
                              linear = T)

x.si <- spline_interpolated$x
y.si <- spline_interpolated$y
z.si <- spline_interpolated$z
plot_ly(x=x.si, y=y.si, z=z.si, type = "surface")

# c -----------------------------------------------------------------------

# NO


# d -----------------------------------------------------------------------
# local volatility

local_vol <- function(Tm, K, price) {
    # Square of local volatility
    S0 <- price0
    deltaT <- 0.001*Tm
    deltaK <- 0.001*K
    sig <- imp_vol(Tm, K, price)
    d1 <- (log(S0/K) + (r0 - div + sig^2/2)*(Tm))/sig/sqrt(Tm)
    dT1 <- (imp_vol(Tm+deltaT, K, price) - imp_vol(Tm, K, price))/deltaT
    dK1 <- (imp_vol(Tm, K+deltaK, price) - imp_vol(Tm, K, price))/deltaK
    dK2 <- (imp_vol(Tm, K+deltaK, price) - 2*imp_vol(Tm, K, price)
            + imp_vol(Tm, K-deltaK, price))/deltaK^2
    numerator <- 2*sig*dT1*Tm + sig^2 + 2*sig*(r0-div)*Tm*K*dK1
    denominator <- (1 +K*d1*dK1*sqrt(Tm))^2 + K^2*Tm*sig*(dK2 - d1*dK2^2*sqrt(Tm))
    lsigma <- numerator/denominator
    return(lsigma)
}

Local_Vol <- dfselect[, 1]
for(i in 1:nrow(dfselect)) {
    # calculate local vols
    price <- dfselect[i, 3]
    Tm <- dfselect[i, 1]
    K <- dfselect[i, 2]
    Local_Vol[i] <- local_vol(Tm, K, price)
}
head(Local_Vol)
dfselect <- cbind(dfselect, Local_Vol)

# e -----------------------------------------------------------------------

# Solve the PDE

# European Call option - Explicit Finite Difference method
Option_Ex <- function(Tm, K, b, N, Nj) {
    # precompute constants
    dt <- Tm/N
    dx <- 0.2
    
    # initialise asset prices at maturity
    St <- seq(1,2*Nj+1)
    St <- price0 + dx*(St-1-Nj)
    
    # initialise option values at maturity
    C <- matrix(0, ncol = (N + 1), nrow = (2*Nj + 1))
    C[, N+1] <- pmax(C[, N+1], St - K)
    
    # step back 
    for (i in N:1) {
        for(j in (2-i+N):(2*Nj+i-N)) {
            dCdS <- (C[j+1, i+1] - 2*C[j, i+1] + C[j-1, i+1])/dx^2
            dCdt <- -b/2*dCdS
            C[j, i] <- C[j, i+1] - dCdt*dt
        }
    }
    ans <- C[Nj+1, 1]
    return(ans)
}

PricebyD <- dfselect[, 1]
N <- 500
Nj <- 300
for(k in 1:nrow(dfselect)) {
    # calculate option price by Dupire’s local vol
    Tm <- dfselect[k, 1]
    K <- dfselect[k, 2]
    b <- dfselect[k, 5]
    PricebyD[i] <- Option_Ex(Tm, K, b, N, Nj)
}
PricebyD

# f -----------------------------------------------------------------------
dfselect <- cbind(dfselect, PricebyD)
write.xlsx(dfselect, "SPXvolatility.xls", row.names = F) # write xlsx file
# write.xlsx(dfselect, "AAPLvolatility.xls", row.names = F)

# g -----------------------------------------------------------------------

# Prepare Data
# Get option chains from Yahoo finance
AAPL.OPTS <- getOptionChain("AAPL", NULL)

D0 <- as.numeric(as.Date("2017/05/15", format="%Y/%m/%d"))

C1 <- AAPL.OPTS$May.19.2017$calls
D1 <- as.numeric(as.Date("2017/05/19", format="%Y/%m/%d"))
T1 <- D1-D0
C1$Date <- D1
C1$Tm <- T1/365

C2 <- AAPL.OPTS$May.26.2017$calls
D2 <- as.numeric(as.Date("2017/05/26", format="%Y/%m/%d"))
T2 <- D2-D0
C2$Date <- D2
C2$Tm <- T2/365

C3 <- AAPL.OPTS$Jun.02.2017$calls
D3 <- as.numeric(as.Date("2017/06/02", format="%Y/%m/%d"))
T3 <- D3-D0
C3$Date <- D3
C3$Tm <- T3/365

C4 <- AAPL.OPTS$Jun.09.2017$calls
D4 <- as.numeric(as.Date("2017/06/09", format="%Y/%m/%d"))
T4 <- D4-D0
C4$Date <- D4
C4$Tm <- T4/365

C5 <- AAPL.OPTS$Jun.16.2017$calls
D5 <- as.numeric(as.Date("2017/06/16", format="%Y/%m/%d"))
T5 <- D5-D0
C5$Date <- D5
C5$Tm <- T5/365

temp <- rbind(C1, C2, C3, C4, C5)
temp$price <- (temp$Bid + temp$Ask)/2

AAPL <- temp[, c(8, 9, 1, 10)]

# Get the actual stock price
todaystock <- getQuote("AAPL") 
S_0 <- todaystock[, 2]

r <- 0.0075
# Treasury bill rate

row0 <- c(D0, S_0, r, NaN)
AAPL <- rbind(row0, c("Date", "T", "K", "Price"), AAPL)
write.xlsx(AAPL, "AAPL.xls", sheetName = "QuoteData-1", row.names = F, col.names = F) # write xlsx file
