# Problem 1 -----------------------------------
library(quantmod)
library(ggplot2)
#library(rgl)

# a)
S0 <- 100       # Stock price
tau <- 30/252   # Time to maturity
K <- 100        # Strike price
r <- 0.075      # Interest rate
sigma <- 0.2    # Volatility

Callprice <- function(S0, tau, K, r, sigma) {
    d1 <- (log(S0/K)+(r+sigma^2/2)*tau)/sigma/sqrt(tau)
    d2 <- d1 - sigma*sqrt(tau)
    c <- S0*pnorm(d1) - K*exp(-r*tau)*pnorm(d2)
    return(c)
}

Putprice <- function(S0, tau, K, r, sigma) {
    d1 <- (log(S0/K)+(r+sigma^2/2)*tau)/sigma/sqrt(tau)
    d2 <- d1 - sigma*sqrt(tau)
    p <- K*exp(-r*tau)*pnorm(-d2) - S0*pnorm(-d1)
    return(p)
}



Callprice(S0, tau, K, r, sigma)
# 3.049619
Putprice(S0, tau, K, r, sigma)
# 2.456149

# b)
# Put-Call parity


# c)
# Get option chains from Yahoo finance
AAPL.OPTS <- getOptionChain("AAPL", NULL)
C1 <- AAPL.OPTS$Mar.17.2017$calls
C1$Ave.Price <- (C1$Bid + C1$Ask)/2
C1 <- C1[, c(1, 8)]
C2 <- AAPL.OPTS$Apr.21.2017$calls
C2$Ave.Price <- (C2$Bid + C2$Ask)/2
C2 <- C2[, c(1, 8)]
C3 <- AAPL.OPTS$Jul.21.2017$calls
C3$Ave.Price <- (C3$Bid + C3$Ask)/2
C3 <- C3[, c(1, 8)]

temp <- merge(C1, C2, by = "Strike")
calls <- merge(temp, C3, by = "Strike")
colnames(calls) <- c("Strike", "Mar.17.2017", "Apr.21.2017", "Oct.20.2017")
calls <- calls[1:20, ]
print(calls)

# Get the actual stock price
todaystock <- getQuote("AAPL") 
S_0 <- todaystock[, 2]

r <- 0.0075
# Treasury bill rate
# https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=billrates
fsigma <- function(sigma, K_i, maturity_i) {
    cc <- calls[K_i, maturity_i + 1]
    K <- calls[K_i, 1]
    if (maturity_i == 1) tau = 23/252 # time to maturity
    else if(maturity_i == 2) tau = 48/252
    else tau = 111/252
    ans <- Callprice(S_0, tau, K, r, sigma) - cc
    return(ans)
}

interval <- calls
interval[, 2:4] <- NaN
# To find a interval [a, a+delta] that makes the secant method converge
delta <- 0.1
for(i in 1:20) {
    for(j in 1:3) {
        a <- seq(1, 5, by=delta)
        for(k in a) {
            if (fsigma(k,i,j)*fsigma(k + delta,i,j) < 0) {interval[i, j+1] <- k}
        }
    }
}

ImpliedVolatility <- calls
ptm <- proc.time()
count <- 0
for(i in 1:20) {
    for(j in 1:3) {
        a <- interval[i, j + 1]
        if (is.na(a) == T) {
            ImpliedVolatility[i, j + 1] <- NaN
            next
            }
        b <- a + delta
        epsilon <- abs(a - b)
        while(epsilon > 1e-4) {
            count <- count +1
            mid <- (a + b)/2
            if(fsigma(a, i, j)*fsigma(mid, i, j) < 0 ) b <- mid
            else a <- mid
            epsilon <- abs(a - b)
        }
        ImpliedVolatility[i, j+1] <- a
    }
}
count
ImpliedVolatility
proc.time() - ptm

# d)
# Secant
ImpliedVolatility2 <- calls
ptm <- proc.time()
count2 <- 0
for(i in 1:20) {
    for(j in 1:3) {
        x1 <- interval[i, j + 1]
        if (is.na(x1) == T) {
            ImpliedVolatility2[i, j + 1] <- NaN
            next
        }
        x2 <- x1 + delta
        dfxn <- (fsigma(x1, i, j) - fsigma(x2, i, j))/(x1 - x2)
        tang <- fsigma(x2, i, j)/dfxn
        epsilon <- abs(tang)
        while(epsilon > 1e-4) {
            count2 <- count2 + 1
            x1 <- x2
            x2 <- x2 - tang
            dfxn <- (fsigma(x1, i, j) - fsigma(x2, i, j))/(x1 - x2)
            tang <- fsigma(x2, i, j)/dfxn
            epsilon <- abs(tang)
        }
        ImpliedVolatility2[i, j+1] <- x2
    }
}
count2
ImpliedVolatility2
proc.time() - ptm

# e)
temp1 <- as.matrix(ImpliedVolatility)
temp2 <- rbind(temp1[, 1:2], temp1[, c(1, 3)], temp1[, c(1, 4)])
maturity <- c(rep("Mar.17.2017", 20), rep("Apr.21.2017", 20), rep("Oct.20.2017", 20))
newdf <- data.frame(df1, maturity)
colnames(newdf)[2] <- "Implied_Volatility"

ggplot(data = newdf, aes(x = Strike, y = Implied_Volatility, colour = maturity)) + geom_point()


plot3d(ImpliedVolatility$Strike, 1, ImpliedVolatility$Mar.17.2017, col = "blue", ylim = c(.5, 3.5))
points3d(ImpliedVolatility$Strike, 2, ImpliedVolatility$Apr.21.2017, col = "green")
points3d(ImpliedVolatility$Strike, 3, ImpliedVolatility$Oct.20.2017, col = "red")

# f)

Delta1 <- function(d, S0, tau, K, r, sigma, ty) {
    if (ty == "b") D <- (Callprice(S0, tau, K, r, sigma) - Callprice(S0 - d, tau, K, r, sigma))/d
    else if (ty == "f") D <- (Callprice(S0 + d, tau, K, r, sigma) - Callprice(S0, tau, K, r, sigma))/d
    else if (ty == "c") D <- (Callprice(S0 + d, tau, K, r, sigma) - Callprice(S0 - d, tau, K, r, sigma))/d/2
    else D <- NaN
    return(D)
}

Vega1 <- function(d, S0, tau, K, r, sigma, ty) {
    if (ty == "b") V <- (Callprice(S0, tau, K, r, sigma) - Callprice(S0, tau, K, r, sigma - d))/d
    else if (ty == "f") V <- (Callprice(S0, tau, K, r, sigma + d) - Callprice(S0, tau, K, r, sigma))/d
    else if (ty == "c") V <- (Callprice(S0, tau, K, r, sigma + d) - Callprice(S0, tau, K, r, sigma - d))/d/2
    else V <- NaN
    return(V)
}

Gamma1 <- function(d, S0, tau, K, r, sigma, ty) {
    if (ty == "b") {G <- (Callprice(S0, tau, K, r, sigma) - 2*Callprice(S0 - d, tau, K, r, sigma) + 
                              Callprice(S0 - 2*d, tau, K, r, sigma))/d^2}
    else if (ty == "f") {G <- (Callprice(S0 + 2*d, tau, K, r, sigma) - 2*Callprice(S0 + d, tau, K, r, sigma) + 
                                   Callprice(S0, tau, K, r, sigma))/d^2}
    else if (ty == "c") {G <- (Callprice(S0 + d, tau, K, r, sigma) - 2*Callprice(S0, tau, K, r, sigma) + 
                                       Callprice(S0 - d, tau, K, r, sigma))/d^2}
    else G <- NaN
    return(G)
}

# for delta
d1 <- seq(1, .01,by=-0.01)
fd_f <- sapply(d1, Delta1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "f")
fd_b <- sapply(d1, Delta1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "b")
fd_c <- sapply(d1, Delta1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "c")

type <- c(rep("forward", length(d1)), rep("backword", length(d1)), rep("central", length(d1)))

delta_dd <- data.frame("d" = rep(d1, 3), "delta" = c(fd_f, fd_b, fd_c), "type" = type)
ggplot(data = delta_dd, aes(x = d, y = delta, colour = type)) + geom_line()

# for vega
fv_f <- sapply(d1, Vega1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "f")
fv_b <- sapply(d1, Vega1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "b")
fv_c <- sapply(d1, Vega1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "c")

vega_dd <- data.frame("d" = rep(d1, 3), "vega" = c(fv_f, fv_b, fv_c), "type" = type)
ggplot(data = vega_dd, aes(x = d, y = vega, colour = type)) + geom_line()

# for gamma
fg_f <- sapply(d1, Gamma1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "f")
fg_b <- sapply(d1, Gamma1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "b")
fg_c <- sapply(d1, Gamma1, S0=S0, tau=tau, K=K, r=r, sigma=sigma, "c")

gamma_dd <- data.frame("d" = rep(d1, 3), "gamma" = c(fg_f, fg_b, fg_c), "type" = type)
ggplot(data = gamma_dd, aes(x = d, y = gamma, colour = type)) + geom_line()

d <- 1e-5
Delta1(d, S0, tau, K, r, sigma, "f") #forward
Delta1(d, S0, tau, K, r, sigma, "b") #backward
Delta1(d, S0, tau, K, r, sigma, "c") #central

Vega1(d, S0, tau, K, r, sigma, "f") #forward
Vega1(d, S0, tau, K, r, sigma, "b") #backword
Vega1(d, S0, tau, K, r, sigma, "c") #central

Gamma1(d, S0, tau, K, r, sigma, "f") #forward
Gamma1(d, S0, tau, K, r, sigma, "b") #backword
Gamma1(d, S0, tau, K, r, sigma, "c") #central
# g)

delta_df <- ImpliedVolatility
vega_df <- ImpliedVolatility
gamma_df <- ImpliedVolatility

for (i in 1:20) {
    for (j in 1:3) {
        sigma <- ImpliedVolatility[i, j+1]
        if (is.na(sigma) == T) {
            delta_df[i, j+1] <- NaN
            vega_df[i, j+1] <- NaN
            gamma_df[i, j+1] <- NaN
            break}
        delta_df[i, j+1] <- Delta1(d, S_0, tau, K, r, sigma, "c")
        vega_df[i, j+1] <- Vega1(d, S_0, tau, K, r, sigma, "c")
        gamma_df[i, j+1] <- Gamma1(d, S_0, tau, K, r, sigma, "c")
    }
}
delta_df
vega_df
gamma_df


