
# Dependencies
library(TMB)
source("utils.R")

TMB::compile("msingarch_linear.cpp", "-O1 -g", DLLFLAGS = "") 
TMB::dynlib("msingarch_linear")

#TMB::compile("msingarch_log_linear.cpp", "-O1 -g", DLLFLAGS = "") 
TMB::dynlib("msingarch_log_linear")



# Specify model
modela <- MSingarch(m = 2, 
                    par = list(a = c(0.1, 0.4),
                               b = c(0.3, 0.4),
                               d = c(2, 3),
                               gamma = matrix(c(0.95, 0.05,0.05, 0.95),nrow = 2)),
                    mean_spec = "linear")

modelb <- MSingarch(m = 2, 
                    par = list(a = c(0.2, 0.4),
                               b = c(0.3, 0.4),
                               d = c(1, 0.5),
                               gamma = matrix(c(0.95, 0.05,0.05, 0.95),nrow = 2)),
                    mean_spec = "log-linear")


# Simulate a sample
n <- 500
set.seed(1)
sima <- rMSingarch(n = n, model = modela);sima$y
simb <- rMSingarch(n = n, model = modelb)

# Estimate data
fita <- fitMSingarch(data = sima$y,
                     model = modela, 
                     init_par = modela$par,
                     init_mean = sima$init_mean)
fita
plot(fita)

fitb <- fitMSingarch(data = simb$y,
                     model = modelb, 
                     init_par = modelb$par,
                     init_mean = simb$init_mean)
fitb
plot(fitb)

# Compare along with smoothing probabilities
library(GGally)

dfa <- data.frame(t = rep(seq(n),2), y = rep(fita$regimeprobs[ ,1],2),
                  lambda = c(fita$lambda_mat[,1], fita$lambda_mat[ ,2]),
                  regime = factor(c(rep("1",n), rep("2",n))))

smootha <- ggplot(dfa, aes(x = t, y = y)) +
  geom_line()
lambdasa <- ggplot(dfa, aes(x = t, y = lambda, color = regime)) +
  geom_line()
  
ggmatrix(list(plot(sima),
              plot(fita),
              smootha,
              lambdasa), nrow = 4, ncol = 1)


#####
dfb <- data.frame(t = rep(seq(n),2), y = rep(fitb$regimeprobs[ ,1],2),
                  lambda = c(fitb$lambda_mat[,1], fitb$lambda_mat[ ,2]),
                  regime = factor(c(rep("1",n), rep("2",n))))

smoothb <- ggplot(dfb, aes(x = t, y = y)) +
  geom_line()
lambdasb <- ggplot(dfb, aes(x = t, y = lambda, color = regime)) +
  geom_line()

ggmatrix(list(plot(simb),
              plot(fitb),
              smoothb,
              lambdasb), nrow = 4, ncol = 1)







ggmatrix(list(plot(simb),
              plot(fitb),
              smoothb), nrow = 3, ncol = 1)




# Test of stationary mean:

# Specify model
m <- 2
a = c(0.3, 0.3)
b = c(0.4, 0.2)
d = c(5, 2)
gamma = matrix(c(0.99, 0.01,0.01, 0.99),nrow = 2)
delta <- solve (t( diag (m)-gamma +1) ,rep (1,m))

delta %*% gamma
delta
sum(delta)

modela <- MSingarch(m = 2, 
                    par = list(a = a,
                               b = b,
                               d = d,
                               gamma = gamma),
                    mean_spec = "linear")



# Teorethical mean
tmumix <- sum((delta*d)/(1 - a))/(1 - sum(delta*b/(1 - a)))
A <-  cbind(b, b)
B1 <-  rbind(d*delta/(1 - a), d*delta/(1 - a))
B2 <-  rbind(b*delta/(1 - a), b*delta/(1 - a))
C <- matrix(nrow = 2, ncol = 2)
g1 <- gamma %*% solve(diag(2) - a[1]*gamma)
g2 <- gamma %*% solve(diag(2) - a[2]*gamma)
for(l in 1:2){
  for(k in 1:2){
    if(k == 1){C[l,k] = g1[l,k]}
    if(k == 2){C[l,k] = g2[l,k]}
  
  }
}
one <- sum(d*delta/(1 - a))
two <- sum(A*B1*C)
three <- sum(A*B2*C)
tmu <- (one + two)/(1 - three)


two <- three <- 0
for(l in 1:2){
  for(k in 1:2){
    if(k == 1){
      two <- two + b[k]*delta[l]*d[l]*g1[l,k]
      three <- three + b[k]*delta[l]*b[l]*g1[l,k]}
    if(k == 2){
      two <- two + b[k]*delta[l]*d[l]*g1[l,k]
      three <- three + b[k]*delta[l]*b[l]*g1[l,k]}
    
  }
}

tmu <- (one + two)/(1 - three)
tmu

# Simulate a sample
n <- 50000
sima <- rMSingarch(n = n, model = modela)

tmu
tmumix
mean(sima$y)

# Figure of intensity process
d/(1 - a - b)
plot(seq(n), sima$lambdavec[,2], type = "l")
lines(seq(n), sima$lambdavec[,1],col = "red")

# Variation in sample mean
R <- 100
n <- 10000
simmu <- rep(NA, R)
for(i in 1:R){
  sima <- rMSingarch(n = n, model = modela)
  
simmu[i] <- mean(sima$y)  
d/(1 - a - b)

colMeans(sima$lambdavec)
print(i)
}

mean(simmu)
sd(simmu)
tmu
hist(simmu)



# Test 2

colMeans(sima$lambdavec)
empmu <- mean(sima$y)

muvec <- d/(1-a) + empmu*b/(1 - a)
muvec
colMeans(sima$lambdavec)

sum(colMeans(sima$lambdavec)*delta)
mean(sima$lambda)
mean(sima$y)
