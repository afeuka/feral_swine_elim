### Title: testing hypergeometric distribution in nimble
### Author: Abbey Feuka
### Date: 14AUG24
### Notes:

require(nimble)

# dhypergeo <- nimbleFunction(
#   run = function(x = double(), #number of successful balls drawn
#                  N = double(), #total number of balls (all kinds)
#                  m = double(), #number of desired ball type in pop
#                  n = double(), #number of draws
#                  log = integer(0, default=0)) {
#     returnType(double())
# 
#     xprob <- (factorial(m)/(factorial(x)*factorial(m-x))) *
#               (factorial(N-m)/(factorial(n-x)*factorial((N-m)-(n-x)))) /
#                 (factorial(N)/(factorial(n)*factorial(N-n)))
# 
#     if(log){
#       return(log(xprob))
#     } else {
#       return(xprob)
#     }
#   }
# )
dhypergeo <- nimbleFunction(
  run = function(x = double(), #number of white balls drawn
                 m = double(), #number of white ball type in pop
                 b = double(), #number of black ball type in pop
                 k = double(), #number of draws
                 log = integer(0, default=0)) {
    returnType(double())
  
    xprob <- (factorial(m)/(factorial(x)*factorial(m-x))) *
                  (factorial(b)/(factorial(k-x)*factorial(b-(k-x)))) /
                    (factorial(m+b)/(factorial(k)*factorial((m+b)-k)))
    
    if(log){
      return(log(xprob))
    } else {
      return(xprob)
    }
  }
)

x<-35
m<-60
b<-100
k<-50

dhyper(x=x,m=m,n=b,k=k)
dhypergeo(x=x,m=m,b=b,k=k)
(factorial(m)/(factorial(x)*factorial(m-x))) *
  (factorial(b)/(factorial(k-x)*factorial(b-(k-x)))) /
  (factorial(m+b)/(factorial(k)*factorial((m+b)-k)))

dat <- rhyper(300,m=m,n=b,k=k)

hist(rnbinom(10000,size=100,prob=0.01))
# fit model ---------------------------
## model specification -------------------------------------------------
hypergeom_code <- nimbleCode({
  
  m ~ dnbinom(1,500)
  b ~ dnbinom(1,500)
  k ~ dnbinom(1,500)
  
  for(i in 1:nsamp){
    y[i] ~ dhypergeo(m=m,b=b,k=k)
  }
  
})

## nimble lists  -------------------------------------------------
modDat <- list(y=dat)

### constants -----------------------------
const <- list(nsamp=length(dat))


### initial values -----------------------------
# inits <- list()

## nimble configuration -------------------------------------------------
mod <- nimbleModel(code = hypergeom_code,
                   data = modDat,
                   constants = const)
# mod$initializeInfo()

mcmc.conf <- configureMCMC(mod,enableWAIC=F)

Cmcmc <- buildMCMC(mcmc.conf,resetFuncitions=T) #uncompiled MCMC
Cmod <- compileNimble(mod,Cmcmc) #compiled model

## run mcmc ----------------------------------- 
niter <- 10000
samples <- runMCMC(Cmod$Cmcmc,
                   niter = niter,
                   thin= 1,
                   nburn = niter*0.2,
                   setSeed=1,
                   nchains = 1)

samples <- as.data.frame(samples)

plot(samples[,"m"],typ="l")
abline(h=m,col="red")
plot(samples[,"b"],typ="l")   
abline(h=b,col="red")
plot(samples[,"k"],typ="l")   
abline(h=k,col="red")
