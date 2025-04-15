### Title:  Nimble function for n area correction
### Author: Abbey Feuka
### Date: 18JUN24
### Notes: 
###################
library(nimble)

#double(0) <- scalar input/output
#double(1) <- vector input/output

avail_fun <- nimbleFunction(
  run = function(pass = double(0), #removal pass number
                 # type = double(0), #removal type index
                 gamma = double(0), #area of removal impact/area of site (time t)
                 gamma_past = double(1), #area of removal impact/area of site (time t-1)
                 theta = double(0), #capture probability (time t)
                 theta_past = double(1)) { #capture probability (time)
    returnType(double(0))
    if(pass!=1) {
      return(gamma*theta*prod((1-gamma_past)+gamma_past*(1-theta_past)))
    } else {
      return(gamma*theta)
    }
  })
