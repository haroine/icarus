Xs=cbind(
  c(1,1,1,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,1,1,1),
  c(1,2,3,4,5,6,7,8,9,10)
)
# inclusion probabilities
piks=rep(0.2,times=10)
# vector of population totals
total=c(24,26,290)


# the g-weights using the truncated method
# g=sampling::calib(Xs,d=1/piks,total, method="logit", bounds=c(0.8,1.2))
# the calibration estimator of X is equal to 'total' vector
# tcal=t(g/piks)%*%Xs
# g

# the g-weights using the truncated method
# g2=calib(Xs,d=1/piks,total, method="logit", bounds=c(0.8,1.2))
# the calibration estimator of X is equal to 'total' vector
# tcal2=t(g/piks)%*%Xs
# g2
