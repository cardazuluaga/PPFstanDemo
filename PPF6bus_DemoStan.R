rm(list=ls()) 

library(ggplot2)
library(rstan)
require(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model <- "
// Inferring the state variables in a IEEE 6 bus test system
data {
  int<lower=1> N; // number of nodes
  int<lower=1> n; // number of unknown variables
  row_vector[n] y;
  matrix[n, n] Sigma;
  real y21; real y22; real y23; real y24; real y25; real y26;
  real y31; real y32; real y33; real y34; real y35; real y36;
  real y41; real y42; real y43; real y44; real y45; real y46;
  real y51; real y52; real y53; real y54; real y55; real y56;
  real y61; real y62; real y63; real y64; real y65; real y66;
  real phi21; real phi22; real phi23; real phi24; real phi25; real phi26;
  real phi31; real phi32; real phi33; real phi34; real phi35; real phi36;
  real phi41; real phi42; real phi43; real phi44; real phi45; real phi46;
  real phi51; real phi52; real phi53; real phi54; real phi55; real phi56;
  real phi61; real phi62; real phi63; real phi64; real phi65; real phi66;
  //real pd4; real pd5; real pd6;
  //real qd4; real qd5; real qd6;
  //real pg2; real pg3;
  real V1; real V2; real V3;
  real th1;
  real rho; real R; real Cp;
  real vci; real vco; real vr;
  real pi;
}
parameters {
  real <lower=-1,upper=1> th2;  real <lower=-1,upper=1> th3;
  real <lower=-1,upper=1> th4;  real <lower=-1,upper=1> th5;
  real <lower=-1,upper=1> th6;
  real <lower=0.2,upper=1.8>  V4;  real <lower=0.2,upper=1.8>  V5;
  real <lower=0.2,upper=1.8>  V6;
  //real <lower=0,upper=1> pd4;
  real <lower=5,upper=50> vwind;
  real <lower=0,upper=1> pd5; 
  real <lower=0,upper=1> pd6;
  real <lower=0,upper=1> qd4; 
  real <lower=0,upper=1> qd5; 
  real <lower=0,upper=1> qd6;
  real <lower=0,upper=1> pg2; 
  real <lower=0,upper=1> pg3;
}
transformed parameters {
  row_vector[n] Gx[1];
  real pd4;

  if (vwind <= vci) {
            pd4 <- 0;
  } else if (vci <= vwind && vwind <= vr) {
            pd4 <- 0.5*rho*pi*R^2*Cp*vwind^3;
  } else if (vr <= vwind && vwind <= vco){
            pd4 <- 0.5*rho*pi*R^2*Cp*vr^3;
  } else {
            pd4 <- 0;
  }

  Gx[1,1] <- V2 * (V1*y21*cos(th2 - th1 - phi21) + V2 * y22 * cos(th2 - th2 - phi22) + V3*y23*cos(th2 - th3 - phi23) + V4*y24*cos(th2 - th4 - phi24) + V5*y25*cos(th2 - th5 - phi25) + V6*y26*cos(th2 - th6 - phi26)) - pg2;
  Gx[1,2] <- V3 * (V1*y31*cos(th3 - th1 - phi31) + V2 * y32 * cos(th3 - th2 - phi32) + V3*y33*cos(th3 - th3 - phi33) + V4*y34*cos(th3 - th4 - phi34) + V5*y35*cos(th3 - th5 - phi35) + V6*y36*cos(th3 - th6 - phi36)) - pg3;
  Gx[1,3] <- V4 * (V1*y41*cos(th4 - th1 - phi41) + V2 * y42 * cos(th4 - th2 - phi42) + V3*y43*cos(th4 - th3 - phi43) + V4*y44*cos(th4 - th4 - phi44) + V5*y45*cos(th4 - th5 - phi45) + V6*y46*cos(th4 - th6 - phi46)) + pd4;
  Gx[1,4] <- V5 * (V1*y51*cos(th5 - th1 - phi51) + V2 * y52 * cos(th5 - th2 - phi52) + V3*y53*cos(th5 - th3 - phi53) + V4*y54*cos(th5 - th4 - phi54) + V5*y55*cos(th5 - th5 - phi55) + V6*y56*cos(th5 - th6 - phi56)) + pd5;
  Gx[1,5] <- V6 * (V1*y61*cos(th6 - th1 - phi61) + V2 * y62 * cos(th6 - th2 - phi62) + V3*y63*cos(th6 - th3 - phi63) + V4*y64*cos(th6 - th4 - phi64) + V5*y65*cos(th6 - th5 - phi65) + V6*y66*cos(th6 - th6 - phi66)) + pd6;
  Gx[1,6] <- V4 * (V1*y41*sin(th4 - th1 - phi41) + V2 * y42 * sin(th4 - th2 - phi42) + V3*y43*sin(th4 - th3 - phi43) + V4*y44*sin(th4 - th4 - phi44) + V5*y45*sin(th4 - th5 - phi45) + V6*y46*sin(th4 - th6 - phi46)) + qd4;
  Gx[1,7] <- V5 * (V1*y51*sin(th5 - th1 - phi51) + V2 * y52 * sin(th5 - th2 - phi52) + V3*y53*sin(th5 - th3 - phi53) + V4*y54*sin(th5 - th4 - phi54) + V5*y55*sin(th5 - th5 - phi55) + V6*y56*sin(th5 - th6 - phi56)) + qd5;
  Gx[1,8] <- V6 * (V1*y61*sin(th6 - th1 - phi61) + V2 * y62 * sin(th6 - th2 - phi62) + V3*y63*sin(th6 - th3 - phi63) + V4*y64*sin(th6 - th4 - phi64) + V5*y65*sin(th6 - th5 - phi65) + V6*y66*sin(th6 - th6 - phi66)) + qd6;
}
model {
  // Priors
  th2 ~ normal(0, sqrt(0.01));  th3 ~ normal(0, sqrt(0.01));
  th4 ~ normal(0, sqrt(0.01));  th5 ~ normal(0, sqrt(0.01));
  th6 ~ normal(0, sqrt(0.01));
  V4 ~ normal(1, sqrt(0.01));  V5 ~ normal(1, sqrt(0.01));  V6 ~ normal(1, sqrt(0.01));
  
  vwind ~ weibull(30,5);
  //pd4 ~ normal(0.7,sqrt(0.03));
  pd5 ~ normal(0.7,sqrt(0.03)); 
  pd6 ~ normal(0.7,sqrt(0.03));
  qd4 ~ normal(0.7,sqrt(0.03)); 
  qd5 ~ normal(0.7,sqrt(0.03));  
  qd6 ~ normal(0.7,sqrt(0.03));
  pg2 ~ normal(0.5,sqrt(0.03)); 
  pg3 ~ normal(0.6,sqrt(0.03));

  // Data Come From A Gaussian
  y ~ multi_normal(Gx, Sigma);
}"

N <- 1000

muy <- matrix(0,2,1)

Sigma <- 0.03 * diag(8)

n <- 8
y <- numeric(n)

y21 = 4.4721; y22 = 25.0497; y23 = 3.9223; y24 = 8.9443; y25 = 3.1623; y26 = 4.7193;
y31 = 0; y32 = y23; y33 = 17.1121; y34 = 0; y35 = 3.4922; y36 = 9.8058;
y41 = 4.8507; y42 = y24; y43 = y34; y44 = 15.9180; y45 = 2.2361; y46 = 0;
y51 = 3.2208; y52 = y25; y53 = y35; y54 = y45; y55 = 15.1641; y56 = 3.1623;
y61 = 0; y62 = y26; y63 = y36; y64 = y46; y65 = y56; y66 = 17.6169;

phi21 = 2.0344; phi22 = -1.1892; phi23 = 1.7682; phi24 = 2.0344; phi25 = 1.8925; phi26 = 1.9075;
phi31 = 0; phi32 = phi23; phi33 = -1.3255; phi34 = 0; phi35 = 2.0032; phi36 = 1.7682;
phi41 = 1.8158; phi42 = phi24; phi43 = phi34; phi44 = -1.1723; phi45 = 2.0344; phi46 = 0;
phi51 = 1.8314; phi52 = phi25; phi53 = phi35; phi54 = phi45; phi55 = -1.2142; phi56 = 1.8925;
phi61 = 0; phi62 = phi26; phi63 = phi36; phi64 = phi46; phi65 = phi56; phi66 = -1.3135;
V1 <- 1.05; V2 <- 1.05; V3 <- 1.07;
th1 <- 0;

rho <- 1.2235; R <- 45; Cp <- 0.473;
vci <- 10; vco <- 50; vr <- 25;

#Deterministic solution
V4t = 0.9858; V5t = 0.9804; V6t = 1.0052;
th2t = -0.0463; th3t = -0.0470; th4t = -0.0599; th5t = -0.0717; th6t = -0.0693; 

data <- list(y = y, 
             N = N, 
             n = n,
             y21 = y21, y22 = y22, y23 = y23, y24 = y24, y25 = y25, y26 = y26,
             y31 = y31, y32 = y32, y33 = y33, y34 = y34, y35 = y35, y36 = y36,
             y41 = y41, y42 = y42, y43 = y43, y44 = y44, y45 = y45, y46 = y46,
             y51 = y51, y52 = y52, y53 = y53, y54 = y54, y55 = y55, y56 = y56,
             y61 = y61, y62 = y62, y63 = y63, y64 = y64, y65 = y65, y66 = y66,
             phi21 = phi21, phi22 = phi22, phi23 = phi23, phi24 = phi24, phi25 = phi25, phi26 = phi26,
             phi31 = phi31, phi32 = phi32, phi33 = phi33, phi34 = phi34, phi35 = phi35, phi36 = phi36,
             phi41 = phi41, phi42 = phi42, phi43 = phi43, phi44 = phi44, phi45 = phi45, phi46 = phi46,
             phi51 = phi51, phi52 = phi52, phi53 = phi53, phi54 = phi54, phi55 = phi55, phi56 = phi56,
             phi61 = phi61, phi62 = phi62, phi63 = phi63, phi64 = phi64, phi65 = phi65, phi66 = phi66,
             rho = rho, R = R, Cp = Cp,
             vci = vci, vco = vco, vr = vr,
             V1 = V1, V2 = V2, V3 = V3,
             th1 = th1,
             Sigma = Sigma) # to be passed on to Stan

myinits <- list(
  list(th2 = 0, th3 = 0, th4 = 0, th5 = 0, th6 = 0, V4 = 1, V5 = 1, V6 = 1, V5 = 1, V6 = 1,pd4 = 0.7, pd5 = 0.7, pd6 = 0.5, qd4 = 0.7, qd5 = 0.7, qd6 = 0.7, pg2 = 0.5, pg3 = 0.6))

# parameters to be monitored: 
parameters <- c("th2", "th3", "th4", "th5", "th6", "V4", "V5", "V6")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=40000, 
                chains=1, 
                thin=1
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples)

th2 <- extract(samples)$th2
th3 <- extract(samples)$th3
th4 <- extract(samples)$th4
th5 <- extract(samples)$th5
th6 <- extract(samples)$th6
V4 <- extract(samples)$V4 
V5 <- extract(samples)$V5 
V6 <- extract(samples)$V6 

hist(th2)
hist(th3)
hist(th4)
hist(th5)
hist(th6)
hist(V4)
hist(V5)
hist(V6)
