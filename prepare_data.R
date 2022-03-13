

nt <- 20
gridPoints <- seq(0, 1, length.out = nt)


nknots <- 10


##Number of facilities:
n <- length(unique(dat$facility_id))


## Time-invariant data, one row per subject:
data.id <- dat[!duplicated(dat$subject_id), ]

##Number of subjects within each facility:
nf <- table(x = data.id$facility_id)

n.subj <- c(as.vector(nf))

##Number of observations for each subject: 
nts <- table(x = dat$subject_id)

nt.subj <- c(0, as.vector(nts))


Time.e <- data.id$time_surv


## 15-point Gauss-Kronrod rule for the integrals within the survival model

wk <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785, 
                0.209482141084728, 0.190350578064785, 0.140653259715526, 
                0.0630920926299786, 0.0229353220105292, 0.10479001032225, 
                0.169004726639268, 0.204432940075299, 0.204432940075299, 
                0.169004726639268, 0.10479001032225, 0.0229353220105292)
sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397, 
        0, 0.405845151377397, 0.741531185599394, 0.949107912342758, 
        -0.991455371120813, -0.864864423359769, -0.586087235467691, 
        -0.207784955007898, 0.207784955007898, 0.586087235467691, 
        0.864864423359769, 0.991455371120813)

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)

P <- Time.e/2
st <- outer(P, sk + 1) 

############################################
# Functions for P-splines

# truncated p-th power function
tpower <- function(x, t, p)
  (x - t) ^ p * (x > t)

# B-spline basis with degree deg
bbase1 <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B 
}

############################################
# P-splines for the coefficients in the longitudinal model 

t.hi = floor(max(Time.e)) +1
t.lo = 0
n.seg = nknots - 2
deg = 2


Betabs = bbase1(gridPoints, t.lo, t.hi, n.seg, deg)
BetaBs <- Betabs

Dbeta <- diag(ncol(BetaBs))

priorTau.betaX0 <- crossprod(diff(Dbeta, diff = 2)) + 1e-06 * Dbeta
priorTau.betaX1 <- crossprod(diff(Dbeta, diff = 2)) + 1e-06 * Dbeta
priorTau.betaX2 <- crossprod(diff(Dbeta, diff = 2)) + 1e-06 * Dbeta


priorTau.betaZ1 <- crossprod(diff(Dbeta, diff = 2)) + 1e-06 * Dbeta
priorTau.betaZ2 <- crossprod(diff(Dbeta, diff = 2)) + 1e-06 * Dbeta

############################################
# P-splines for the coefficients in the survival model 

t.hi = floor(max(Time.e)) +1
t.lo = 0
n.seg = nknots-2
deg = 2

Gammabs = bbase1(Time.e, t.lo, t.hi, n.seg, deg)
GammaBs <- Gammabs

Gammas = bbase1(c(t(st)), t.lo, t.hi, n.seg, deg)
GammaBs.surv <- Gammas


Dgamma <- diag(ncol(GammaBs))

priorTau.Bs.h0 <- crossprod(diff(Dgamma, diff = 2)) + 1e-06 * Dgamma

priorTau.gammaX1 <- crossprod(diff(Dgamma, diff = 2)) + 1e-06 * Dgamma
priorTau.gammaX2 <- crossprod(diff(Dgamma, diff = 2)) + 1e-06 * Dgamma

priorTau.gammaZ1 <- crossprod(diff(Dgamma, diff = 2)) + 1e-06 * Dgamma
priorTau.gammaZ2 <- crossprod(diff(Dgamma, diff = 2)) + 1e-06 * Dgamma

priorTau.alphas <- crossprod(diff(Dgamma, diff = 2)) + 1e-06 * Dgamma

############################################


###Initials:

Bs.betaX0 <- rep(0, nknots)

Bs.betaX1 <- rep(0, nknots)

Bs.betaX2 <- rep(0, nknots)

Bs.betaZ1 <- rep(0, nknots)

Bs.betaZ2 <- rep(0, nknots)


Bs.h0 <- rep(0, nknots)

Bs.gammaX1 <- rep(0, nknots)
Bs.gammaX2 <- rep(0, nknots)

Bs.gammaZ1 <- rep(0, nknots)
Bs.gammaZ2 <- rep(0, nknots)

Bs.alpha <- rep(0, nknots)


############################################
# Data for the MCMC

Data <- list(n = n, 
             K = K,
             event = data.id$survival, 
             T1 = data.id$time_surv,
             y = dat$y, 
             x1 = data.id$x1,
             x2 = data.id$x2,
             z1 = data.id$z1,
             z2 = data.id$z2,
             zeros = rep(0, dim(data.id)[1]),
             wk = wk,
             GammaBs.surv = GammaBs.surv,
             GammaBs = GammaBs,
             priorMean.Bs.betaX0 = Bs.betaX0,
             priorTau.betaX0 = priorTau.betaX0,
             BetaBs = BetaBs,
             priorMean.Bs.betaX1 = Bs.betaX1,
             priorTau.betaX1 = priorTau.betaX1,
             priorMean.Bs.betaX2 = Bs.betaX2,
             priorTau.betaX2 = priorTau.betaX2,
             priorMean.Bs.betaZ1 = Bs.betaZ1,
             priorTau.betaZ1 = priorTau.betaZ1,
             priorMean.Bs.betaZ2 = Bs.betaZ2,
             priorTau.betaZ2 = priorTau.betaZ2,
             priorA = 2,
             priorB = 0.5,
             mu0 = 0,
             mu1 = 0,
             priorMean.Bs.alpha = Bs.alpha,
             priorTau.alpha = priorTau.alphas, 
             priorMean.Bs.gammaX1 = Bs.gammaX1,
             priorTau.gammaX1 = priorTau.gammaX1, 
             priorMean.Bs.gammaX2 = Bs.gammaX2,
             priorTau.gammaX2 = priorTau.gammaX2, 
             priorMean.Bs.gammaZ1 = Bs.gammaZ1,
             priorTau.gammaZ1 = priorTau.gammaZ1, 
             priorMean.Bs.gammaZ2 = Bs.gammaZ2,
             priorTau.gammaZ2 = priorTau.gammaZ2, 
             priorMean.Bs.h0 = Bs.h0,
             priorTau.Bs.h0 = priorTau.Bs.h0,
             ncF = nknots,
             nt = nt,
             nt.subj = nt.subj,
             fac.id = rep(1:n, n.subj),
             dim_data.id = dim(data.id)[1]
)

############################################


### Parameters to keep track of during the MCMC,

parms <- c( "Bs.betaX0", "Bs.betaX1", "Bs.betaX2", "Bs.betaZ1",  "Bs.betaZ2", "Bs.alpha", "Bs.h0", "Bs.gammaX1", "Bs.gammaX2", "Bs.gammaZ1", "Bs.gammaZ2", "sigmasq_S", "sigmasq_F")

