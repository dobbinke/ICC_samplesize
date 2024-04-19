# B is the number of biological replicates
# L is the number of labs
# R is the number of reps per cell in balanced design
# MCreps is the number of columns of the big matrix
myb0volfun <- function(B,L,R,MCreps,randseed) {
# This function returns the mean and SEM of the estimated volume.

 cols <- MCreps;
 lambda1sq <- sigmaesq+B*R*sigmalsq;
 lambda2sq <- sigmaesq+L*R*sigmabsq;
 lambda3sq <- sigmaesq;
 df1 <- L-1;
 df2 <- B-1;
 df3 <- B*L*R-B-L+1;

 # These should always be the same (?)
 set.seed(randseed);
 SLSQS <- lambda1sq * rchisq(rows,df=L-1) / (L-1);

 # These should be different to enhance reproducibility
 set.seed(runif(1,0,100000))
 SBSQS <- lambda2sq * rchisq(rows,df=B-1) / (B-1);
 SESQS <- lambda3sq * rchisq(rows,df=B*L*R-B-L+1)/ (B*L*R-B-L+1);

 Mmat <- matrix(rep(NA,rows*cols),ncol=cols);

 for (myr in 1:rows) {

  # These should also be the same, but not same as above
  set.seed(randseed+myr);
  W1 <- rchisq(cols,L-1);

  # these should be different
  set.seed(runif(1,0,1000000));
  W2 <- rchisq(cols,B-1);
  W3 <- rchisq(cols,B*L*R-B-L+1);

  thisMnum <-  ((B-1)*SBSQS[myr])/(L*R*W2) - ((B*L*R-B-L+1)*SESQS[myr])/(L*R*W3);
  belowzeros <- which(thisMnum<0);
  thisMnum[belowzeros] <- 0;
  thisMden <- ((L-1)*SLSQS[myr])/(B*R*W1) + ((B-1)*SBSQS[myr])/(L*R*W2) + ((B*L*R-B-L)*(B*L*R-B-L+1)*SESQS[myr])/(B*L*R*W3);
  Mmat[myr,] <- thisMnum/thisMden;

 }

 Lowers <- rep(NA,rows);
 Uppers <- rep(NA,rows);
 for (i in 1:rows) {
  Lowers[i] <- quantile(as.numeric(Mmat[i,]),alpha/2);
  Uppers[i] <- quantile(as.numeric(Mmat[i,]),1-alpha/2);
 }
 
 myret <- rep(NA,2);
 myret[1] <- mean(Uppers-Lowers);
 myret[2] <- sqrt(var(Uppers-Lowers)/rows);
 myret;
}


myl0volfun <- function(B,L,R,MCreps,randseed) {
 
 cols <- MCreps;
 lambda1sq <- sigmaesq+B*R*sigmalsq;
 lambda2sq <- sigmaesq+L*R*sigmabsq;
 lambda3sq <- sigmaesq;
 df1 <- L-1;
 df2 <- B-1;
 df3 <- B*L*R-B-L+1;

 # These should always be the same (?)
 set.seed(randseed);
 SBSQS <- lambda2sq * rchisq(rows,df=B-1) / (B-1);

 # These should be different to enhance reproducibility
 set.seed(runif(1,0,100000))
 SLSQS <- lambda1sq * rchisq(rows,df=L-1) / (L-1);
 SESQS <- lambda3sq * rchisq(rows,df=B*L*R-B-L+1)/ (B*L*R-B-L+1);

 Mmat <- matrix(rep(NA,rows*cols),ncol=cols);

 for (myr in 1:rows) {

  # These should also be the same, but not same as above
  set.seed(randseed+myr);
  W2 <- rchisq(cols,B-1);

  # these should be different
  set.seed(runif(1,0,1000000));
  W1 <- rchisq(cols,L-1);
  W3 <- rchisq(cols,B*L*R-B-L+1);

  thisMnum <-  ((B-1)*SBSQS[myr])/(L*R*W2) - ((B*L*R-B-L+1)*SESQS[myr])/(L*R*W3);
  belowzeros <- which(thisMnum<0);
  thisMnum[belowzeros] <- 0;
  thisMden <- ((L-1)*SLSQS[myr])/(B*R*W1) + ((B-1)*SBSQS[myr])/(L*R*W2) + ((B*L*R-B-L)*(B*L*R-B-L+1)*SESQS[myr])/(B*L*R*W3);
  Mmat[myr,] <- thisMnum/thisMden;

 }

 Lowers <- rep(NA,rows);
 Uppers <- rep(NA,rows);
 for (i in 1:rows) {
  Lowers[i] <- quantile(as.numeric(Mmat[i,]),alpha/2);
  Uppers[i] <- quantile(as.numeric(Mmat[i,]),1-alpha/2);
 }
 
 myret <- rep(NA,2);
 myret[1] <- mean(Uppers-Lowers);
 myret[2] <- sqrt(var(Uppers-Lowers)/rows);
 myret;
}



myr0volfun <- function(B,L,R,MCreps,randseed) {
 
 cols <- MCreps;
 lambda1sq <- sigmaesq+B*R*sigmalsq;
 lambda2sq <- sigmaesq+L*R*sigmabsq;
 lambda3sq <- sigmaesq;
 df1 <- L-1;
 df2 <- B-1;
 df3 <- B*L*R-B-L+1;

 # These should always be the same (?)
 set.seed(randseed);
 SBSQS <- lambda2sq * rchisq(rows,df=B-1) / (B-1);
 SLSQS <- lambda1sq * rchisq(rows,df=L-1) / (L-1);

 # These should be different to enhance reproducibility
 set.seed(runif(1,0,100000))
 SESQS <- lambda3sq * rchisq(rows,df=B*L*R-B-L+1)/ (B*L*R-B-L+1);

 Mmat <- matrix(rep(NA,rows*cols),ncol=cols);

 for (myr in 1:rows) {

  # These should also be the same, but not same as above
  set.seed(randseed+myr);
  W1 <- rchisq(cols,L-1);
  W2 <- rchisq(cols,B-1);

  # these should be different
  set.seed(runif(1,0,1000000));
  W3 <- rchisq(cols,B*L*R-B-L+1);

  thisMnum <-  ((B-1)*SBSQS[myr])/(L*R*W2) - ((B*L*R-B-L+1)*SESQS[myr])/(L*R*W3);
  belowzeros <- which(thisMnum<0);
  thisMnum[belowzeros] <- 0;
  thisMden <- ((L-1)*SLSQS[myr])/(B*R*W1) + ((B-1)*SBSQS[myr])/(L*R*W2) + ((B*L*R-B-L)*(B*L*R-B-L+1)*SESQS[myr])/(B*L*R*W3);
  Mmat[myr,] <- thisMnum/thisMden;

 }

 Lowers <- rep(NA,rows);
 Uppers <- rep(NA,rows);
 for (i in 1:rows) {
  Lowers[i] <- quantile(as.numeric(Mmat[i,]),alpha/2);
  Uppers[i] <- quantile(as.numeric(Mmat[i,]),1-alpha/2);
 }
 
 myret <- rep(NA,2);
 myret[1] <- mean(Uppers-Lowers);
 myret[2] <- sqrt(var(Uppers-Lowers)/rows);
 myret;
}





Bisectforbiol <- function(myB,myL,myR,MCreps,randseed) {

 myret <- myb0volfun(myB,myL,myR,MCreps,randseed);
 myret;

} 



Bisectforlab <- function(myB,myL,myR,MCreps,randseed) {

 myret <- myl0volfun(myB,myL,myR,MCreps,randseed);
 myret;

} 

Bisectfortechrep <- function(myB,myL,myR,MCreps,randseed) {

 myret <- myr0volfun(myB,myL,myR,MCreps,randseed);
 myret;

} 




mybisection <- function(f,lower,upper,tol=1e-4) {
# THIS FUNCTION RETURNS A TWO-ELEMENT VECTOR: 
#  ELEMENT 1 IS THE ESTIMATED VALUE AT WHICH THE FUNCTION IS ZERO
#  ELEMENT 2 IS THE NUMBER OF EVALUATIONS TO GET TO THE VALUE (SEARCH STEPS)
# For an arbitrary function "f", this function applies a bisection algorithm
# to find a zero between "lower" and "upper", assuming they are of different signs.

 flow <- f(lower);
 fupper <- f(upper);
 diff <- upper - lower;
 feval <- 2;

 if (flow*fupper>0) {stop("Interval does not contain zero.\n"); }

 while ( abs(diff)>1 & feval < 100) {
  newpoint <- round( (lower+upper)/2 );
  newf <- f(newpoint);
  if (abs(newf)<= tol) break;
  if (newpoint==lower) break;
  if (newpoint==upper) break;
  if (flow*newf < 0) {upper <- newpoint; }
  if (fupper*newf < 0) {lower <- newpoint; }
  diff <- upper-lower;
  feval <- feval+1; 
 } 
 c(newpoint,feval);
}


minWidthB0 <- function(alpha,targetw,L,R,sigmabsq,sigmalsq,sigmaesq,MCReps = 1000000) {

  W1 <- rchisq(MCReps,L-1);
  W2 <- rchisq(MCReps,L-1);
  Gtb0inf <- sigmabsq/(sigmabsq+sigmaesq+sigmalsq * (W2/(L-1))/(W1/(L-1)));
  Right <- quantile(Gtb0inf,1-alpha/2);
  Left <- quantile(Gtb0inf,alpha/2);
  minW <- as.numeric(Right)-as.numeric(Left);
  minW;

}

b0estfun <- function(initial=3,stepsize,target,L,R,MCreps) {

 if (target < minwidthb0varies(alpha,sigmabsq,sigmalsq,sigmaesq,L)) 
  { myret <- Inf; }
 else { 
  myb0s <- rep(NA,40);
  for (i in 1:40) {
   myb0s[i] <- b0fun(initial,stepsize,target,L,R,MCreps,randseed=i)[1] 
  }
  myret <- myb0s; # c(mean(myb0s),sqrt(var(myb0s)));
 }
 myret;
}

b0fun <- function(initial=3,stepsize,target,L,R,MCreps,randseed) {
# Initial is the starting value for number of biological replicates (smallest possible)
# Steps is the size of the steps to increase by until obtaining a width below target
# Target is the targeted mean interval width

 # First, check if the targeted width is achievable. 
  
 curB <- initial;
 nexB <- initial;
 myret <- initial;
 nexmu <- Bisectforbiol(nexB,L,R,MCreps,randseed)[1];

 if (Bisectforbiol(initial,L,R,MCreps,randseed)[1]<target) {
  return(initial); break;
 }

 # Create the function to which the bisection algorithm will be applied
 tempfun <- function(x) Bisectforbiol(x,L,R,MCreps,randseed)[1]-target;
 numsteps <- 0;

 while (tempfun(nexB)>0 & numsteps < 30) {
  nexB <- nexB + stepsize;
  numsteps <- numsteps+1;
 }
 if (numsteps==30) { myret <- c(Inf,Inf); }
 else {
  curB <- nexB-stepsize; 
  mybisret <- mybisection(tempfun,curB,nexB);
  myret <- c(1 + mybisret[1],mybisret[2]);
 }
 myret;
}

l0estfun <- function(initial=3,stepsize,target,B,R,MCreps) {

 if (target < minwidthl0varies(alpha,sigmabsq,sigmalsq,sigmaesq,B)) 
  { myret <- Inf; }
 else { 
  myl0s <- rep(NA,40);
  for (i in 1:40) {
   myl0s[i] <- l0fun(initial,stepsize,target,B,R,MCreps,randseed=i)[1] 
  }
  myret <- c(mean(myl0s),sqrt(var(myl0s)));
 }
 myret;
}

l0fun <- function(initial=3,stepsize,target,B,R,MCreps,randseed) {
# Initial is the starting value for number of biological replicates (smallest possible)
# Steps is the size of the steps to increase by until obtaining a width below target
# Target is the targeted mean interval width

 curL <- initial;
 nexL <- initial;
 myret <- initial;
 nexmu <- Bisectforlab(B,nexL,R,MCreps,randseed)[1];

 if (Bisectforlab(initial,B,R,MCreps,randseed)[1]<target) {
  return(initial); break;
 }

 # Create the function to which the bisection algorithm will be applied
 tempfun <- function(x) Bisectforlab(x,B,R,MCreps,randseed)[1]-target;
 numsteps <- 0;

 while (tempfun(nexL)>0 & numsteps < 30) {
  nexL <- nexL + stepsize;
  numsteps <- numsteps+1;
 }
 if (numsteps==30) { myret <- c(Inf,Inf); }
 else {
  curB <- nexL-stepsize; 
  mybisret <- mybisection(tempfun,curL,nexL);
  myret <- c(1 + mybisret[1],mybisret[2]);
 }
 myret;
}

r0estfun <- function(initial=1,stepsize,target,B,L,MCreps) {

 if (target < minwidthr0varies(alpha,sigmabsq,sigmalsq,sigmaesq,B,L)) 
  { myret <- Inf; }
 else { 
  myr0s <- rep(NA,40);
  for (i in 1:40) {
   myr0s[i] <- r0fun(initial,stepsize,target,B,L,MCreps,randseed=i)[1] 
  }
  myret <- c(mean(myr0s),sqrt(var(myr0s)));
 }
 myret;
}


r0fun <- function(initial=1,stepsize,target,B,L,MCreps,randseed) {
# Initial is the starting value for number of biological replicates (smallest possible)
# Steps is the size of the steps to increase by until obtaining a width below target
# Target is the targeted mean interval width

 curR <- initial;
 nexR <- initial;
 myret <- initial;
 nexmu <- Bisectfortechrep(B,L,nexR,MCreps,randseed)[1];

 if (Bisectfortechrep(B,L,initial,MCreps,randseed)[1]<target) {
  return(initial); break;
 }

 # Create the function to which the bisection algorithm will be applied
 tempfun <- function(x) Bisectfortechrep(B,L,x,MCreps,randseed)[1]-target;
 numsteps <- 0;

 while (tempfun(nexR)>0 & numsteps < 30) {
  nexR <- nexR + stepsize;
  numsteps <- numsteps+1;
 }
 if (numsteps==30) { myret <- c(Inf,Inf); }
 else {
  curB <- nexR-stepsize; 
  mybisret <- mybisection(tempfun,curR,nexR);
  myret <- c(1 + mybisret[1],mybisret[2]);
 }
 myret;
}


### 
## MINIMUM WIDTH FUNCTIONS
###


oldminwidthb0varies <- function(alpha,sigmabsq,sigmalsq,sigmaesq,L,MCreps = 1000000) {

 myw1s <- rchisq(MCreps,L-1);
 mygs <- sigmabsq/(sigmabsq+sigmaesq+sigmalsq/(myw1s/(L-1)));
 Q3 <- as.numeric(quantile(mygs,1-alpha/2));
 Q1 <- as.numeric(quantile(mygs,alpha/2));
 myret <- Q3-Q1;
 myret;

}


minwidthb0varies <- function(alpha,sigmabsq,sigmalsq,sigmaesq,L,MCreps = 1000000) {

 tablerows <- 500;
 myslsqsOverb0r0 <- (sigmalsq/(L-1))*rchisq(tablerows,L-1);

 mytable <- matrix(rep(NA,tablerows*2),ncol=2);
 for (i in 1:tablerows) {
  myw1s <- rchisq(MCreps,L-1);
  thisrow <-  sigmabsq/(sigmabsq+sigmaesq+myslsqsOverb0r0[i]/(myw1s/(L-1)));
  thisQ3 <- as.numeric(quantile(thisrow,1-alpha/2));
  thisQ1 <- as.numeric(quantile(thisrow,alpha/2));
  mytable[i,] <- c(thisQ1,thisQ3);
 }
 myret <- mean(c(mytable[,2]-mytable[,1]));
 myret;
}


minwidthl0varies <- function(alpha,sigmabsq,sigmalsq,sigmaesq,B,MCreps = 1000000) {

 tablerows <- 500;
 mysbsqsOverl0r0 <- (sigmabsq/(B-1))*rchisq(tablerows,B-1);

 mytable <- matrix(rep(NA,tablerows*2),ncol=2);
 for (i in 1:tablerows) {
  myw2s <- rchisq(MCreps,B-1);
  thisrow <-  (mysbsqsOverl0r0[i]/(myw2s/(B-1)))/(sigmalsq+sigmaesq+mysbsqsOverl0r0[i]/(myw2s/(B-1)));
  thisQ3 <- as.numeric(quantile(thisrow,1-alpha/2));
  thisQ1 <- as.numeric(quantile(thisrow,alpha/2));
  mytable[i,] <- c(thisQ1,thisQ3);
 }
 myret <- mean(c(mytable[,2]-mytable[,1]));
 myret;
}

minwidthr0varies <- function(alpha,sigmabsq,sigmalsq,sigmaesq,B,L,MCreps = 1000000) {

 tablerows <- 500;
 myslsqsOverb0r0 <- (sigmalsq/(L-1))*rchisq(tablerows,L-1);
 mysbsqsOverl0r0 <- (sigmabsq/(B-1))*rchisq(tablerows,B-1);

 mytable <- matrix(rep(NA,tablerows*2),ncol=2);
 for (i in 1:tablerows) {
  myw1s <- rchisq(MCreps,L-1);
  myw2s <- rchisq(MCreps,B-1);
  thisrow <-  
  (mysbsqsOverl0r0[i]/(myw2s/(B-1)))/
   (sigmaesq+ myslsqsOverb0r0[i]/(myw1s/(L-1)) + mysbsqsOverl0r0[i]/(myw2s/(B-1)));
  thisQ3 <- as.numeric(quantile(thisrow,1-alpha/2));
  thisQ1 <- as.numeric(quantile(thisrow,alpha/2));
  mytable[i,] <- c(thisQ1,thisQ3);
 }
 myret <- mean(c(mytable[,2]-mytable[,1]));
 myret;
}


# Rao Blackwellization functions

mypigamma <- function(myinput,alpha,beta) {
 # This fixes a problem in pigamma that it returns errors for 0 and smaller
 myret = rep(NA,length(myinput));
 
 for (i in 1:length(myinput)) {
  if (myinput[i] <= 0) { myret[i] = 0; }
  if (myinput[i] > 0) { 
   myret[i] = pigamma(myinput[i], alpha,beta);
  }
 }
 myret;
}


mycircrb <- function(g0,M1,M2,M3,a,b,c,d) {
 # This calculates an estimate of the cdf for g given MC W's and one set of S's
 for (i in 1:length(M1)) { 
  input23[i] <- (max(M2[i]*a-M3[i]*b,0)- M2[i]*a*g0 - M3[i]*d*g0)/(c*g0) ;
 }
 
 f23 <- 
  1-mypigamma( input23
  ,alpha=(l0-1)/2,beta=1/2);
  myret <- f23;
  myret <- mean(myret);
  myret;
}

quantilecircrb <- function(myquant,M1,M2,M3,a,b,c,d) {
 # This calculates the quantile function (inverse cdf) given MC W's and one set of S's

 myfun <- function(x) { mycircrb(x,M1,M2,M3,a,b,c,d)-myquant; }
 rootfind <- uniroot(myfun,lower=0.0000001,upper=0.999999);
 myres <- rootfind$root;
 myres;

}






