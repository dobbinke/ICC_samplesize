

# This source file has the mybisection function
#source("C:\\Documents and Settings\\dobbinke\\My Documents\\My Dropbox\\CurrentFiles\\JournalsInProgress\\SampSizeGCI\\Simulations\\Rprog.Core.SampSizes.txt")

source("/Users/kevindobbin/Desktop/git/ICC_samplesize/Rprog.Core.SampSizes.txt")
 


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





myMLSvolfun <- function(alpha,b0,l0,r0,sigmabsq,sigmalsq,sigmaesq,MCruns=100000,randseed=5) {

 set.seed(randseed);
 bigreps <- MCruns;
 bigdf <- 10000000;
 
 # Calculate vectors of observed mean squares
 mysbsqs <- ((sigmaesq+l0*r0*sigmabsq)/(b0-1))*rchisq(bigreps,b0-1);
 myslsqs <- ((sigmaesq+b0*r0*sigmalsq)/(l0-1))*rchisq(bigreps,l0-1);
 mysesqs <- (sigmaesq/(b0*l0*r0-b0-l0+1))*rchisq(bigreps,b0*l0*r0-b0-l0+1);

 # Constants from  Burdick book
 myG2 <- 1-qf(alpha/2,bigdf,b0-1);
 myF5 <- qf(1-alpha/2,b0-1,b0*l0*r0-b0-l0+1);
 myF4 <- qf(alpha/2,l0-1,b0-1);
 myH2 <- qf(1-alpha/2,bigdf,b0-1)-1;
 myF6 <- qf(alpha/2,b0-1,b0*l0*r0-b0-l0+1);
 myF3 <- qf(1-alpha/2,l0-1,b0-1);

 # Calculate the lower bounds for the bigreps intervals
 Lstarnums <- 
 b0*(1-myG2)*mysbsqs^2-b0*mysbsqs*mysesqs+b0*(myF5-(1-myG2)*myF5^2)*mysesqs^2;

 Lstardenoms <- 
 l0*(b0*r0-1)*mysbsqs*mysesqs+l0*(1-myG2)*mysbsqs*myslsqs/myF4;

 Lstars <- Lstarnums/Lstardenoms;
 for (i in 1:length(Lstars)) { Lstars[i] <- max(Lstars[i],0); }

 Ls <- Lstars/(1+Lstars);  

 # Calculate the upper bounds for the bigreps interval
 Ustarnums <- 
 b0*(1+myH2)*mysbsqs^2-b0*mysbsqs*mysesqs+b0*(myF6-(1+myH2)*myF6^2)*mysesqs^2;

 Ustardenoms <- 
 l0*(b0*r0-1)*mysbsqs*mysesqs+l0*(1+myH2)*mysbsqs*myslsqs/myF3;

 Ustars <- Ustarnums/Ustardenoms;
 for (i in 1:length(Ustars)) { Ustars[i] <- max(Ustars[i],0); }

 Us <- Ustars/(1+Ustars);
 # Returns the mean of the interval widths
 meandiffs <- mean(Us-Ls);

 myret <- meandiffs;
 myret;
}

now <- Sys.time();
temp <- myMLSvolfun(alpha=0.05,b0=100,l0=100,r0=100,sigmabsq=9,sigmalsq=0.5,sigmaesq=0.5,MCruns=100000,randseed=5); 
later <- Sys.time();

later-now;


myMLSvolfunRBization <- function(alpha,b0,l0,r0,sigmabsq,sigmalsq,sigmaesq,MCruns=100000,randseed=5) {

 set.seed(randseed);
 bigreps <- MCruns;
 bigdf <- 10000000;
 
 # Uses Rao-Blackwellization with RBreps on sigmalsq
 RBreps <- 2;
 myslsqs <- ((sigmaesq+b0*r0*sigmalsq)/(l0-1))*rchisq(RBreps,l0-1);


 # Calculate vectors of observed mean squares
 mysbsqs <- ((sigmaesq+l0*r0*sigmabsq)/(b0-1))*rchisq(bigreps,b0-1);
 mysesqs <- (sigmaesq/(b0*l0*r0-b0-l0+1))*rchisq(bigreps,b0*l0*r0-b0-l0+1);

 # Constants from  Burdick book
 myG2 <- 1-qf(alpha/2,bigdf,b0-1);
 myF5 <- qf(1-alpha/2,b0-1,b0*l0*r0-b0-l0+1);
 myF4 <- qf(alpha/2,l0-1,b0-1);
 myH2 <- qf(1-alpha/2,bigdf,b0-1)-1;
 myF6 <- qf(alpha/2,b0-1,b0*l0*r0-b0-l0+1);
 myF3 <- qf(1-alpha/2,l0-1,b0-1);

 mymeandiffs <- rep(NA,RBreps);
 for (RBstep in 1:RBreps) {
  # Calculate the lower bounds for the bigreps intervals
  Lstarnums <- 
  b0*(1-myG2)*mysbsqs^2-b0*mysbsqs*mysesqs+b0*(myF5-(1-myG2)*myF5^2)*mysesqs^2;
 
  Lstardenoms <- 
  l0*(b0*r0-1)*mysbsqs*mysesqs+l0*(1-myG2)*mysbsqs*myslsqs[RBstep]/myF4;
 
  Lstars <- Lstarnums/Lstardenoms;
  for (i in 1:length(Lstars)) { Lstars[i] <- max(Lstars[i],0); }
 
  Ls <- Lstars/(1+Lstars);  
 
  # Calculate the upper bounds for the bigreps interval
  Ustarnums <- 
  b0*(1+myH2)*mysbsqs^2-b0*mysbsqs*mysesqs+b0*(myF6-(1+myH2)*myF6^2)*mysesqs^2;
 
  Ustardenoms <- 
  l0*(b0*r0-1)*mysbsqs*mysesqs+l0*(1+myH2)*mysbsqs*myslsqs[RBstep]/myF3;
 
  Ustars <- Ustarnums/Ustardenoms;
  for (i in 1:length(Ustars)) { Ustars[i] <- max(Ustars[i],0); }
 
  Us <- Ustars/(1+Ustars);
  # Returns the mean of the interval widths
  mymeandiffs[RBstep] <- mean(Us-Ls);
 }
 myret <- mean(mymeandiffs);
 myret;
}


now <- Sys.time();
temp <- myMLSvolfunRBization(alpha=0.05,b0=100,l0=100,r0=100,sigmabsq=9,sigmalsq=0.5,sigmaesq=0.5,MCruns=100000,randseed=5); 
later <- Sys.time();

later-now;
# Rao-Blackwellization actually seems to make the program run slower
# It seems to be that the vector calculations are quick, but the loop repeats are slow



MLSb0fun <- function(alpha=0.05,initial=3,stepsize,target,L,R,sigmabsq,sigmalsq,sigmaesq,MCreps=100000,randseed=5) {
# Initial is the starting value for number of biological replicates (smallest possible)
# Steps is the size of the steps to increase by until obtaining a width below target
# Target is the targeted mean interval width

 # First, check if the targeted width is achievable. 
  
 curB <- initial;
 nexB <- initial;
 myret <- initial;
 nexmu <- myMLSvolfun(alpha,nexB,L,R,sigmabsq,sigmalsq,sigmaesq,MCreps,randseed)[1];

 if (myMLSvolfun(alpha,curB,L,R,sigmabsq,sigmalsq,sigmaesq)[1]<target) {
  return(initial); break;
 }

 # Create the function to which the bisection algorithm will be applied
 tempfun <- function(x) myMLSvolfun(alpha,x,L,R,sigmabsq,sigmalsq,sigmaesq,MCreps,randseed)[1]-target;
 numsteps <- 0;

 while (tempfun(nexB)>0 & numsteps < 30) {
  nexB <- nexB + stepsize;
  numsteps <- numsteps+1;
 }
 if (numsteps==30) { myret <- c(Inf,Inf); }
 else {
  curB <- nexB-stepsize; 
  mybisret <- mybisection(tempfun,curB,nexB);
  myret <- c(1 + mybisret[1]);
 }
 myret;
}




#MLSb0fun(alpha=0.05,initial=3,stepsize=30,target=0.5,L=5,R=1,sigmabsq=9.9,sigmalsq=0.05,sigmaesq=0.05) 


MLSb0estfun <- function(alpha,initial,stepsize,target,L,R,sigmabsq,sigmalsq,sigmaesq,MCreps) {

# if (target < minwidthb0varies(alpha,sigmabsq,sigmalsq,sigmaesq,L)) 
#  { myret <- Inf; }
# else { 
  myb0s <- rep(NA,40);
  for (i in 1:40) {
   myb0s[i] <- MLSb0fun(alpha,initial,stepsize,target,L,R,sigmabsq,sigmalsq,sigmaesq,MCreps,randseed=i)[1] 
  }
  myret <- c(mean(myb0s),sqrt(var(myb0s)/length(myb0s)),min(myb0s),max(myb0s));
# }
 myret;
}

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.5,L=5,R=1,sigmabsq=9.9,sigmalsq=0.05,sigmaesq=0.05,MCreps=1000) 

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.4,L=5,R=1,sigmabsq=9.9,sigmalsq=0.05,sigmaesq=0.05,MCreps=1000) 

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.3,L=5,R=1,sigmabsq=9.9,sigmalsq=0.05,sigmaesq=0.05,MCreps=1000) 

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.2,L=5,R=1,sigmabsq=9.9,sigmalsq=0.05,sigmaesq=0.05,MCreps=1000) 

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.1,L=5,R=1,sigmabsq=9.9,sigmalsq=0.05,sigmaesq=0.05,MCreps=1000) 


MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.5,L=5,R=1,sigmabsq=9.,sigmalsq=0.5,sigmaesq=0.5,MCreps=1000) 

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.4,L=5,R=1,sigmabsq=9.,sigmalsq=0.5,sigmaesq=0.5,MCreps=1000) 

MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.3,L=5,R=1,sigmabsq=9.,sigmalsq=0.5,sigmaesq=0.5,MCreps=1000) 


MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.5,L=5,R=1,sigmabsq=7.0,sigmalsq=1.5,sigmaesq=1.5,MCreps=1000) 


MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.5,L=5,R=1,sigmabsq=3.0,sigmalsq=3.5,sigmaesq=3.5,MCreps=1000) 
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.4,L=5,R=1,sigmabsq=3.0,sigmalsq=3.5,sigmaesq=3.5,MCreps=1000) 
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.3,L=5,R=1,sigmabsq=3.0,sigmalsq=3.5,sigmaesq=3.5,MCreps=1000) 


MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.5,L=5,R=1,sigmabsq=0.0,sigmalsq=5,sigmaesq=5,MCreps=1000) 
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.4,L=5,R=1,sigmabsq=0.0,sigmalsq=5,sigmaesq=5,MCreps=1000) 
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.3,L=5,R=1,sigmabsq=0.0,sigmalsq=5,sigmaesq=5,MCreps=1000) 
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.2,L=5,R=1,sigmabsq=0.0,sigmalsq=5,sigmaesq=5,MCreps=1000) 
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.1,L=5,R=1,sigmabsq=0.0,sigmalsq=5,sigmaesq=5,MCreps=1000) 


# focus on just target=0.3

# This is the ICC=0.7 case
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.3,L=10,R=5,sigmabsq=7,sigmalsq=1.5,sigmaesq=1.5,MCreps=1000) 


# This is the ICC=0.5 case
MLSb0estfun(alpha=0.05,initial=3,stepsize=30,target=0.3,L=10,R=5,sigmabsq=5,sigmalsq=2.5,sigmaesq=2.5,MCreps=1000) 




