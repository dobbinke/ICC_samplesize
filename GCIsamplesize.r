mcrunsMS <- 100;
mcrunsW = 100;
alpha = 0.05;


####
##  Beginning of code cut-and-pasted from RB file
####
install.packages("pscl");


library(pscl);


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


mycircrb <- function(g0,M1,M2,M3,a,b,c,d,l0) {

 input23 = rep(NA,length(M1));
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

quantilecircrb <- function(myquant,M1,M2,M3,a,b,c,d,l0) {

 myfun <- function(x) { mycircrb(x,M1,M2,M3,a,b,c,d,l0)-myquant; }
 if (myfun(0.0000001)>= 0) { myres = 0; }
 else {
  rootfind <- uniroot(myfun,lower=0.0000001,upper=0.999999);
  myres <- rootfind$root;
 }
 myres;

}


####
##  Ending of code cut-and-pasted from RB file
####


####
##  Beginning of code cut-and-pasted from simplemcwithrbvscvwithrb.txt file
####


CtrlVoverMS = function(sigmabsq,sigmalsq,sigmaesq,b0,l0,r0,randomseed=2) {

  set.seed(randomseed);

  iccb = sigmabsq/10;

  lambda1sq <- sigmaesq+b0*r0*sigmalsq;
  lambda2sq <- sigmaesq+l0*r0*sigmabsq;

 #Generate ws, using same set of w's for each mean square

 sbsq = (sigmaesq+sigmabsq*l0*r0)*rchisq(mcrunsMS,b0-1)/(b0-1);
 slsq = (sigmaesq+sigmalsq*b0*r0)*rchisq(mcrunsMS,l0-1)/(l0-1);
 sesq = sigmaesq*rchisq(mcrunsMS,b0*l0*r0-b0-l0+1)/(b0*l0*r0-b0-l0+1);
 a = (b0-1)*sbsq/(l0*r0);
 b = (b0*l0*r0-b0-l0+1)*sesq/(l0*r0);
 c = (l0-1)*slsq/(b0*r0);
 d = b*(b0*l0*r0-b0-l0)/b0;
 theselowers = rep(NA,mcrunsMS);
 theseuppers = rep(NA,mcrunsMS);

 retwidths = rep(NA,mcrunsMS);
 # Generate mean squares and CI's
 for (i in 1:mcrunsMS) {
  W1 <- rchisq(mcrunsW,l0-1);
  W2 <- rchisq(mcrunsW,b0-1);
  W3 <- rchisq(mcrunsW,b0*l0*r0-b0-l0+1);
  M1 = 1/W1;
  M2 = 1/W2;
  M3 = 1/W3;
  theselowers[i] = quantilecircrb(0.025,M1,M2,M3,a[i],b[i],c[i],d[i],l0);
  theseuppers[i] = quantilecircrb(0.975,M1,M2,M3,a[i],b[i],c[i],d[i],l0);
 }
 myrats = slsq/sbsq;
 myFvar = (lambda1sq/lambda2sq)^2*2*(b0-1)^2*(b0-1+l0-1-2)/((l0-1)*(b0-1-2)^2*(b0-1-4));
 mycovLower = cov(theselowers,myrats);
 CtrlVLower = theselowers - (mycovLower/myFvar)*((slsq/sbsq) - (lambda1sq/lambda2sq)*(b0-1)/(b0-1-2));
 mycovUpper = cov(theseuppers,myrats);
 CtrlVUpper = theseuppers - (mycovUpper/myFvar)*((slsq/sbsq) - (lambda1sq/lambda2sq)*(b0-1)/(b0-1-2));
 retwidths = CtrlVUpper - CtrlVLower;

 for (i in 1:mcrunsMS) { retwidths[i] = max(0,retwidths[i]); }

 retwidths;
}


####
##  Ending of code cut-and-pasted from simplemcwithrbvscvwithrb.txt file
####





myGCIwidthfun = function(sigmabsq,sigmalsq,sigmaesq,b0,l0,r0,randomseed=2) {

  widthvec = CtrlVoverMS(sigmabsq,sigmalsq,sigmaesq,b0,l0,r0);

  myret = mean(widthvec);

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

b0estfun<- function(initial=3,biggest=1000,target,sigmabsq,sigmalsq,sigmaesq,l0,r0,randomseed=2) {

  myinitial = initial;

  tempfun = function(x) myGCIwidthfun(sigmabsq,sigmalsq,sigmaesq,x,l0,r0,randomseed) - target;

  if (tempfun(myinitial) > 0) {
   myres = mybisection(tempfun,myinitial,biggest)
  }
  else {
   myres = myinitial;
  }
  myres;

}

samp.size.ests = rep(NA,5);
starttime = Sys.time();
for (i in 1:5) {
 samp.size.ests[i] = b0estfun(initial=30,biggest=500,target=.10,sigmabsq=0.85888,sigmalsq=0.0135125,sigmaesq=0.0722,10,5,randomseed=i)[1];
}
endtime = Sys.time();
timeper = (endtime-starttime)/2;
timeper

#myoutput = c(mean(samp.size.ests),sqrt(var(samp.size.ests)),sigmabsq,sigmalsq,sigmaesq);

samp.size.ests;


#write(file="GCImethodresults.txt",myoutput,append=T)
