#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>


#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

using namespace Rcpp;
using namespace sugar;


// [[Rcpp::export]]
IntegerVector oneMultinomCall(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    return(ans);
}

// [[Rcpp::export]]
int whichMultinom(NumericVector probs) {
  //    int k = probs.size();
    IntegerVector v = oneMultinomCall(probs);
    int ret = which_max(v);
    //    Rprintf("ret %i \n",ret);
    return(ret);
}

// [[Rcpp::export]]
NumericMatrix submat(NumericMatrix X, LogicalVector condition) { 
    int n=X.nrow(), k=X.ncol();
    NumericMatrix out(sum(condition),k);
    for (int i = 0, j = 0; i < n; i++) {
        if(condition[i]) {
            out(j,_) = X(i,_);
            j = j+1;
        }
    }
    return(out);
}

// [[Rcpp::export]]
IntegerVector getsrcC(NumericMatrix tmat, NumericVector nv) { 
  int i;
  int k=tmat.nrow();
  int choice;
  NumericVector mprb(k) ;
  double tmp=0;
  NumericVector prbs(k);
  IntegerVector out(k);

  for (i=0;i<k;i++)
    {
      prbs=tmat(i,_);
      prbs(i)=0;
      
      mprb=prbs * nv;

      if (sum(mprb)>1) //prob of dispersing out of cell needt to be normalized; should be really uncommon
	{
	  mprb=mprb/sum(mprb);
	  mprb[i]=0;
	} else { //more normal process
	mprb[i]=1-sum(mprb);
	mprb=mprb/sum(mprb);
      }
      
      //      Rprintf("i %i, out[i], %i, mprb[i], %g\n",i,out[i],mprb[i]);
      
      choice = whichMultinom(mprb);
      
      //      Rprintf("out[i] %i \n",out[i]);
      
      out[i]=choice+1;
      
      //      Rprintf("out[i] %i \n",out[i]);
      
    }
  
  return(out);
}

//not used currently
// [[Rcpp::export]]
IntegerVector getsrcC2(NumericMatrix tmat, NumericVector nv) { 
  int i;
  int k=tmat.nrow();
  int choice;
  double home =0.0;
  NumericVector mprb(k) ;
  NumericVector prbs(k);
  IntegerVector out(k);

  for (i=0;i<k;i++)
    {
      if (nv[i]==0) //only work on rows that could be colonized (popsize==0)
	{
	  prbs=tmat(i,_);
	  home=(prbs/sum(prbs))[i];
	  mprb=prbs * nv;

	  if (sum(mprb)>0)
	    {
	      mprb(i)= sum(mprb)*(home/(1-home));
	      mprb=mprb/sum(mprb);
	      choice = whichMultinom(mprb);
	      //choice = gsl_mult1(mprb)[0];
	      out[i]=choice+1;
	    } else {out[i]=i+1;}  //there are no source populations
	} else //this population has been colonized
	{
	  out[i]=i+1;
	}
    }
  
  return(out);
}

// [[Rcpp::export]]
NumericMatrix makeRmatC(NumericMatrix pops, NumericMatrix dm, NumericVector cent)
{
  int npop=pops.nrow();
  int offy;
  int offx;
  int dmrow=dm.nrow();
  int dmcol=dm.ncol();
  NumericMatrix rmat(npop,npop);
  int i=0;
  int j=0;
  Rprintf("not failed yet");
  for (i=0;i<npop;i++)    //rows
    for (j=0;j<npop;j++)  //cols
      {
	offx=pops(i,2)-pops(j,2);
	offy=pops(i,1)-pops(j,1);
	/*
	if (((offx+cent[1])>0) & ((offx+cent[1])<dmcol))
	  if (((offy+cent[0])>0) & ((offy+cent[0])<dmrow))
	    rmat(i,j)=dm(offx+cent[1],offy+cent[0]);
	*/

	if ((abs(offx)+cent[1])<dmrow) 
	  if ((abs(offy)+cent[0])<dmcol)
	    rmat(i,j)=dm(offx+cent[1],offy+cent[0]);
      }
  return(rmat);
}
