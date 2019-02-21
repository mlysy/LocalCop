// Reimplementation of PDF in BiCopPDF

#include <Rcpp.h>
using namespace Rcpp;
using namespace R;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define pi 3.14159265358979

double StableGammaDivision(double x1, double x2)
{
  int i;
  double a1, a2, b1, b2, sum=1.0;
  a1 = fmod(MAX(x1,x2),1.0);
  a2 = MAX(x1,x2)-a1;
  b1 = fmod(MIN(x1,x2),1.0);
  b2 = MIN(x1,x2)-b1;
  if(a1==0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=b2 ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
  }
  else if(a1>0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=(int)b2 ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1);
  }
  else if(a1==0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum /= gammafn(b1);
  }
  else if(a1>0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1)/gammafn(b1);
  }
  if(x2 > x1) sum = 1.0/sum;
  return sum;
}

// [[Rcpp::export]]
NumericVector StudentLL(NumericVector u, NumericVector v, NumericVector rho, double nu) {
  int n = u.size();
  double t1, t2, f;
  NumericVector ll(n);
  
  for(int j=0;j<n;j++)
    {
      // dat[0] = u[j]; dat[1] = v[j];
      t1 = Rf_qt(u[j],nu,1,0); t2 = Rf_qt(v[j],nu,1,0);
      f = StableGammaDivision((nu+2.0)/2.0,nu/2.0)/(nu*pi*sqrt(1.0-pow(rho[j],2.0))*Rf_dt(t1,nu,0)*Rf_dt(t2,nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho[j]*t1*t2)/(nu*(1.0-pow(rho[j],2.0))),-(nu+2.0)/2.0);
      ll[j] = f;
      // if(log(f)>XINFMAX) ll += log(XINFMAX);
      // else if(f < DBL_MIN) ll += log(DBL_MIN);
      // else ll += log(f);
    }
  return ll;
}
