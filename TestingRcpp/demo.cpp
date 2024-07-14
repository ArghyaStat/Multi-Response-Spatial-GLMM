#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


NumericVector ar1Cpp(int n, double rho, double start){
  NumericVector out(n);
  NumericVector foo(n);
  foo = rnorm(n, 0, 1);
  out[0] = start;
  
  for(int i = 1; i < n; i++)
  {
    out[i] = rho*out[i-1] + foo[i];
  }
  return out;
}

// 
// int candles(int age){
//   int remain = age;
//   int att = 0;
//   IntegerVector samp;  
//   IntegerVector temp(1);
//   while(remain > 0){
//     samp = seq(1, remain);
//     temp = sample(samp, 1); // sample returns a vector
//     remain = remain - temp[0]; //picking up the first element
//     att++;
//   }
//   return att;
// }
// 
// 
// IntegerVector repCandlesC(int age, int n) {
//   IntegerVector reps(n);
//   
//   for(int i = 0; i < n; i++)
//   {
//     reps[i] = candles(age);
//   }
//   return reps;
// }
// 
// 
