#include <Rcpp.h>
#include <string>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector gcContent(std::string x, int windowSize) {
  int size = 0;
  if(x.size() % windowSize == 0)
    size = x.size() / windowSize;
  else
    size = std::floor(x.size() / windowSize) + 1;
  NumericVector gcCount(size);
  for(unsigned int i = 0; i < x.size(); i += windowSize){
    int count = 0;
    for(int j = 0; j < windowSize; j++){
      if(x[i+j] == 'G' || x[i+j] == 'C')
        count++;
    }
    gcCount[i/windowSize] = count;
  }
  return gcCount;
}
