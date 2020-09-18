#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

//' Calculates it-kendall-tau
//' 
//' @param x numeric vector
//' @param y numeric vector
//' @param perspective should we consider the "local" or "global" perspective?
//' 
//' @importFrom Rcpp sourceCpp
//' @name kendallc
//' @export
//' @useDynLib visualizationQualityControl
//' @return kendall tau correlation
// [[Rcpp::export]]
double kendallt(NumericVector x, NumericVector y, String perspective = "local") {
  
  double sum_concordant = 0;
  double sum_discordant = 0;
  double sum_x_ties = 0;
  double sum_y_ties = 0;
  double sum_tied_x = 0;
  double sum_tied_y = 0;
  double sum_tied_x_na = 0;
  double sum_tied_y_na = 0;
  double sum_all_na = 0;
  bool is_concordant;
  bool is_discordant;
  bool is_tied_x;
  bool is_tied_y;
  bool is_tied_x_na;
  bool is_tied_y_na;
  bool is_all_na;
  double k_numerator;
  double k_denominator;
  double k_tau;
  
  LogicalVector matching_na;
  double n_matching_na;
  
  if (perspective == "local") {
    matching_na = is_na(x) & is_na(y);
    //n_matching_na = sum(matching_na);
    x = x[!matching_na];
    y = y[!matching_na];
  }
  
  double n_entry = x.size();
  
  if (n_entry < 2) {
    return NA_REAL;
  }
  
  for (double i = 0; i < (n_entry - 1); i++) {
    for (double j = (i+1); j < n_entry; j++) {
      if (perspective == "global") {
        is_concordant = (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||     //3 1
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||   // ## 2
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & !NumericVector::is_na(y[i]) & NumericVector::is_na(y[j])) ||                             // ## 3
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j])) ||                            //  ## 4
          (!NumericVector::is_na(x[i]) & NumericVector::is_na(x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||                            //  ## 5
          (NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||                            //  ## 6
          (!NumericVector::is_na(x[i]) & NumericVector::is_na(x[j]) & !NumericVector::is_na(y[i]) & NumericVector::is_na(y[j])) ||                                             //            ## 7
          (NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]));
        
        is_discordant = (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & !NumericVector::is_na(y[i]) & NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) & NumericVector::is_na(x[j])  & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||
          (NumericVector::is_na(x[i])  & !NumericVector::is_na(x[j])  & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||
          (!NumericVector::is_na(x[i]) & NumericVector::is_na(x[j])  & NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j])) ||
          (NumericVector::is_na(x[i])  & !NumericVector::is_na(x[j])  & !NumericVector::is_na(y[i]) & NumericVector::is_na(y[j]));
      } else {
        is_concordant = 
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||    // #3 1
        (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||  // ## 2
        (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & !NumericVector::is_na(y[i]) & NumericVector::is_na(y[j])) ||                             // ## 3
        (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j])) ||                             // ## 4
        (!NumericVector::is_na(x[i]) & NumericVector::is_na(x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||                             // ## 5
        (NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])); 
        
        is_discordant = 
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j])) ||
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] > x[j]) & NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] < x[j]) & !NumericVector::is_na(y[i]) & NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) & NumericVector::is_na(x[j])  & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] < y[j])) ||
          (NumericVector::is_na(x[i])  & !NumericVector::is_na(x[j])  & !NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] > y[j]));
      }
      
      is_tied_x = (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] == x[j]) & (!NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] != y[j])));
      is_tied_y = (!NumericVector::is_na(x[i]) & !NumericVector::is_na(x[j]) & (x[i] != x[j]) & (!NumericVector::is_na(y[i]) & !NumericVector::is_na(y[j]) & (y[i] == y[j])));
      is_tied_x_na = NumericVector::is_na(x[i]) & NumericVector::is_na(x[j]) & (!NumericVector::is_na(y[i]) | !NumericVector::is_na(y[j]));
      is_tied_y_na = (!NumericVector::is_na(x[i]) | !NumericVector::is_na(x[j])) & NumericVector::is_na(y[i]) & NumericVector::is_na(y[j]);
      is_all_na = NumericVector::is_na(x[i]) & NumericVector::is_na(x[j]) & NumericVector::is_na(y[i]) & NumericVector::is_na(y[j]);

      if (is_concordant) {
        sum_concordant++;
      }
      if (is_discordant) {
        sum_discordant++;
      }
      if (is_tied_x) {
        sum_tied_x++;
      }
      if (is_tied_y) {
        sum_tied_y++;
      }
      if (is_tied_x_na) {
        sum_tied_x_na++;
      }
      if (is_tied_y_na) {
        sum_tied_y_na++;
      }
      if (is_all_na) {
        sum_all_na++;
      }
    }
  }
  
  double half_sum_na_ties = sum_all_na / 2;
  
  if (perspective == "global") {
    sum_x_ties = sum_tied_x + sum_tied_x_na + half_sum_na_ties;
    sum_y_ties = sum_tied_y + sum_tied_y_na + half_sum_na_ties;
  } else {
    sum_x_ties = sum_tied_x;
    sum_y_ties = sum_tied_y;
  }
  
  k_numerator = sum_concordant - sum_discordant;
  k_denominator = sum_discordant + sum_concordant + sum_x_ties + sum_y_ties;
  
  
  // debugging
  // Rprintf("n_entry: %f \n", n_entry);
  // Rprintf("x_ties: %f \n", sum_tied_x);
  // Rprintf("x_na_ties: %f \n", sum_tied_x_na);
  // Rprintf("sum_x_ties: %f \n", sum_x_ties);
  // Rprintf("sum_y_ties: %f \n", sum_y_ties);
  // Rprintf("sum_na_ties: %f \n", sum_all_na);
  // Rprintf("half_sum_na_ties: %f \n", half_sum_na_ties);
  // Rprintf("sum_concordant: %f \n", sum_concordant);
  // Rprintf("sum_discordant: %f \n", sum_discordant);
  // Rprintf("k_numerator: %f \n", k_numerator);
  // Rprintf("k_denominator: %f \n", k_denominator);
  // 
  k_tau = k_numerator / k_denominator;
  
  return k_tau;
}

