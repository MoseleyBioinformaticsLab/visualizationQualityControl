#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

//' Returns sorted index in Rcpp
//' 
//' @param x the vector to be sorted
//' 
//' @importFrom Rcpp sourceCpp
//' @export
//' @useDynLib visualizationQualityControl
//' @return ordered indices
// [[Rcpp::export]]
IntegerVector sortedIndex(NumericVector x){
  IntegerVector idx = seq_along(x) - 1;
  
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});
  
  return idx;
}

// [[Rcpp::export]]
IntegerVector compare_self(NumericVector x){
  int n_entry = x.size();
  IntegerVector match_self (n_entry);
  match_self[0] = 1;
  
  int idx = 1;
  
  for (int i = 1; i < (n_entry); i++) {
    if (x[i] != x[(i - 1)]) {
      match_self[idx] = 1;
    } else {
      match_self[idx] = 0;
    }
    idx++;
  }
  return match_self;
}

// [[Rcpp::export]]
IntegerVector compare_both(IntegerVector x, IntegerVector y){
  int n_entry = x.size();
  IntegerVector match_self (n_entry);
  match_self[0] = 1;
  
  int idx = 1;
  
  for (int i = 1; i < (n_entry); i++) {
    if ((x[i] != x[(i - 1)]) | (y[i] != y[(i - 1)])) {
      match_self[idx] = 1;
    } else {
      match_self[idx] = 0;
    }
    idx++;
  }
  match_self.push_back(1);
  return match_self;
}

// [[Rcpp::export]]
IntegerVector which_notzero(IntegerVector x){
  IntegerVector notzero (x.size());
  int idx = 0;
  
  for (int i = 0; i < x.size(); i++) {
    if (x[i] != 0) {
      notzero[idx] = i;
      idx++;
    }
  }
  IntegerVector keep_loc = seq(0, (idx - 1));
  notzero = notzero[keep_loc];
  return notzero;
}

// [[Rcpp::export]]
int kendall_discordant(IntegerVector x, IntegerVector y){
  double sup = 1 + max(y);

  IntegerVector arr(sup, 0);
  double i = 0;
  double k = 0;
  int n = x.size();
  int idx = 0;
  int dis = 0;
  
  while (i < n){
    while ((k < n) && (x[i] == x[k])) {
      dis = dis + i;
      idx = y[k];
      while (idx != 0) {
        dis = dis - arr[idx];
        idx = idx & (idx - 1);
      }
      k++;
    }
    while (i < k) {
      idx = y[i];
      while (idx < sup) {
        arr[idx] = arr[idx] + 1;
        idx = idx + (idx & (-1*idx));
      }
      i++;
    }
  }
  return dis;
}

inline double signC(double x) {
  if (x > 0) {
    return 1.0;
  } else if (x == 0) {
    return 0.0;
  } else {
    return -1.0;
  }
}

//' Calculates ici-kendall-tau
//' 
//' @param x numeric vector
//' @param y numeric vector
//' @param perspective should we consider the "local" or "global" perspective?
//' 
//' @details Calculates the information-content-informed Kendall-tau correlation measure.
//'   This correlation is based on concordant and discordant ranked pairs, like Kendall-tau,
//'   but also includes missing values (as NA). Missing values are assumed to be *primarily* due
//'   to lack of detection due to instrumental sensitivity, and therefore encode *some* information.
//'   
//'   For more details see the ICI-Kendall-tau vignette:
//'   \code{vignette("ici-kendalltau", package = "visualizationQualityControl")}
//' 
//' @examples 
//' data("grp_cor_data")
//' exp_data = grp_cor_data$data
//' x = exp_data[, 1]
//' y = exp_data[, 2]
//' kendallt(x, y)
//' cor(x, y, method = "kendall") 
//' 
//' x = sort(rnorm(100))
//' y = x + 1
//' y2 = y
//' y2[1:10] = NA
//' kendallt(x, y)
//' kendallt(x, y2, "global")
//' kendallt(x, y2)
//' 
//' @importFrom Rcpp sourceCpp
//' @export
//' @useDynLib visualizationQualityControl
//' @return kendall tau correlation
// [[Rcpp::export]]
int ici_kendallt_mergesort(NumericVector x, NumericVector y, String perspective = "local", String alternative = "two.sided", String output = "simple") {
  
  if (x.length() != y.length()) {
    throw std::range_error("X and Y are not the same length!");
    exit(-1);
  }
  
  double sum_concordant = 0;
  double sum_discordant = 0;
  double sum_x_ties = 0;
  double sum_y_ties = 0;
  double sum_tied_x = 0;
  double sum_tied_y = 0;
  double sum_tied_x_na = 0;
  double sum_tied_y_na = 0;
  double sum_all_na = 0;
  double k_numerator;
  double k_denominator;
  double k_tau;
  bool reject_concordant;
  bool reject_discordant;
  
  // for generating the p-value
  //double xties = 0;
  //double yties = 0;
  double t_0 = 0;
  double s_adjusted = 0;
  double x_tied_sum_t1;
  double y_tied_sum_t2;
  double v_0_sum = 0;
  double v_t_sum = 0;
  double v_u_sum = 0;
  double v_t1_sum = 0;
  double v_t2_sum = 0;
  double s_adjusted_variance = 0;
  NumericVector z_b (1);
  NumericVector p_value (1);
  
  LogicalVector matching_na;
  //double n_matching_na;
  
  if (perspective == "local") {
    matching_na = is_na(x) & is_na(y);
    //n_matching_na = sum(matching_na);
    x = x[!matching_na];
    y = y[!matching_na];
  }
  
  NumericVector x2 = clone(x);
  NumericVector y2 = clone(y);
  
  int n_na_x = sum(is_na(x));
  int n_na_y = sum(is_na(y));
  
  if ((n_na_x == x.size()) || (n_na_y == y.size())) {
    return 0.0;
  }
  
  x2 = x[!is_na(x)];
  y2 = y[!is_na(y)];
  
  double min_value = min(NumericVector::create(min(x2), min(y2)));
  double na_value = min_value - 0.1;
  
  x2 = clone(x);
  y2 = clone(y);
  x2[is_na(x)] = na_value;
  y2[is_na(y)] = na_value;
  
  
  int n_entry = x2.size();
  Rprintf("n_entry: %i\n", n_entry);
  
  if (n_entry < 2) {
    return 0.0;
  }
  
  IntegerVector low_subset = seq(1, (n_entry - 1));
  //Rprintf("n_low: %i\n", low_subset.size());
  IntegerVector hi_subset = seq(0, (n_entry - 2));
  //Rprintf("n_hi: %i\n", hi_subset.size());
  
  IntegerVector perm_y = sortedIndex(y2);
  x2 = x2[perm_y];
  y2 = y2[perm_y];
  IntegerVector y3 = compare_self(y2);
  IntegerVector y4 = cumsum(y3);
  //return y4;
  
  IntegerVector perm_x = sortedIndex(x2);
  x2 = x2[perm_x];
  y4 = y4[perm_x];
  IntegerVector x3 = compare_self(x2);
  IntegerVector x4 = cumsum(x3);

  //return x4;
  IntegerVector obs = compare_both(x4, y4);
  IntegerVector cnt = diff(which_notzero(obs));
  int dis = kendall_discordant(x4, y4);
  
  return dis;
}


/*** R
x = c(12, 2, 1, 12, 2)
y = c(1, 4, 7, 1, 0)

t1 = ici_kendallt_mergesort(x, y)
t1
#compare_self(c(0, 1, 1, 4, 7))
*/
