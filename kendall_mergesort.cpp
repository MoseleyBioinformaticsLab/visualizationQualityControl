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

// [[Rcpp::export]]
NumericVector count_rank_tie(IntegerVector ranks){
  
  LogicalVector dup_ranks(ranks.size());
  dup_ranks = duplicated(ranks);
  IntegerVector ranks2 = ranks[dup_ranks];
  IntegerVector number_tied;
  number_tied = table(ranks2) + 1;
  
  NumericVector counts(3);
  counts(0) = sum(number_tied * (number_tied - 1)) / 2;
  counts(1) = sum(number_tied * (number_tied - 1) * (number_tied - 2)) / 2;
  counts(2) = sum(number_tied * (number_tied - 1) * (2 * number_tied + 5));
  counts.names() = CharacterVector({"ntie", "t0", "t1"});
  
  return counts;
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
NumericVector ici_kendallt_mergesort(NumericVector x, NumericVector y, String perspective = "local", String alternative = "two.sided", String output = "simple") {
  
  if (x.length() != y.length()) {
    throw std::range_error("X and Y are not the same length!");
    exit(-1);
  }
  
  
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
  //Rprintf("n_entry: %i\n", n_entry);
  
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
  
  double ntie = sum(cnt * (cnt - 1)) / 2;
  // three values should be read as:
  // xtie, x0, and x1, and then same for y
  NumericVector x_counts = count_rank_tie(x4);
  double xtie = x_counts[0];
  double x0 = x_counts[1];
  double x1 = x_counts[2];
  
  NumericVector y_counts = count_rank_tie(y4);
  double ytie = y_counts[0];
  double y0 = y_counts[1];
  double y1 = y_counts[2];
  
  int tot = (n_entry * (n_entry - 1)) / 2;
  
  //Note that tot = con + dis + (xtie - ntie) + (ytie - ntie) + ntie
  //              = con + dis + xtie + ytie - ntie
  
  NumericVector k_res(2);
  k_res.names() = CharacterVector({"tau", "pvalue"});
  if ((xtie == tot) || (ytie == tot)) {
    return k_res;
  }
  
  double con_minus_dis = tot - xtie - ytie + ntie - 2 * dis;
  double tau = con_minus_dis / sqrt(tot - xtie) / sqrt(tot - ytie);
  if (tau > 1) {
    tau = 1;
  } else if (tau < -1) {
    tau = -1;
  }
  
  double m = n_entry * (n_entry - 1);
  //Rprintf("m: %f\n", m);
  double var = ((m * (2 * n_entry + 5) - x1 - y1) / 18 +
                (2 * xtie * ytie) / m + x0 * y0 / (9 * m * (n_entry - 2)));
  //Rprintf("var: %f\n", var);
  double s_adjusted = tau * sqrt(((m / 2) - xtie) * ((m / 2) - ytie));
  //Rprintf("s_adjusted: %f\n", s_adjusted);
  double s_adjusted2 = signC(s_adjusted) * (std::abs(s_adjusted) - 1);
  //Rprintf("s_adjusted2: %f\n", s_adjusted2);
  z_b[0] = s_adjusted2 / sqrt(var);

  if (alternative == "less") {
    k_res[1] = pnorm(z_b, 0.0, 1.0)[0];
  } else if (alternative == "greater") {
    k_res[1] = pnorm(z_b, 0.0, 1.0, false, false)[0];
  } else if (alternative == "two.sided") {
    NumericVector p_res (2);
    p_res[0] = pnorm(z_b, 0.0, 1.0)[0];
    p_res[1] = pnorm(z_b, 0.0, 1.0, false)[0];
    k_res[1] = 2 * min(p_res);
  }
  k_res[0] = tau;
  
  
  return k_res;
}


/*** R
x = c(12, 2, 1, 12, 2)
y = c(1, 4, 7, 1, 0)

t1 = ici_kendallt_mergesort(x, y)
t1

t2 = visualizationQualityControl::ici_kendallt(x, y, output = "crap")
*/
