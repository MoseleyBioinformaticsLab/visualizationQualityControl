#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

double signC(double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}

//' Calculates ici-kendall-tau
//' 
//' @param x numeric vector
//' @param y numeric vector
//' @param perspective should we consider the "local" or "global" perspective?
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
//' @name kendallt
//' @export
//' @useDynLib visualizationQualityControl
//' @return kendall tau correlation
// [[Rcpp::export]]
double ici_kendallt(NumericVector x, NumericVector y, String perspective = "local", String output = "simple") {
  
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
  
  x2 = x[!is_na(x)];
  y2 = y[!is_na(y)];
  
  double min_value = min(NumericVector::create(min(x2), min(y2)));
  double na_value = min_value - (0.1 * abs(min_value));
  //Rprintf("%f \n", min_value);
  //Rprintf("%f \n", na_value);
  
  x2 = clone(x);
  y2 = clone(y);
  x2[is_na(x)] = na_value;
  y2[is_na(y)] = na_value;
  
  
  double n_entry = x.size();
  
  if (n_entry < 2) {
    return NA_REAL;
  }
  
  if (perspective == "global") {
    for (int i = 0; i < (n_entry - 1); i++) {
      for (int j = (i+1); j < n_entry; j++) {
        sum_concordant+= (signC(x2[i] - x2[j]) * signC(y2[i] - y2[j])) > 0;
        sum_discordant+= (signC(x2[i] - x2[j]) * signC(y2[i] - y2[j])) < 0;
        
        sum_tied_x+= ((x2[i] != na_value) && (x2[j] != na_value) && (x2[i] == x2[j]) && ((y2[i] != na_value) && (y2[j] != na_value) && (y2[i] != y2[j])));
        sum_tied_y+= ((x2[i] != na_value) && (x2[j] != na_value) && (x2[i] != x2[j]) && ((y2[i] != na_value) && (y2[j] != na_value) && (y2[i] == y2[j])));
        sum_tied_x_na+= (x2[i] == na_value) && (x2[j] == na_value) && ((y2[i] != na_value) | (y2[j] != na_value));
        sum_tied_y_na+= ((x2[i] != na_value) | (x2[j] != na_value)) && (y2[i] == na_value) && (y2[j] == na_value);
        sum_all_na+= (x2[i] == na_value) && (x2[j] == na_value) && (y2[i] == na_value) && (y2[j] == na_value);
      }
    }
  } else {
    for (int i = 0; i < (n_entry - 1); i++) {
      for (int j = (i+1); j < n_entry; j++) {
        reject_concordant = ((x2[i] != na_value) && (x2[j] == na_value) && (y2[i] != na_value) && (y2[j] == na_value)) ||                                             //            ## 7
          ((x[i] == na_value) && (x[j] != na_value) && (y[i] == na_value) && (y[j] != na_value));
        reject_discordant = ((x2[i] != na_value) && (x2[j] == na_value)  && (y2[i] == na_value) && (y2[j] != na_value)) ||
          ((x2[i] == na_value)  && (x2[j] != na_value)  && (y2[i] != na_value) && (y2[j] == na_value));
        
        if (!reject_concordant) {
          sum_concordant+= (signC(x2[i] - x2[j]) * signC(y2[i] - y2[j])) > 0;
        }
        
        if (!reject_discordant) {
          sum_discordant+= (signC(x2[i] - x2[j]) * signC(y2[i] - y2[j])) < 0;
        }
        
        
        sum_tied_x+= ((x2[i] != na_value) && (x2[j] != na_value) && (x2[i] == x2[j]) && ((y2[i] != na_value) && (y2[j] != na_value) && (y2[i] != y2[j])));
        sum_tied_y+= ((x2[i] != na_value) && (x2[j] != na_value) && (x2[i] != x2[j]) && ((y2[i] != na_value) && (y2[j] != na_value) && (y2[i] == y2[j])));
        sum_tied_x_na+= (x2[i] == na_value) && (x2[j] == na_value) && ((y2[i] != na_value) | (y2[j] != na_value));
        sum_tied_y_na+= ((x2[i] != na_value) | (x2[j] != na_value)) && (y2[i] == na_value) && (y2[j] == na_value);
        sum_all_na+= (x2[i] == na_value) && (x2[j] == na_value) && (y2[i] == na_value) && (y2[j] == na_value);
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
  if (output != "simple") {
    Rprintf("n_entry: %f \n", n_entry);
    Rprintf("x_ties: %f \n", sum_tied_x);
    Rprintf("x_na_ties: %f \n", sum_tied_x_na);
    Rprintf("sum_x_ties: %f \n", sum_x_ties);
    Rprintf("sum_y_ties: %f \n", sum_y_ties);
    Rprintf("sum_na_ties: %f \n", sum_all_na);
    Rprintf("half_sum_na_ties: %f \n", half_sum_na_ties);
    Rprintf("sum_concordant: %f \n", sum_concordant);
    Rprintf("sum_discordant: %f \n", sum_discordant);
    Rprintf("k_numerator: %f \n", k_numerator);
    Rprintf("k_denominator: %f \n", k_denominator);
  }
  
  // 
  k_tau = k_numerator / k_denominator;
  
  return k_tau;
}


//' Calculates ici-kendall-tau
//' 
//' @param x numeric vector
//' @param y numeric vector
//' @param perspective should we consider the "local" or "global" perspective?
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
//' @name kendallt
//' @useDynLib visualizationQualityControl
//' @return kendall tau correlation
// [[Rcpp::export]]
double ici_ref_kendallt(NumericVector x, NumericVector y, String perspective = "local", String output = "simple") {
  
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
  
  LogicalVector matching_na;
  //double n_matching_na;
  
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
  
  if (perspective == "global") {
    for (int i = 0; i < (n_entry - 1); i++) {
      for (int j = (i+1); j < n_entry; j++) {
        sum_concordant+= (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||     //3 1
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||   // ## 2
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && !NumericVector::is_na(y[i]) && NumericVector::is_na(y[j])) ||                             // ## 3
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j])) ||                            //  ## 4
          (!NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||                            //  ## 5
          (NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||                            //  ## 6
          (!NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && !NumericVector::is_na(y[i]) && NumericVector::is_na(y[j])) ||                                             //            ## 7
          (NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]));
    
        sum_discordant+= (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && !NumericVector::is_na(y[i]) && NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) && NumericVector::is_na(x[j])  && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||
          (NumericVector::is_na(x[i])  && !NumericVector::is_na(x[j])  && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||
          (!NumericVector::is_na(x[i]) && NumericVector::is_na(x[j])  && NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j])) ||
          (NumericVector::is_na(x[i])  && !NumericVector::is_na(x[j])  && !NumericVector::is_na(y[i]) && NumericVector::is_na(y[j]));
        
        sum_tied_x+= (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] == x[j]) && (!NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] != y[j])));
        sum_tied_y+= (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] != x[j]) && (!NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] == y[j])));
        sum_tied_x_na+= NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && (!NumericVector::is_na(y[i]) | !NumericVector::is_na(y[j]));
        sum_tied_y_na+= (!NumericVector::is_na(x[i]) | !NumericVector::is_na(x[j])) && NumericVector::is_na(y[i]) && NumericVector::is_na(y[j]);
        sum_all_na+= NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && NumericVector::is_na(y[i]) && NumericVector::is_na(y[j]);
      }
    }
  } else {
    for (int i = 0; i < (n_entry - 1); i++) {
      for (int j = (i+1); j < n_entry; j++) {
        sum_concordant+= 
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||    // #3 1
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||  // ## 2
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && !NumericVector::is_na(y[i]) && NumericVector::is_na(y[j])) ||                             // ## 3
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j])) ||                             // ## 4
          (!NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||                             // ## 5
          (NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])); 
        
        sum_discordant+= 
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j])) ||
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] > x[j]) && NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] < x[j]) && !NumericVector::is_na(y[i]) && NumericVector::is_na(y[j])) ||
          (!NumericVector::is_na(x[i]) && NumericVector::is_na(x[j])  && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] < y[j])) ||
          (NumericVector::is_na(x[i])  && !NumericVector::is_na(x[j])  && !NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] > y[j]));
        
        sum_tied_x+= (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] == x[j]) && (!NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] != y[j])));
        sum_tied_y+= (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[j]) && (x[i] != x[j]) && (!NumericVector::is_na(y[i]) && !NumericVector::is_na(y[j]) && (y[i] == y[j])));
        sum_tied_x_na+= NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && (!NumericVector::is_na(y[i]) | !NumericVector::is_na(y[j]));
        sum_tied_y_na+= (!NumericVector::is_na(x[i]) | !NumericVector::is_na(x[j])) && NumericVector::is_na(y[i]) && NumericVector::is_na(y[j]);
        sum_all_na+= NumericVector::is_na(x[i]) && NumericVector::is_na(x[j]) && NumericVector::is_na(y[i]) && NumericVector::is_na(y[j]);
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
  
  
  if (output != "simple") {
    Rprintf("n_entry: %f \n", n_entry);
    Rprintf("x_ties: %f \n", sum_tied_x);
    Rprintf("x_na_ties: %f \n", sum_tied_x_na);
    Rprintf("sum_x_ties: %f \n", sum_x_ties);
    Rprintf("sum_y_ties: %f \n", sum_y_ties);
    Rprintf("sum_na_ties: %f \n", sum_all_na);
    Rprintf("half_sum_na_ties: %f \n", half_sum_na_ties);
    Rprintf("sum_concordant: %f \n", sum_concordant);
    Rprintf("sum_discordant: %f \n", sum_discordant);
    Rprintf("k_numerator: %f \n", k_numerator);
    Rprintf("k_denominator: %f \n", k_denominator);
  }
  // 
  k_tau = k_numerator / k_denominator;
  
  return k_tau;
}

//' Calculates ici-kendall-tau matrix
//' 
//' @param x data matrix
//' @param perspective should we consider the "local" or "global" perspective?
//' 
//'
//' @name kendallt_matrix
//' @export
//' @useDynLib visualizationQualityControl
//' @return kendall tau correlation
// [[Rcpp::export]]
NumericMatrix ici_kendall_matrix(NumericMatrix &x, String perspective = "local") {
  int x_col = x.ncol();
  
  NumericMatrix cor_matrix(x_col);
  
  for (int i = 0; i < x_col; i++) {
    NumericVector x_i = x(_ , i);
    for (int j = i; j < x_col; j++) {
      NumericVector x_j = x(_ , j);
      double kendall_cor = ici_kendallt(x_i, x_j, perspective);
      cor_matrix(i, j) = kendall_cor;
      cor_matrix(j, i) = kendall_cor;
    }
  } 
  
  return cor_matrix;
  
}
