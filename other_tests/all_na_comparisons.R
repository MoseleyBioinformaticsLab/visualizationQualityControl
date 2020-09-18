library(visualizationQualityControl)
library(furrr)
plan(multiprocess)

create_na_indices <- function(x = 20) {
  
  n_na = seq(1, x)
  
  where_na = purrr::map(n_na, function(in_na){
    na_comb = combn(x, in_na)
    asplit(na_comb, 2)
  })
  where_na = unlist(where_na, recursive = FALSE)
  
}

compare_positive_kt <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    visualizationQualityControl:::ref_kendallt(tmp_x, tmp_y, perspective = perspective)
  }, .progress = TRUE)
  
  tmp
}

compare_positive_kt_c <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dbl(where_na, function(use_na){
    #purrr::imap_dbl(where_na, function(use_na, .y){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    kendallt(tmp_x, tmp_y, perspective = perspective)
    #})
  }, .progress = TRUE)
  
  tmp
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_negative_kt <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  tmp = furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    y_na = n_entry - y_na + 1
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    visualizationQualityControl:::ref_kendallt(tmp_x, tmp_y, perspective = perspective)
  }, .progress = TRUE)
  tmp
  
}

compare_negative_kt_c <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  tmp = furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    y_na = n_entry - y_na + 1
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    kendallt(tmp_x, tmp_y, perspective = perspective)
  }, .progress = TRUE)
  
  tmp
}

x = seq(1, 10)
y = seq(1, 10)
y2 = seq(10, 1)

where_na = create_na_indices()

positive_kt = compare_positive_kt(x, y, where_na, perspective = "global")
positive_kt_c = compare_positive_kt_c(x, y, where_na, perspective = "global")
negative_kt = compare_negative_kt(x, y2, where_na, perspective = "global")
negative_kt_c = compare_negative_kt_c(x, y2, where_na, perspective = "global")

all.equal(positive_kt, positive_kt_c)
all.equal(negative_kt, negative_kt_c)


positive_kt_local = compare_positive_kt(x, y, where_na, perspective = "local")
positive_kt_c_local = compare_positive_kt_c(x, y, where_na, perspective = "local")
negative_kt_local = compare_negative_kt(x, y2, where_na, perspective = "local")
negative_kt_c_local = compare_negative_kt_c(x, y2, where_na, perspective = "local")

all.equal(positive_kt_local, positive_kt_c_local)
all.equal(negative_kt_local, negative_kt_c_local)
