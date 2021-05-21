library(reticulate)
use_python("/usr/local/bin/python3")
#py_config()

run_np_kt = function(x, y, py_stats, type = "global"){
  if (type == "local") {
    x_na = is.na(x)
    y_na = is.na(y)
    match_na = x_na & y_na
    x = x[!match_na]
    y = y[!match_na]
  }
  min_xy = min(c(x, y), na.rm = TRUE)
  na_replace = min_xy - 0.1
  x2 = x
  y2 = y
  x2[is.na(x)] = na_replace
  y2[is.na(y)] = na_replace
  
  np_out = py_stats$kendalltau(x2, y2)
  np_out
}

run_cor_test = function(x, y, type = "global"){
  if (type == "local") {
    x_na = is.na(x)
    y_na = is.na(y)
    match_na = x_na & y_na
    x = x[!match_na]
    y = y[!match_na]
  }
  min_xy = min(c(x, y), na.rm = TRUE)
  na_replace = min_xy - 0.1
  x2 = x
  y2 = y
  x2[is.na(x)] = na_replace
  y2[is.na(y)] = na_replace
  
  r_out = cor.test(x2, y2, method = "kendall", continuity = TRUE)
  r_out
}

py_scipy = import("scipy")
py_stats = import("scipy.stats")
x = rnorm(100)
y = x + 1
#x = round(x, 2)
#y = round(y, 2)

np1 = run_np_kt(x, y, py_stats)
np1$correlation
np1$pvalue

r1 = run_cor_test(x, y)
r1

library(visualizationQualityControl)
ici1 = ici_kendallt(x, y, perspective = "global")
ici1

x2 = x
x2[1:5] = x[100]
y2 = y
y2[10:15] = y[100]

np2 = run_np_kt(x2, y, py_stats)
np2$correlation
np2$pvalue

r2 = run_cor_test(x2, y)
r2

ici2 = ici_kendallt(x2, y, perspective = "global")
ici2

np3 = run_np_kt(x2, y2, py_stats)
np3$correlation
np3$pvalue

r3 = run_cor_test(x2, y2)
r3

ici3 = ici_kendallt(x2, y2, perspective = "global")
ici3

x3 = x
x3[1:20] = NA
y3 = y
y3[40:60] = NA

np4 = run_np_kt(x3, y, py_stats)
np4$correlation
np4$pvalue

r4 = run_cor_test(x3, y)
r4

ici4 = ici_kendallt(x3, y, perspective = "global")
ici4


np5 = run_np_kt(x3, y3, py_stats)
np5$correlation
np5$pvalue

r5 = run_cor_test(x3, y3)
r5

ici5 = ici_kendallt(x3, y3, perspective = "global")
ici5

# very good thus far, everything tracks nicely
# Now lets try putting some in common to verify local vs global
y4 = y
y4[c(1:5, 35:40)] = NA

np4 = run_np_kt(x3, y4, py_stats)
np4$correlation
np4$pvalue


r4 = run_cor_test(x3, y4)
r4

ici4 = ici_kendallt(x3, y4, perspective = "global")
ici4

np5 = run_np_kt(x3, y4, py_stats, "local")
np5$correlation
np5$pvalue

r4 = run_cor_test(x3, y4, "local")
r4

ici5 = ici_kendallt(x3, y4, "local")
ici5
