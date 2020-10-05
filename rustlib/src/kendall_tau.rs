use crate::Numeric;
use itertools::Itertools;

#[derive(Debug, PartialEq)]
pub enum Perspective {
    Local,
    Global,
}

/// Calculates ici-kendall-tau
///
/// ### Note
///
pub fn ici_kendall_tau2(x: Vec<Numeric>, y: Vec<Numeric>, perspective: Perspective) -> Numeric {
    let sum_x_ties: u64;
    let sum_y_ties: u64;
    // let mut sum_concordant: u64 = 0;
    // let mut sum_discordant: u64 = 0;
    // let mut sum_tied_x: u64 = 0;
    // let mut sum_tied_y: u64 = 0;
    // let mut sum_tied_x_na: u64 = 0;
    // let mut sum_tied_y_na: u64 = 0;
    // let mut sum_all_na: u64 = 0;
    let _is_concordant: bool;
    let _is_discordant: bool;
    let _is_tied_x: bool;
    let _is_tied_y: bool;
    let _is_tied_x_na: bool;
    let _is_tied_y_na: bool;
    let _is_all_na: bool;
    let k_numerator: f64;
    let k_denominator: f64;
    let k_tau: f64;

    let mut x: Vec<f64> = x;
    let mut y: Vec<f64> = y;

    let matching_na = (&x)
        .iter()
        .zip((&y).iter())
        .map(|(x, y)| x.is_nan() & y.is_nan());

    if let Perspective::Local = perspective {
        // remove the values, that are missing across both vectors.
        let xy = x
            .iter()
            .zip(y.iter())
            .zip(matching_na)
            .filter(|((_x, _y), flag)| !*flag)
            .map(|((x, y), _)| (*x, *y));
        let (xx, yy) = xy.unzip();
        x = xx;
        y = yy;
    }

    let n_entry = x.len();
    debug_assert_eq!(x.len(), y.len());
    if n_entry < 2 {
        return f64::NAN;
    }

    use wide::{u64x2, u64x4};

    #[derive(derive_more::Add, derive_more::Sum)]
    struct Term {
        // sum_concordant: u64,
        // sum_discordant: u64,
        sum_con_dis_cordant: u64x2,

        // sum_tied_x: u64,
        // sum_tied_y: u64,
        // sum_tied_x_na: u64,
        // sum_tied_y_na: u64,
        sum_tied_xy_na: u64x4,

        sum_all_na: u64,
    }
    // sum_concordant: u64,
    // sum_discordant: u64,
    let mut sum_con_dis_cordant: u64x2 = Default::default();

    // sum_tied_x: u64,
    // sum_tied_y: u64,
    // sum_tied_x_na: u64,
    // sum_tied_y_na: u64,
    let mut sum_tied_xy_na: u64x4 = Default::default();

    let mut sum_all_na: u64x2 = Default::default();

    match perspective {
        Perspective::Global => todo!(),
        Perspective::Local => {
            (0..n_entry - 1)
                // .into_par_iter()
                .map(|x| ((x + 1)..n_entry).into_iter().map(move |xx| (x, xx)))
                .flatten()
                .map(|(i,j)| ( x[i],x[j],y[i],y[j],))
                .for_each(|(xi, xj, yi, yj)| {
                    sum_con_dis_cordant += (&[((!(xi.is_nan()) && !(xj.is_nan()) && (xi > xj) && !(yi.is_nan()) && !(yj.is_nan()) && (yi > yj)) ||    // #3 1
                        (!(xi.is_nan()) && !(xj.is_nan()) && (xi < xj) && !(yi.is_nan()) && !(yj.is_nan()) && (yi < yj)) ||  // ## 2
                        (!(xi.is_nan()) && !(xj.is_nan()) && (xi > xj) && !(yi.is_nan()) && (yj.is_nan())) ||                             // ## 3
                        (!(xi.is_nan()) && !(xj.is_nan()) && (xi < xj) && (yi.is_nan()) && !(yj.is_nan())) ||                             // ## 4
                        (!(xi.is_nan()) && (xj.is_nan()) && !(yi.is_nan()) && !(yj.is_nan()) && (yi > yj)) ||                             // ## 5
                        ((xi.is_nan()) && !(xj.is_nan()) && !(yi.is_nan()) && !(yj.is_nan()) && (yi < yj)))
                        as u64,
                        ((!(xi.is_nan())
                            && !(xj.is_nan())
                            && (xi > xj)
                            && !(yi.is_nan())
                            && !(yj.is_nan())
                            && (yi < yj))
                            || (!(xi.is_nan())
                            && !(xj.is_nan())
                            && (xi < xj)
                            && !(yi.is_nan())
                            && !(yj.is_nan())
                            && (yi > yj))
                            || (!(xi.is_nan())
                            && !(xj.is_nan())
                            && (xi > xj)
                            && (yi.is_nan())
                            && !(yj.is_nan()))
                            || (!(xi.is_nan())
                            && !(xj.is_nan())
                            && (xi < xj)
                            && !(yi.is_nan())
                            && (yj.is_nan()))
                            || (!(xi.is_nan())
                            && (xj.is_nan())
                            && !(yi.is_nan())
                            && !(yj.is_nan())
                            && (yi < yj))
                            || ((xi.is_nan())
                            && !(xj.is_nan())
                            && !(yi.is_nan())
                            && !(yj.is_nan())
                            && (yi > yj))) as u64].into());

                    sum_tied_xy_na += (&[
                        // sum_tied_x:
                        (!(xi.is_nan())
                            && !(xj.is_nan())
                            // && (xi == xj)
                            && float_eq::float_eq!(xi, xj, ulps <= 4)
                            // && (!(yi.is_nan()) && !(yj.is_nan()) && (yi != yj)))
                            && (!(yi.is_nan()) && !(yj.is_nan()) && float_eq::float_ne!(yi, yj, ulps <= 4)))
                            as u64,
                        // sum_tied_y:
                        (!(xi.is_nan())
                            && !(xj.is_nan())
                            // && (xi != xj)
                            && float_eq::float_ne!(xi, xj, ulps <= 4)
                            // && (!(yi.is_nan()) && !(yj.is_nan()) && (yi == yj)))
                            && (!(yi.is_nan()) && !(yj.is_nan()) && float_eq::float_eq!(yi, yj, ulps <= 4)))
                            as u64,
                        // sum_tied_x_na:
                        ((xi.is_nan()) && (xj.is_nan()) && (!(yi.is_nan()) | !(yj.is_nan()))) as u64,
                        // sum_tied_y_na:
                        ((!(xi.is_nan()) | !(xj.is_nan())) && (yi.is_nan()) && (yj.is_nan())) as u64,
                    ].into());

                    sum_all_na += (&(
                        ((xi.is_nan()) && (xj.is_nan()) && (yi.is_nan()) && (yj.is_nan())) as u64).into());
                });
        }
    }
    let [sum_concordant, sum_discordant]: [u64; 2] = sum_con_dis_cordant.into();
    let [sum_tied_x, sum_tied_y, sum_tied_x_na, sum_tied_y_na]: [u64; 4] = sum_tied_xy_na.into();

    let [half_sum_na_ties, _]: [u64; 2] = sum_all_na.into();
    let half_sum_na_ties = half_sum_na_ties / 2;
    match perspective {
        Perspective::Global => {
            sum_x_ties = sum_tied_x + sum_tied_x_na + half_sum_na_ties;
            sum_y_ties = sum_tied_y + sum_tied_y_na + half_sum_na_ties;
        }
        Perspective::Local => {
            sum_x_ties = sum_tied_x;
            sum_y_ties = sum_tied_y;
        }
    }

    k_numerator = sum_concordant as f64 - sum_discordant as f64;
    k_denominator = (sum_discordant + sum_concordant + sum_x_ties + sum_y_ties) as f64;

    k_tau = k_numerator / k_denominator;
    k_tau
}

///
/// Mimics the sign function in R.
/// Note: this is not equivalent to `f64::signum`
#[must_use]
fn sign_c(x: Numeric) -> isize {
    if x > 0. {
        1
    } else if float_eq::float_eq!(x, 0., ulps <= 6) {
        0
    } else {
        -1
    }
}

///
/// `perspective` is to have either local or global perspective.
/// Note: There is a missing an `Output` parameter.
///
/// Panics if `x` and `y` are of unequal length.
#[must_use]
pub fn ici_kendall_tau(x: Vec<Numeric>, y: Vec<Numeric>, perspective: Perspective) -> Numeric {
    assert_eq!(x.len(), y.len());

    let mut sum_concordant: u64 = 0;
    let mut sum_discordant: u64 = 0;
    let sum_x_ties: Numeric;
    let sum_y_ties: Numeric;
    let mut sum_tied_x: u64 = 0;
    let mut sum_tied_y: u64 = 0;
    let mut sum_tied_x_na: u64 = 0;
    let mut sum_tied_y_na: u64 = 0;
    let mut sum_all_na: u64 = 0;

    // remove the values that are NA on both
    let mut x = x;
    let mut y = y;

    if let Perspective::Local = perspective {
        let (xx, yy) = x
            .iter()
            .zip_eq(y.iter())
            .filter(|(x, y)| !(x.is_nan() && y.is_nan()))
            .unzip();
        x = xx;
        y = yy;
    }

    let n_na_x = x.iter().map(|x| x.is_nan() as usize).sum::<usize>();
    let n_na_y = y.iter().map(|x| x.is_nan() as usize).sum::<usize>();

    if (x.len() == n_na_x) || (y.len() == n_na_y) {
        return 0.;
    }
    let na_value = {
        let x = x.iter().filter(|x| !x.is_nan());
        let y = y.iter().filter(|x| !x.is_nan());

        let min_value = *x.chain(y).min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        min_value - (0.1 * min_value.abs())
    };
    let x2 = x
        .iter()
        .map(|x| if x.is_nan() { na_value } else { *x })
        .collect_vec();
    let y2 = y
        .iter()
        .map(|y| if y.is_nan() { na_value } else { *y })
        .collect_vec();

    let n_entry = x.len(); // `x` is equal length to `y` and
    if n_entry < 2 {
        return 0.;
    }

    #[warn(clippy::float_cmp)]
    for i in 0..(n_entry - 1) {
        for j in (i + 1)..n_entry {
            // let x2i = x2[i];
            // let x2j = x2[j];
            // let y2i = y2[i];
            // let y2j = y2[j];
            let x2i = *unsafe { x2.get_unchecked(i) };
            let x2j = *unsafe { x2.get_unchecked(j) };
            let y2i = *unsafe { y2.get_unchecked(i) };
            let y2j = *unsafe { y2.get_unchecked(j) };

            let reject_concordant = ((x2i != na_value) && (x2j == na_value) && (y2i != na_value) && (y2j == na_value)) ||                                             //            ## 7
                ((x2i == na_value) && (x2j != na_value) && (y2i == na_value) && (y2j != na_value));
            let reject_discordant =
                ((x2i != na_value) && (x2j == na_value) && (y2i == na_value) && (y2j != na_value))
                    || ((x2i == na_value)
                        && (x2j != na_value)
                        && (y2i != na_value)
                        && (y2j == na_value));

            if perspective == Perspective::Global || !reject_concordant {
                sum_concordant += ((sign_c(x2i - x2j) * sign_c(y2i - y2j)) > 0) as u64;
            }

            if perspective == Perspective::Global || !reject_discordant {
                sum_discordant += ((sign_c(x2i - x2j) * sign_c(y2i - y2j)) < 0) as u64;
            }

            sum_tied_x += ((x2i != na_value)
                && (x2j != na_value)
                // && (x2i == x2j)
                && float_eq::float_eq!(x2i, x2j, ulps <= 6)
                && ((y2i != na_value) && (y2j != na_value) &&
                // (y2i != y2j)
                float_eq::float_ne!(y2i, y2j, ulps<=6)
            )) as u64;
            sum_tied_y += ((x2i != na_value)
                && (x2j != na_value)
                // && (x2i != x2j)
                && float_eq::float_ne!(x2i, x2j, ulps<=6)
                && ((y2i != na_value) && (y2j != na_value) &&
                // (y2i == y2j)
                float_eq::float_eq!(y2i, y2j, ulps<=6)
            )) as u64;
            sum_tied_x_na += ((x2i == na_value)
                && (x2j == na_value)
                && ((y2i != na_value) | (y2j != na_value))) as u64;
            sum_tied_y_na += (((x2i != na_value) | (x2j != na_value))
                && (y2i == na_value)
                && (y2j == na_value)) as u64;
            sum_all_na +=
                ((x2i == na_value) && (x2j == na_value) && (y2i == na_value) && (y2j == na_value))
                    as u64;
        }
    }

    let half_sum_na_ties: Numeric = sum_all_na as f64 / 2.;

    match perspective {
        Perspective::Local => {
            sum_x_ties = sum_tied_x as f64;
            sum_y_ties = sum_tied_y as f64;
        }
        Perspective::Global => {
            sum_x_ties = sum_tied_x as f64 + sum_tied_x_na as f64 + half_sum_na_ties;
            sum_y_ties = sum_tied_y as f64 + sum_tied_y_na as f64 + half_sum_na_ties;
        }
    };

    let k_numerator = sum_concordant as f64 - sum_discordant as f64;
    let k_denominator = sum_discordant as f64 + sum_concordant as f64 + sum_x_ties + sum_y_ties;

    // debugging
    // if (output != "simple") {
    //     Rprintf("min_value: %f \n", min_value);
    //     Rprintf("na_value: %f \n", na_value);
    //
    //     Rprintf("n_entry: %f \n", n_entry);
    //     Rprintf("x_ties: %f \n", sum_tied_x);
    //     Rprintf("x_na_ties: %f \n", sum_tied_x_na);
    //     Rprintf("sum_x_ties: %f \n", sum_x_ties);
    //     Rprintf("sum_y_ties: %f \n", sum_y_ties);
    //     Rprintf("sum_na_ties: %f \n", sum_all_na);
    //     Rprintf("half_sum_na_ties: %f \n", half_sum_na_ties);
    //     Rprintf("sum_concordant: %f \n", sum_concordant);
    //     Rprintf("sum_discordant: %f \n", sum_discordant);
    //     Rprintf("k_numerator: %f \n", k_numerator);
    //     Rprintf("k_denominator: %f \n", k_denominator);
    // }

    // `k_tau`
    if float_eq::float_eq!(k_denominator, 0., ulps <= 6) {
        0.
    } else {
        k_numerator / k_denominator
    }
}

/*
double ici_kendallt(NumericVector x, NumericVector y, String perspective = "local", String output = "simple") {

  if (perspective == "global") {
    for (int i = 0; i < (n_entry - 1); i++) {
      for (int j = (i+1); j < n_entry; j++) {
        sum_concordant+= (signC(x2i - x2[j]) * signC(y2[i] - y2[j])) > 0;
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
          ((x2i == na_value) && (x[j] != na_value) && (y[i] == na_value) && (y[j] != na_value));
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
    Rprintf("min_value: %f \n", min_value);
    Rprintf("na_value: %f \n", na_value);

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

  if (k_denominator == 0) {
    k_tau = 0;
  } else {
    k_tau = k_numerator / k_denominator;
  }

  return k_tau;
}
*/

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use rand_distr::StandardNormal;

    #[test]
    fn test_speed_usecase() {
        // let mut rng = StdRng::seed_from_u64(20200929);
        let rng = thread_rng();
        let n = 1000;
        let x: Vec<f64> = StandardNormal.sample_iter(rng).take(n).collect();
        let y: Vec<f64> = StandardNormal.sample_iter(rng).take(n).collect();

        println!("{:?}\n{:?}", &x[..10], &y[..10]);

        // println!("{}", ici_kendall_tau(x, y, Perspective::Local));
        // The benchmark to beat (on this machine) is 36 ms (or preferably 11 ms)
    }
    #[test]
    fn test_standard_normal_many_times() {
        // let mut rng = StdRng::seed_from_u64(20200929);
        let rng = thread_rng();
        let reps = 10;
        let n = 1000;
        for _rep in 0..reps {
            let x: Vec<f64> = StandardNormal.sample_iter(rng).take(n).collect();
            let y: Vec<f64> = StandardNormal.sample_iter(rng).take(n).collect();

            println!("{}", ici_kendall_tau(x, y, Perspective::Local));
        }
    }

    #[test]
    fn test_sign_c_against_default_sign_function() {
        // assert_eq!(sign_c(0_f64), 0_f64.signum() as isize); // this fails
        // assert_eq!(sign_c(-0_f64), (-0.0_f64).signum() as isize); // this fails
        assert_eq!(sign_c(1_f64), 1_f64.signum() as isize);
        assert_eq!(sign_c(2_f64), 2_f64.signum() as isize);
        assert_eq!(sign_c(-1_f64), -1_f64.signum() as isize);
        assert_eq!(sign_c(-2_f64), -2_f64.signum() as isize);
    }
}
