use crate::Numeric;
use itertools::Itertools;
use ndarray::{Array1, ArrayView1};
#[cfg(feature = "rayon")]
use rayon::prelude::*;

use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Debug, PartialEq, Clone)]
pub enum Perspective {
    Local,
    Global,
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
pub fn ici_kendall_tau(
    x: ArrayView1<Numeric>,
    y: ArrayView1<Numeric>,
    perspective: Perspective,
) -> Numeric {
    assert_eq!(x.len(), y.len());

    let sum_x_ties: Numeric;
    let sum_y_ties: Numeric;
    // let mut sum_concordant: u64 = 0;
    // let mut sum_discordant: u64 = 0;
    // let mut sum_tied_x: u64 = 0;
    // let mut sum_tied_y: u64 = 0;
    // let mut sum_tied_x_na: u64 = 0;
    // let mut sum_tied_y_na: u64 = 0;
    // let mut sum_all_na: u64 = 0;
    let sum_concordant = AtomicU64::new(0);
    let sum_discordant = AtomicU64::new(0);
    let sum_tied_x = AtomicU64::new(0);
    let sum_tied_y = AtomicU64::new(0);
    let sum_tied_x_na = AtomicU64::new(0);
    let sum_tied_y_na = AtomicU64::new(0);
    let sum_all_na = AtomicU64::new(0);

    // remove the values that are NA on both
    // let mut x = x;
    // let mut y = y;

    // if let Perspective::Local = perspective {
    //     let matching_na = ndarray::Zip::from(&x)
    //         .and(&y)
    //         .apply_collect(|x, y| !(x.is_nan() && y.is_nan()));
    //     x = Zip::from(&x).and(&matching_na).apply_c
    //     // let (xx, yy) = x
    //     //     .into_iter()
    //     //     .zip_eq(y.into_iter())
    //     //     .filter(|(x, y)| !(x.is_nan() && y.is_nan()))
    //     //     .unzip();
    //     // x = xx;
    //     // y = yy;
    // }

    let mut n_na_x = x.iter().map(|x| x.is_nan() as usize).sum::<usize>();
    let mut n_na_y = y.iter().map(|x| x.is_nan() as usize).sum::<usize>();
    let mut matching_na: Array1<bool> = ndarray::Array1::<bool>::from_elem(x.dim(), false);
    let mut n_matching_na: usize = 0;
    if let Perspective::Local = perspective {
        matching_na = ndarray::Zip::from(&x)
            .and(&y)
            .apply_collect(|x, y| x.is_nan() && y.is_nan());
        n_matching_na = matching_na.map(|&x| x as usize).sum();
        n_na_x = n_na_x.saturating_sub(n_matching_na);
        n_na_y = n_na_y.saturating_sub(n_matching_na);

        if (x.len().saturating_sub(n_na_x) == n_na_x) || (y.len().saturating_sub(n_na_y) == n_na_y)
        {
            dbg!("this?");
            return 0.;
        }
    } else if (x.len() == n_na_x) || (y.len() == n_na_y) {
        dbg!("that?");
        return 0.;
    }

    let na_value = {
        let min_value = x
            .iter()
            .chain(y.iter())
            .min_by(|x, y| match x.partial_cmp(y) {
                None => match (x.is_nan(), y.is_nan()) {
                    (true, true) => std::cmp::Ordering::Equal,
                    (false, _) => std::cmp::Ordering::Less,
                    (_, false) => std::cmp::Ordering::Greater,
                },
                Some(a) => a,
            })
            .cloned()
            .unwrap();
        min_value - (0.1 * min_value.abs())
    };
    let x2 = x.iter().zip_eq(&matching_na).filter_map(|(x, flag)| {
        if !*flag {
            if x.is_nan() {
                Some(na_value)
            } else {
                Some(*x)
            }
        } else {
            None
        }
    });
    // .filter(|(x, flag)| !(**flag))
    // .map(|x| x.0)
    // .map(|x| if x.is_nan() { na_value } else { *x });
    let y2 = y.iter().zip_eq(&matching_na).filter_map(|(x, flag)| {
        if !*flag {
            if x.is_nan() {
                Some(na_value)
            } else {
                Some(*x)
            }
        } else {
            None
        }
    });
    // .filter(|(x, flag)| !(**flag))
    // .map(|x| x.0)
    // .map(|y| if y.is_nan() { na_value } else { *y });

    let n_entry = x.len(); // `x` is equal length to `y` and
    if n_entry < 2 {
        return 0.;
    }

    let pair_iter = x2.zip(y2).combinations(2);

    // #[cfg(feature = "rayon")]
    // let pair_iter = pair_iter.collect_vec().into_par_iter();
    // let pair_iter = pair_iter.par_chunks(250);

    #[warn(clippy::float_cmp)]
    pair_iter.for_each(|v| {
        let (x2i, y2i) = v[0];
        let (x2j, y2j) = v[1];

        let reject_concordant = ((x2i != na_value) && (x2j == na_value) && (y2i != na_value) && (y2j == na_value)) ||                                             //            ## 7
                ((x2i == na_value) && (x2j != na_value) && (y2i == na_value) && (y2j != na_value));
        let reject_discordant = ((x2i != na_value)
            && (x2j == na_value)
            && (y2i == na_value)
            && (y2j != na_value))
            || ((x2i == na_value) && (x2j != na_value) && (y2i != na_value) && (y2j == na_value));

        let sum_concordant_term = if perspective == Perspective::Global || !reject_concordant {
            ((sign_c(x2i - x2j) * sign_c(y2i - y2j)) > 0) as u64
        } else {
            Default::default()
        };

        let sum_discordant_term = if perspective == Perspective::Global || !reject_discordant {
            ((sign_c(x2i - x2j) * sign_c(y2i - y2j)) < 0) as u64
        } else {
            Default::default()
        };

        let sum_tied_x_term = ((x2i != na_value)
                && (x2j != na_value)
                // && (x2i == x2j)
                && float_eq::float_eq!(x2i, x2j, ulps <= 6)
                && ((y2i != na_value) && (y2j != na_value) &&
                // (y2i != y2j)
                float_eq::float_ne!(y2i, y2j, ulps<=6)
            )) as u64;
        let sum_tied_y_term = ((x2i != na_value)
                && (x2j != na_value)
                // && (x2i != x2j)
                && float_eq::float_ne!(x2i, x2j, ulps<=6)
                && ((y2i != na_value) && (y2j != na_value) &&
                // (y2i == y2j)
                float_eq::float_eq!(y2i, y2j, ulps<=6)
            )) as u64;
        let sum_tied_x_na_term = ((x2i == na_value)
            && (x2j == na_value)
            && ((y2i != na_value) | (y2j != na_value))) as u64;
        let sum_tied_y_na_term = (((x2i != na_value) | (x2j != na_value))
            && (y2i == na_value)
            && (y2j == na_value)) as u64;
        let sum_all_na_term =
            ((x2i == na_value) && (x2j == na_value) && (y2i == na_value) && (y2j == na_value))
                as u64;

        sum_concordant.fetch_add(sum_concordant_term, Ordering::Relaxed);
        sum_discordant.fetch_add(sum_discordant_term, Ordering::Relaxed);
        sum_tied_x.fetch_add(sum_tied_x_term, Ordering::Relaxed);
        sum_tied_y.fetch_add(sum_tied_y_term, Ordering::Relaxed);
        sum_tied_x_na.fetch_add(sum_tied_x_na_term, Ordering::Relaxed);
        sum_tied_y_na.fetch_add(sum_tied_y_na_term, Ordering::Relaxed);
        sum_all_na.fetch_add(sum_all_na_term, Ordering::Relaxed);
    });

    let sum_all_na = sum_all_na.load(Ordering::Relaxed);
    let half_sum_na_ties: Numeric = sum_all_na as f64 / 2.;

    let sum_tied_x = sum_tied_x.load(Ordering::Relaxed);
    let sum_tied_y = sum_tied_y.load(Ordering::Relaxed);
    let sum_tied_x_na = sum_tied_x_na.load(Ordering::Relaxed);
    let sum_tied_y_na = sum_tied_y_na.load(Ordering::Relaxed);
    let sum_concordant = sum_concordant.load(Ordering::Relaxed);
    let sum_discordant = sum_discordant.load(Ordering::Relaxed);

    //     println!(
    //         "sum_concordant = {}\n\
    // sum_discordant = {}\n\
    // sum_tied_x = {}\n\
    // sum_tied_y = {}\n\
    // sum_tied_x_na = {}\n\
    // sum_tied_y_na = {}\n\
    // sum_all_na = {}\n",
    //         sum_concordant,
    //         sum_discordant,
    //         sum_tied_x,
    //         sum_tied_y,
    //         sum_tied_x_na,
    //         sum_tied_y_na,
    //         sum_all_na
    //     );

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

    if float_eq::float_eq!(k_denominator, 0., ulps <= 6) {
        0.
    } else {
        k_numerator / k_denominator
    }
}

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
    #[ignore]
    fn test_standard_normal_many_times() {
        // let mut rng = StdRng::seed_from_u64(20200929);
        let rng = thread_rng();
        let reps = 10;
        let n = 1000;
        for _rep in 0..reps {
            let x: Vec<f64> = StandardNormal.sample_iter(rng).take(n).collect();
            let y: Vec<f64> = StandardNormal.sample_iter(rng).take(n).collect();

            println!(
                "{}",
                ici_kendall_tau(
                    Array1::from(x).view(),
                    Array1::from(y).view(),
                    Perspective::Local
                )
            );
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
