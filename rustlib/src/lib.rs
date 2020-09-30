// use itertools::Itertools;
use rayon::prelude::*;

pub type Numeric = f64;

#[derive(Debug)]
pub enum Perspective {
    Local,
    Global,
}

/// Calculates ici-kendall-tau
///
/// ### Note
///
pub fn ici_kendall_tau(x: Vec<Numeric>, y: Vec<Numeric>, perspective: Perspective) -> Numeric {
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

    #[derive(derive_more::Add, derive_more::Sum)]
    struct Term {
        sum_concordant: u64,
        sum_discordant: u64,
        sum_tied_x: u64,
        sum_tied_y: u64,
        sum_tied_x_na: u64,
        sum_tied_y_na: u64,
        sum_all_na: u64,
    }

    let terms = match perspective {
        Perspective::Global => todo!(),
        Perspective::Local => {
            (0..n_entry - 1)
                .into_par_iter()
                .map(|x| ((x + 1)..n_entry).into_par_iter().map(move |xx| (x, xx)))
                .flatten()
                .map(|(i, j)|
                    {
                        let xi = x[i];
                        let xj = x[j];
                        let yi = y[i];
                        let yj = y[j];
                        // {
                        Term {
                            sum_concordant: ((!(xi.is_nan()) && !(xj.is_nan()) && (xi > xj) && !(yi.is_nan()) && !(yj.is_nan()) && (yi > yj)) ||    // #3 1
                                (!(xi.is_nan()) && !(xj.is_nan()) && (xi < xj) && !(yi.is_nan()) && !(yj.is_nan()) && (yi < yj)) ||  // ## 2
                                (!(xi.is_nan()) && !(xj.is_nan()) && (xi > xj) && !(yi.is_nan()) && (yj.is_nan())) ||                             // ## 3
                                (!(xi.is_nan()) && !(xj.is_nan()) && (xi < xj) && (yi.is_nan()) && !(yj.is_nan())) ||                             // ## 4
                                (!(xi.is_nan()) && (xj.is_nan()) && !(yi.is_nan()) && !(yj.is_nan()) && (yi > yj)) ||                             // ## 5
                                ((xi.is_nan()) && !(xj.is_nan()) && !(yi.is_nan()) && !(yj.is_nan()) && (yi < yj)))
                                as u64,

                            sum_discordant: ((!(xi.is_nan())
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
                                && (yi > yj))) as u64,

                            sum_tied_x: (!(xi.is_nan())
                                && !(xj.is_nan())
                                // && (xi == xj)
                                && float_eq::float_eq!(xi, xj, ulps <= 4)
                                // && (!(yi.is_nan()) && !(yj.is_nan()) && (yi != yj)))
                                && (!(yi.is_nan()) && !(yj.is_nan()) && float_eq::float_ne!(yi, yj, ulps <= 4)))
                                as u64,
                            sum_tied_y: (!(xi.is_nan())
                                && !(xj.is_nan())
                                // && (xi != xj)
                                && float_eq::float_ne!(xi, xj, ulps <= 4)
                                // && (!(yi.is_nan()) && !(yj.is_nan()) && (yi == yj)))
                                && (!(yi.is_nan()) && !(yj.is_nan()) && float_eq::float_eq!(yi, yj, ulps <= 4)))
                                as u64,
                            sum_tied_x_na:
                            ((xi.is_nan()) && (xj.is_nan()) && (!(yi.is_nan()) | !(yj.is_nan()))) as u64,
                            sum_tied_y_na:
                            ((!(xi.is_nan()) | !(xj.is_nan())) && (yi.is_nan()) && (yj.is_nan())) as u64,
                            sum_all_na:
                            ((xi.is_nan()) && (xj.is_nan()) && (yi.is_nan()) && (yj.is_nan())) as u64,
                        }
                    })
                .sum::<Term>()
        }
    };

    let half_sum_na_ties = terms.sum_all_na / 2;
    match perspective {
        Perspective::Global => {
            sum_x_ties = terms.sum_tied_x + terms.sum_tied_x_na + half_sum_na_ties;
            sum_y_ties = terms.sum_tied_y + terms.sum_tied_y_na + half_sum_na_ties;
        }
        Perspective::Local => {
            sum_x_ties = terms.sum_tied_x;
            sum_y_ties = terms.sum_tied_y;
        }
    }

    k_numerator = terms.sum_concordant as f64 - terms.sum_discordant as f64;
    k_denominator = (terms.sum_discordant + terms.sum_concordant + sum_x_ties + sum_y_ties) as f64;

    k_tau = k_numerator / k_denominator;
    k_tau
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

        println!("{}", ici_kendall_tau(x, y, Perspective::Local));
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
}
