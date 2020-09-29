use itertools::Itertools;

pub type Numeric = f64;

#[derive(Debug)]
pub enum Perspective {
    Local,
    Global,
}

/// Calculates ici-kendall-tau
pub fn ici_kendall_tau(x: Vec<Numeric>, y: Vec<Numeric>, perspective: Perspective) -> Numeric {
    let mut sum_concordant: u64 = 0;
    let mut sum_discordant: u64 = 0;
    let mut sum_x_ties: u64 = 0;
    let mut sum_y_ties: u64 = 0;
    let mut sum_tied_x: u64 = 0;
    let mut sum_tied_y: u64 = 0;
    let mut sum_tied_x_na: u64 = 0;
    let mut sum_tied_y_na: u64 = 0;
    let mut sum_all_na: u64 = 0;
    let is_concordant: bool;
    let is_discordant: bool;
    let is_tied_x: bool;
    let is_tied_y: bool;
    let is_tied_x_na: bool;
    let is_tied_y_na: bool;
    let is_all_na: bool;
    let k_numerator: f64;
    let k_denominator: f64;
    let k_tau: f64;

    let mut x: Vec<f64> = x;
    let mut y: Vec<f64> = y;

    let mut matching_na = (&x)
        .iter()
        .zip((&y).iter())
        .map(|(x, y)| x.is_nan() & y.is_nan())
        // .collect_vec()
        .into_iter();

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

    match perspective {
        Perspective::Global => todo!(),
        Perspective::Local => {
            for ((xi, yi), (xj, yj)) in x.iter().zip(y.iter()).tuple_combinations() {
                sum_concordant += ((!(xi.is_nan()) && !(xj.is_nan()) && (xi > xj) && !(yi.is_nan()) && !(yj.is_nan()) && (yi > yj)) ||    // #3 1
                    (!(xi.is_nan()) && !(xj.is_nan()) && (xi < xj) && !(yi.is_nan()) && !(yj.is_nan()) && (yi < yj)) ||  // ## 2
                    (!(xi.is_nan()) && !(xj.is_nan()) && (xi > xj) && !(yi.is_nan()) && (yj.is_nan())) ||                             // ## 3
                    (!(xi.is_nan()) && !(xj.is_nan()) && (xi < xj) && (yi.is_nan()) && !(yj.is_nan())) ||                             // ## 4
                    (!(xi.is_nan()) && (xj.is_nan()) && !(yi.is_nan()) && !(yj.is_nan()) && (yi > yj)) ||                             // ## 5
                    ((xi.is_nan()) && !(xj.is_nan()) && !(yi.is_nan()) && !(yj.is_nan()) && (yi < yj)))
                    as u64;

                sum_discordant += ((!(xi.is_nan())
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
                        && (yi > yj))) as u64;

                sum_tied_x += (!(xi.is_nan())
                    && !(xj.is_nan())
                    // && (xi == xj)
                    && float_eq::float_eq!(xi, xj, ulps <= 6)
                    // && (!(yi.is_nan()) && !(yj.is_nan()) && (yi != yj)))
                    && (!(yi.is_nan()) && !(yj.is_nan()) && float_eq::float_ne!(yi, yj, ulps <= 6)))
                    as u64;
                sum_tied_y += (!(xi.is_nan())
                    && !(xj.is_nan())
                    // && (xi != xj)
                    && float_eq::float_eq!(xi, xj, ulps <= 6)
                    // && (!(yi.is_nan()) && !(yj.is_nan()) && (yi == yj)))
                    && (!(yi.is_nan()) && !(yj.is_nan()) && float_eq::float_ne!(yi, yj, ulps <= 6)))
                    as u64;
                sum_tied_x_na +=
                    ((xi.is_nan()) && (xj.is_nan()) && (!(yi.is_nan()) | !(yj.is_nan()))) as u64;
                sum_tied_y_na +=
                    ((!(xi.is_nan()) | !(xj.is_nan())) && (yi.is_nan()) && (yj.is_nan())) as u64;
                sum_all_na +=
                    ((xi.is_nan()) && (xj.is_nan()) && (yi.is_nan()) && (yj.is_nan())) as u64;
            }
        }
    }

    let half_sum_na_ties = sum_all_na / 2;
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use rand_distr::StandardNormal;

    #[test]
    fn test_speed_usecase() {
        let mut rng = StdRng::seed_from_u64(20200929);
        let n = 1000;
        let x = StandardNormal.sample_iter(&mut rng).take(n).collect();
        let y = StandardNormal.sample_iter(rng).take(n).collect();
        println!("{}", ici_kendall_tau(x, y, Perspective::Local));
        // The benchmark to beat (on this machine) is 36 ms (or preferably 11 ms)
    }
}
