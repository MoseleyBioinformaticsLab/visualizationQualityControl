use crate::kendall_tau::{ici_kendall_tau, Perspective};
use crate::Numeric;
use itertools::Itertools;
use ndarray::{Array2, Axis};
use rand_distr::num_traits::Zero;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use std::collections::HashMap;

/// information-content-informed kendall tau
/// Given a data-matrix, computes the information-content-informed (ICI) Kendall-tau-b between
/// all samples.
///
/// - `data_matrix`: samples are rows, features are columns
/// - `exclude_na`: should NA values be treated as NA?
/// - `exclude_inf`: should Inf values be treated as NA?
/// - `exclude_0`: should zero values be treated as NA?
/// - `zero_value`: what is the actual zero value?
/// - `perspective`: how to treat missing data in denominator and ties, see details
/// - `scale_max`: should everything be scaled compared to the maximum correlation?
/// - `diag_not_na`: should the diagonal entries reflect how many entries in the sample were "good"?
///
/// Default call: `function(data_matrix, exclude_na = TRUE, exclude_inf = TRUE,
///     exclude_0 = TRUE, zero_value = 0, perspective = "global", scale_max = TRUE,
///     diag_good = TRUE, progress = FALSE)`
pub fn visqc_ici_kendall_tau(
    data_matrix: ndarray::Array2<Numeric>,
    exclude_na: bool,
    exclude_inf: bool,
    exclude_0: bool,
    _zero_value: f64,
    perspective: Perspective,
    scale_max: bool,
    diag_good: bool,
    _progress: bool,
) -> Array2<f64> {
    // it is assumed now that the data-matrix is features x observations (atleast accoring to the code)
    //FIXME: should panic if `data_matrix` doesn't have atleast 2 variables

    // these are all initially `false`
    let exclude_nothing = Array2::from_elem(data_matrix.dim(), false);

    let na_loc: Array2<bool> = if exclude_na {
        //FIXME: this should be a check for NA values.
        data_matrix.map(|v| v.is_nan())
    } else {
        exclude_nothing.clone()
    };
    let inf_loc: Array2<bool> = if exclude_inf {
        data_matrix.map(|x| x.is_infinite())
    } else {
        exclude_nothing.clone()
    };
    let zero_loc: Array2<bool> = if exclude_0 {
        data_matrix.map(|v| v.is_zero())
    } else {
        exclude_nothing
    };

    let exclude_loc: Array2<bool> = na_loc | zero_loc | inf_loc;
    //TODO: set these to be NA as well.
    let mut data_matrix = data_matrix;
    (ndarray::Zip::from(&mut data_matrix).and(&exclude_loc)).apply(|d, &f| {
        if f {
            //FIXME: this should be NA
            *d = std::f64::NAN;
        }
    });
    //TODO: maybe implement a progress bar {pbr} here?
    // Only reason this is not done immediately, is that the pbr-crate is made for io-stuff, etc.
    // and there is no readily available trait/impl to use with rayon/ndarray.
    let ncols = data_matrix.ncols();
    let cor_map_iter = (0..ncols)
        .map(|c| (c..ncols).map(move |x| (c, x)))
        .flatten();

    #[cfg(feature = "rayon")]
    let cor_map_iter = cor_map_iter.collect_vec().into_par_iter();

    let cor_map = cor_map_iter
        .map(|(icol, jcol)| {
            (
                (icol, jcol),
                ici_kendall_tau(
                    data_matrix.column(icol),
                    data_matrix.column(jcol),
                    perspective.clone(),
                ),
            )
        })
        .collect::<HashMap<(usize, usize), Numeric>>();

    let cor_matrix: Array2<Numeric> = Array2::from_shape_fn((ncols, ncols), |(i, j)| {
        // ensure that i<=j when accessing the correlation records
        let (i, j) = (i.min(j), i.max(j));
        // dbg!(i, j);
        *cor_map.get(&(i, j)).unwrap()
    });

    // calculate the max-cor value for use in scaling across multiple comparisons
    let n_observations = data_matrix.nrows() as u64;
    let n_na = exclude_loc.mapv(|x| x as u64).sum_axis(Axis(0));
    debug_assert_eq!(n_na.len(), data_matrix.ncols());
    let mut n_na = n_na.into_raw_vec();
    n_na.sort_unstable();
    let m_value = (n_na[0] + n_na[1]) / 2;
    let n_m = n_observations - m_value;
    let max_cor_denominator =
        statrs::function::factorial::binomial(n_m, 2) + (n_observations * m_value) as f64;
    let max_cor_numerator = statrs::function::factorial::binomial(n_m, 2)
        + (n_observations * m_value) as f64
        + statrs::function::factorial::binomial(m_value, 2);
    let max_cor = max_cor_denominator / max_cor_numerator;

    let mut out_matrix = if scale_max {
        cor_matrix / max_cor
    } else {
        cor_matrix
    };

    if diag_good {
        let n_good = exclude_loc.mapv(|x| (!x) as u64).sum_axis(Axis(0));
        let max_n_good = *n_good.into_iter().max().unwrap() as f64;
        out_matrix
            .diag_mut()
            .zip_mut_with(&n_good, |v, good_v| *v = *good_v as f64 / max_n_good);
    }

    out_matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use csv::WriterBuilder;
    use ndarray_csv::{Array2Reader, Array2Writer};
    use std::env::current_dir;
    use std::fs::File;

    #[test]
    fn test_wine_example() {
        let wine_file = File::open(format!(
            "{}/data/wine.csv",
            current_dir().unwrap().display()
        ))
        .unwrap();
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b',')
            .flexible(false)
            .from_reader(wine_file);
        let wine_dataset: Array2<Numeric> = reader.deserialize_array2((178, 14)).unwrap();
        println!("Dims: {:?}", wine_dataset.dim());
        // println!("{:#?}", wine_dataset);

        let ici_kendall_tau_mat = visqc_ici_kendall_tau(
            wine_dataset,
            true,
            true,
            true,
            0.,
            Perspective::Global,
            true,
            true,
            false,
        );
        let kendall_wine_mat_file = File::create(format!(
            "{}/result/wine_kendall_tau.csv",
            current_dir().unwrap().display()
        ))
        .unwrap();
        let mut writer = WriterBuilder::new()
            .has_headers(false)
            .from_writer(kendall_wine_mat_file);
        writer.serialize_array2(&ici_kendall_tau_mat).unwrap();
    }

    #[test]
    fn test_big_matrix2() {
        let wine_file = File::open(format!(
            "{}/data/big_matrix2.csv",
            current_dir().unwrap().display()
        ))
        .unwrap();
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b',')
            .flexible(false)
            .from_reader(wine_file);
        let big_matrix2_dataset: Array2<Numeric> = reader.deserialize_array2((875, 179)).unwrap();
        // println!("Dims: {:?}", big_matrix2_dataset.dim());
        // println!("{:#?}", wine_dataset);

        let ici_kendall_tau_mat = visqc_ici_kendall_tau(
            big_matrix2_dataset,
            true,
            true,
            true,
            0.,
            Perspective::Global,
            true,
            true,
            false,
        );
        let kendall_big_matrix2_mat_file = File::create(format!(
            "{}/result/big_matrix2_kendall_tau.csv",
            current_dir().unwrap().display()
        ))
        .unwrap();
        let mut writer = WriterBuilder::new()
            .has_headers(false)
            .from_writer(kendall_big_matrix2_mat_file);
        writer.serialize_array2(&ici_kendall_tau_mat).unwrap();
    }
}
