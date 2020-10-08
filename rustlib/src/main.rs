use ndarray::Array2;
use ndarray_csv::Array2Reader;
use rustlib::kendall_tau::Perspective;
use rustlib::kendall_tau_matrix::visqc_ici_kendall_tau;
use rustlib::Numeric;
use std::env::current_dir;
use std::fs::File;

fn main() {
    // Benchmarking the ici-kendall tau matrix algorithm.
    // We cannot benchmark the concurrent version, thus this must be shut-off somehow

    #[cfg(feature = "rayon")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();

    let big_matrix2_file = File::open(format!(
        "{}/data/big_matrix2.csv",
        current_dir().unwrap().display()
    ))
    .unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b',')
        .flexible(false)
        .from_reader(big_matrix2_file);
    let big_matrix2_dataset: Array2<Numeric> = reader.deserialize_array2((875, 179)).unwrap();
    // println!("Dims: {:?}", big_matrix2_dataset.dim());
    // println!("{:#?}", wine_dataset);

    let times = 5;
    for _ in 0..times {
        let ici_kendall_tau_mat = visqc_ici_kendall_tau(
            big_matrix2_dataset.clone(),
            true,
            true,
            true,
            0.,
            Perspective::Global,
            true,
            true,
            false,
        );
    }

    // Write the result to disk.
    // let kendall_big_matrix2_mat_file = File::create(format!(
    //     "{}/result/big_matrix2_kendall_tau.csv",
    //     current_dir().unwrap().display()
    // ))
    // .unwrap();
    // let mut writer = WriterBuilder::new()
    //     .has_headers(false)
    //     .from_writer(kendall_big_matrix2_mat_file);
    // writer.serialize_array2(&ici_kendall_tau_mat).unwrap();
}
