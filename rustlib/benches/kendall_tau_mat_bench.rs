use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ndarray::Array2;
use ndarray_csv::Array2Reader;
use rustlib::kendall_tau::Perspective;
use rustlib::kendall_tau_matrix::visqc_ici_kendall_tau;
use rustlib::Numeric;
use std::env::current_dir;
use std::fs::File;
fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("ici-kendall-tau mat wine", |b| {
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
        // println!("Dims: {:?}", wine_dataset.dim());
        // println!("{:#?}", wine_dataset);

        b.iter(|| {
            visqc_ici_kendall_tau(
                black_box(wine_dataset.t().into_owned().clone()),
                true,
                true,
                true,
                0.,
                Perspective::Global,
                true,
                true,
                false,
            );
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
