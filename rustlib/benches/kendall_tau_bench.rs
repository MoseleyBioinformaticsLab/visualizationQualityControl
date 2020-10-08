use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ndarray::Array1;
use rand::prelude::*;
use rand_distr::StandardNormal;
use rustlib::kendall_tau::{ici_kendall_tau, Perspective};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("ici-kendall-tau 1000", |b| {
        let mut rng = StdRng::seed_from_u64(20200929);
        let n = 1000;
        let x: Array1<f64> = StandardNormal.sample_iter(&mut rng).take(n).collect();
        let y: Array1<f64> = StandardNormal.sample_iter(rng).take(n).collect();
        // println!("{}", ici_kendall_tau(x, y, Perspective::Local));

        b.iter(|| ici_kendall_tau(black_box(x.view()), black_box(y.view()), Perspective::Local));
        // b.iter(|| ici_kendall_tau2(black_box(x.view()), black_box(y.view()), Perspective::Local));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
