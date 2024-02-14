# linsolve  : System of linear equations solver

Provides solvers for systems of linear equations frequently encountered in engineering problems, such as circuit analysis, structural analysis, and fluid dynamics.


## Usage


```toml

[dependencies]
linsolve = { git = "https://github.com/hash451722/linsolve.git" }

```



```rust

use linsolve::{CooSquareMatrix, Solver};

fn main() {
    // Specify the number of dimensions of the coefficient matrix.
    let nth: usize = 5;

    // Construct a square matrix in COO format.
    let mut a: CooSquareMatrix = CooSquareMatrix::new(nth);

    // Input the non-zero value. a.push(row, col, value)
    a.push(0, 0, 45.0);
    a.push(0, 1, 6.0);
    a.push(0, 2, 7.0);
    a.push(1, 0, 10.0);
    a.push(1, 1, 20.0);
    a.push(1, 2, 23.0);
    a.push(2, 0, 15.0);
    a.push(2, 1, 50.0);
    a.push(2, 2, 67.0);
    a.push(3, 3, 27.0);
    a.push(0, 3, 0.2);
    a.push(4, 4, 17.0);

    // Define the constant vector b
    let b: Vec<f64> = vec![6.0, 16.0, 24.0, 2.5, -0.2];

    // Calculations with various solvers
    let x_lu: Vec<f64> = Solver::lu(&a, &b);
    let x_gs: Vec<f64> = Solver::gaussseidel(&a, &b, 1.0e-4, 100);
    let x_sor: Vec<f64> = Solver::sor(&a, &b, 1.1, 1.0e-4, 100);
    let x_bicg: Vec<f64> = Solver::bicg(&a, &b, 1.0e-4, 100);
    let x_bicgstab: Vec<f64> = Solver::bicgstab(&a, &b, 1.0e-4, 100);

    // Display of calculation results
    println!("{:?}", &a);
    println!("b= {:?}", &b);
    println!("---");
    println!(" x_lu       = {:?}", &x_lu);
    println!(" x_gs       = {:?}", &x_gs);
    println!(" x_sor      = {:?}", &x_sor);
    println!(" x_bicg     = {:?}", &x_bicg);
    println!(" x_bicgstab = {:?}", &x_bicgstab);
    println!("---");

    // Elimination of the coefficient matrix
    a.clear();
    println!("{:?}", &a);

}


fn main() {
    // Specify the number of dimensions of the coefficient matrix.
    let nth: usize = 5;

    // Construct a square matrix in COO format.
    let mut a: CooSquareMatrix = CooSquareMatrix::new(nth);

    // Input the non-zero value. a.push(row, col, value)
    a.push(0, 0, 45.0);
    a.push(0, 1, 6.0);
    a.push(0, 2, 7.0);
    a.push(1, 0, 10.0);
    a.push(1, 1, 20.0);
    a.push(1, 2, 23.0);
    a.push(2, 0, 15.0);
    a.push(2, 1, 50.0);
    a.push(2, 2, 67.0);
    a.push(3, 3, 27.0);
    a.push(0, 3, 0.2);
    a.push(4, 4, 17.0);

    // Define the constant vector b
    let b: Vec<f64> = vec![6.0, 16.0, 24.0, 2.5, -0.2];

    // Calculations with various solvers
    let x_lu: Vec<f64> = Solver::lu(&a, &b);
    let x_gs: Vec<f64> = Solver::gaussseidel(&a, &b, 1.0e-4, 100);
    let x_sor: Vec<f64> = Solver::sor(&a, &b, 1.1, 1.0e-4, 100);
    let x_bicg: Vec<f64> = Solver::bicg(&a, &b, 1.0e-4, 100);
    let x_bicgstab: Vec<f64> = Solver::bicgstab(&a, &b, 1.0e-4, 100);

    // Display of calculation results
    println!("{:?}", &a);
    println!("b= {:?}", &b);
    println!("---");
    println!(" x_lu       = {:?}", &x_lu);
    println!(" x_gs       = {:?}", &x_gs);
    println!(" x_sor      = {:?}", &x_sor);
    println!(" x_bicg     = {:?}", &x_bicg);
    println!(" x_bicgstab = {:?}", &x_bicgstab);
    println!("---");

    // Elimination of the coefficient matrix
    a.clear();
    println!("{:?}", &a);

}

```


## References

