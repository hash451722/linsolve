#[derive(Debug)]
pub struct CooSquareMatrix {
    nth: usize,  // n-th order
    row_indices: Vec<usize>,
    col_indices: Vec<usize>,
    values: Vec<f64>,
}

impl CooSquareMatrix {
    pub fn new(nth: usize) -> CooSquareMatrix {
        CooSquareMatrix {
            nth,
            row_indices: Vec::new(),
            col_indices: Vec::new(),
            values: Vec::new(),
        }
    }

    pub fn push(&mut self, row: usize, col: usize, value: f64) {
        self.row_indices.push(row);
        self.col_indices.push(col);
        self.values.push(value);
    }

    pub fn clear(&mut self) {
        self.row_indices = Vec::new();
        self.col_indices = Vec::new();
        self.values = Vec::new();
    }

    fn ij(&self, i: usize, j: usize) -> f64 {
        let mut a_ij: f64 = 0.0;
        for n in 0..self.row_indices.len() {
            if self.row_indices[n] == i && self.col_indices[n] == j {
                a_ij = self.values[n];
                break;
            }
        }
        a_ij
    }

    // Multiplication
    fn mul(&self, vector: &Vec<f64>) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; vector.len()];
        for n in 0..self.row_indices.len() {
            x[self.row_indices[n]] += self.values[n] * vector[self.col_indices[n]];
        }
        x
    }

    // Traspose > Multiplication
    fn tmul(&self, vector: &Vec<f64>) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; vector.len()];
        for n in 0..self.row_indices.len() {
            x[self.col_indices[n]] += self.values[n] * vector[self.row_indices[n]];
        }
        x
    }
}


struct Vector {}
impl Vector {
    fn add(v0: &Vec<f64>, v1: &Vec<f64>) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; v0.len()]; 
        for i in 0..v0.len() {
            x[i] = v0[i] + v1[i];
        }
        x
    }

    fn sub(v0: &Vec<f64>, v1: &Vec<f64>) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; v0.len()]; 
        for i in 0..v0.len() {
            x[i] = v0[i] - v1[i];
        }
        x
    }

    fn mul(s: f64, v: &Vec<f64>) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; v.len()];
        for i in 0..v.len() {
            x[i] = s * v[i];
        }
        x
    }

    fn dot(v0: &Vec<f64>, v1: &Vec<f64>) -> f64 {
        let mut x: f64 = 0.0; 
        for i in 0..v0.len() {
            x += v0[i]*v1[i];
        }
        x
    }
}


pub struct Solver {}
impl Solver {
    fn lu_decomposition(a: &CooSquareMatrix) -> Vec<Vec<f64>> {
        let mut lu_matrix: Vec<Vec<f64>> = vec![vec![0.0; a.nth]; a.nth];
        for i in 0..a.nth {
            for j in 0..a.nth {
                if i > j {  // Lower
                    for k in 0..j {
                        lu_matrix[i][j] -= lu_matrix[i][k] * lu_matrix[k][j];
                    }
                    lu_matrix[i][j] += a.ij(i, j);
                    lu_matrix[i][j] /= lu_matrix[j][j];
                } else {  // Upper
                    for k in 0..i {
                        lu_matrix[i][j] -= lu_matrix[i][k]*lu_matrix[k][j];
                    }
                    lu_matrix[i][j] += a.ij(i, j);
                }
            }
        }
        lu_matrix
    }

    // Solve the linear system Ax = b using LU factorization
    pub fn lu(a: &CooSquareMatrix, b: &Vec<f64>) -> Vec<f64> {
        let lu_matrix: Vec<Vec<f64>> = Solver::lu_decomposition(&a);
        // Ly=b
        let mut y: Vec<f64> = vec![0.0; b.len()];
        for i in 0..b.len() {
            for j in 0..i {
                y[i] -= lu_matrix[i][j] * y[j];
            }
            y[i] += b[i];
        }
        //Ux=y
        let mut x: Vec<f64> = vec![0.0; b.len()];
        for i in (0..b.len()).rev() {
            for j in (i+1)..b.len() {
                x[i] -= lu_matrix[i][j] * x[j];
            }
            x[i] = (y[i] + x[i]) / lu_matrix[i][i];
        }
        x
    }


    // 2-norm
    fn norm2(x: &Vec<f64>) -> f64 {
        let mut s: f64 = 0.0;
        for n in 0..x.len() {
            s += (x[n].abs()).powi(2);
        }
        s.sqrt()
    }
    fn err(a: &CooSquareMatrix, x: &Vec<f64>, b: &Vec<f64>) -> f64 {
        let mut e: Vec<f64> = a.mul(&x);
        for i in 0..b.len() {
            e[i] = b[i] - e[i];
        }
        Solver::norm2(&e) / Solver::norm2(&b)
    }

    pub fn gaussseidel(a: &CooSquareMatrix, b: &Vec<f64>, tol: f64, maxiter: u32) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; b.len()];
        // println!("{}:  x_gs= {:?}", 0, &x);
        for _iter in 0..maxiter {
            for i in 0..b.len() {
                let mut s: f64 = 0.0;
                for j in 0..b.len() {
                    if i != j {
                        s -= a.ij(i, j) * x[j];
                    }
                } 
                x[i] = (b[i] + s) / a.ij(i, i);
            }
            // println!("{}:  x_gs= {:?}", iter+1, &x);
            if Solver::err(&a, &x, &b) < tol {
                break;
            }
        }
        x
    }

    pub fn sor(a: &CooSquareMatrix, b: &Vec<f64>, omega: f64, tol: f64, maxiter: u32) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; b.len()];
        // println!("{}:  x_sor= {:?}", 0, &x);
        for _iter in 0..maxiter {
            for i in 0..b.len() {
                let mut s: f64 = 0.0;
                for j in 0..b.len() {
                    if i != j {
                        s -= a.ij(i, j) * x[j];
                    }
                } 
                x[i] = (1.0 - omega) * x[i] +  omega * (b[i] + s) / a.ij(i, i);
            }
            // println!("{}:  x_sor= {:?}", iter+1, &x);
            if Solver::err(&a, &x, &b) < tol {
                break;
            }
        }
        x
    }

    // Bi-Conjugate Gradient method
    pub fn bicg(a: &CooSquareMatrix, b: &Vec<f64>, tol: f64, maxiter: u32) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; b.len()];
        let mut r: Vec<f64> = Vector::sub(b, &a.mul(&x));  // residual error
        let mut rs: Vec<f64> = r.clone();  // shadow residual error
        let mut p: Vec<f64> = r.clone();
        let mut ps: Vec<f64> = rs.clone();

        for _iter in 0..maxiter {
            let alpha: f64 = Vector::dot(&r, &rs) / Vector::dot(&a.mul(&p), &rs);
            // println!("{}:  alpha= {:?}", _iter, alpha);

            x = Vector::add(&x, &Vector::mul(alpha, &p));
            // println!("{}:  x= {:?}", _iter+1, x);

            if Solver::err(&a, &x, &b) < tol {
                break;
            }

            let mut beta: f64 = 1.0 / Vector::dot(&r, &rs);

            r = Vector::sub(&r, &Vector::mul(alpha, &a.mul(&p)));
            // println!("{}:  r= {:?}", _iter+1, r);

            rs = Vector::sub(&rs, &Vector::mul(alpha, &a.tmul(&ps)));
            // println!("{}:  r*= {:?}", _iter+1, rs);

            // println!("{}:  beta-1= {:?}", _iter, beta);
            beta *= Vector::dot(&r, &rs);
            // println!("{}:  beta= {:?}", _iter, beta);

            p = Vector::add(&r, &Vector::mul(beta, &p));
            // println!("{}:  p= {:?}", _iter+1, p);
            ps = Vector::add(&rs, &Vector::mul(beta, &ps));
            // println!("{}:  p*= {:?}", _iter+1, ps);

        }
        x
    }

    // Bi-Conjugate Gradient Stabilized method
    pub fn bicgstab (a: &CooSquareMatrix, b: &Vec<f64>, tol: f64, maxiter: u32) -> Vec<f64> {
        let mut x: Vec<f64> = vec![0.0; b.len()];
        let mut r: Vec<f64> = Vector::sub(b, &a.mul(&x));  // residual error
        let rs0: Vec<f64> = r.clone();  // shadow residual error
        let mut p: Vec<f64> = r.clone();

        for _iter in 0..maxiter {
            let alpha: f64 = Vector::dot(&r, &rs0) / Vector::dot(&a.mul(&p), &rs0);
            // println!("{}:  alpha= {:?}", _iter, alpha);

            let s: Vec<f64> = Vector::sub(&r, &Vector::mul(alpha, &a.mul(&p)));
            // println!("{}:  s= {:?}", _iter, &s);

            let omega: f64 = Vector::dot(&a.mul(&s), &s) / Vector::dot(&a.mul(&s), &a.mul(&s));
            // println!("{}:  omega= {:?}", _iter, omega);

            x = Vector::add(&x, &Vector::mul(alpha, &p));
            x = Vector::add(&x, &Vector::mul(omega, &s));
            // println!("{}:  x= {:?}", _iter+1, &x);

            if Solver::err(&a, &x, &b) < tol {
                break;
            }

            let mut beta: f64 =  Vector::dot(&r, &rs0);

            r = Vector::sub(&s, &Vector::mul(omega, &a.mul(&s)));
            // println!("{}:  r= {:?}", _iter+1, &r);

            beta = Vector::dot(&r, &rs0) / beta * alpha/omega;
            // println!("{}:  omega= {:?}", _iter, beta);

            p = Vector::sub(&p, &Vector::mul(omega, &a.mul(&p)));
            p = Vector::add(&r, &Vector::mul(beta, &p));
            // println!("{}:  p= {:?}", _iter+1, &p);
        }
        x
    }

}




// fn main() {
//     // Specify the number of dimensions of the coefficient matrix.
//     let nth: usize = 5;

//     // Construct a square matrix in COO format.
//     let mut a: CooSquareMatrix = CooSquareMatrix::new(nth);

//     // Input the non-zero value. a.push(row, col, value)
//     a.push(0, 0, 45.0);
//     a.push(0, 1, 6.0);
//     a.push(0, 2, 7.0);
//     a.push(1, 0, 10.0);
//     a.push(1, 1, 20.0);
//     a.push(1, 2, 23.0);
//     a.push(2, 0, 15.0);
//     a.push(2, 1, 50.0);
//     a.push(2, 2, 67.0);
//     a.push(3, 3, 27.0);
//     a.push(0, 3, 0.2);
//     a.push(4, 4, 17.0);

//     // Define the constant vector b
//     let b: Vec<f64> = vec![6.0, 16.0, 24.0, 2.5, -0.2];

//     // Calculations with various solvers
//     let x_lu: Vec<f64> = Solver::lu(&a, &b);
//     let x_gs: Vec<f64> = Solver::gaussseidel(&a, &b, 1.0e-4, 100);
//     let x_sor: Vec<f64> = Solver::sor(&a, &b, 1.1, 1.0e-4, 100);
//     let x_bicg: Vec<f64> = Solver::bicg(&a, &b, 1.0e-4, 100);
//     let x_bicgstab: Vec<f64> = Solver::bicgstab(&a, &b, 1.0e-4, 100);

//     // Display of calculation results
//     println!("{:?}", &a);
//     println!("b= {:?}", &b);
//     println!("---");
//     println!(" x_lu       = {:?}", &x_lu);
//     println!(" x_gs       = {:?}", &x_gs);
//     println!(" x_sor      = {:?}", &x_sor);
//     println!(" x_bicg     = {:?}", &x_bicg);
//     println!(" x_bicgstab = {:?}", &x_bicgstab);
//     println!("---");

//     // Elimination of the coefficient matrix
//     a.clear();
//     println!("{:?}", &a);

// }
