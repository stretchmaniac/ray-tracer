use rand::prelude::*;

pub const EPSILON:f64 = 0.00000001;

// 3D dot product. Returns x.y
pub fn dot(x: &[f64; 3], y: &[f64; 3]) -> f64 {
    x[0]*y[0] + x[1]*y[1] + x[2]*y[2]
}

// 3D cross product. Returns x (cross) y
pub fn cross(x: &[f64; 3], y: &[f64; 3]) -> [f64; 3] {
    [ x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0] ]
}

// addition of two 3D vectors. Returns x + y
pub fn add(x: &[f64; 3], y: &[f64; 3]) -> [f64; 3] {
    [ x[0]+y[0], x[1]+y[1], x[2]+y[2] ]
}

// subtraction of two 3D vectors. Returns x - y
pub fn sub(x: &[f64; 3], y: &[f64; 3]) -> [f64; 3] {
    [ x[0]-y[0], x[1]-y[1], x[2]-y[2] ]
}

// returns the euclidean norm of x
pub fn mag(x: &[f64; 3]) -> f64 {
    (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]).sqrt()
}

// returns the euclidean distance between two vectors
pub fn d(x: &[f64; 3], y: &[f64; 3]) -> f64 {
    mag(&sub(&x, &y))
}

// gives c * x, where x \in \mathbb{R}^3, c\in\mathbb{R}
pub fn scale(c: f64, x: &[f64; 3]) -> [f64; 3] {
    [ c*x[0], c*x[1], c*x[2] ]
}

pub fn normalize(x: &[f64; 3]) -> [f64; 3] {
    let m = mag(&x);
    if m == 0.0 {
        panic!("Division by zero error in normalize. Input vector: [{}, {}, {}]", x[0], x[1], x[2]);
    }
    scale(1.0/m, &x)
}

pub fn min(x:f64, y:f64) -> f64 {
    if x < y {x} else {y}
}
pub fn max(x:f64, y:f64) -> f64 {
    if x < y {y} else {x}
}

// returns a random point in the unit disk
pub fn sample_disk(radius: f64) -> [f64; 2] {
    // we sample the rectangle [0,1] x [0, 2pi], sampling 
    // [0, 2pi] uniformly but [0,1] with x \in [0, 1] \mapsto sqrt(x)
    let mut rng = rand::thread_rng();
    let r_raw: f64 = rng.gen();
    let r = radius * r_raw.sqrt();
    let theta_raw: f64 = rng.gen();
    let theta: f64 = std::f64::consts::PI * 2.0 * theta_raw;
    [r * theta.cos(), r * theta.sin()]
}

// returns (A x^T)^T, where T is the transpose operator
pub fn mat_mult_3(A: &[[f64; 3]; 3] , x: &[f64; 3]) -> [f64; 3] {
    [
        A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2], 
        A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2],
        A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2]
    ]
}

pub fn poly_inclusion_2d(points: &Vec<[f64; 2]>, test: &[f64; 2]) -> bool {
    // typical algorithm
    let l = points.len();
    let mut crossings = 0;
    for i in 0..l {
        let p1 = points[i];
        let mut p2 = points[0];
        if i < l - 1 {
            p2 = points[i + 1];
        }
        if p1[1] > test[1] && p2[1] < test[1] {
            let d = p1[1] - p2[1];
            let w1 = d - (p1[1] - test[1]);
            let w2 = d - (test[1] - p2[1]);
            let x = (w1 * p1[0] + w2 * p2[0]) / d;
            if x >= test[0] {
                crossings += 1;
            }
        }
        if p1[1] < test[1] && p2[1] > test[1] {
            let d = p2[1] - p1[1];
            let w1 = d - (test[1] - p1[1]);
            let w2 = d - (p2[1] - test[1]);
            let x = (w1 * p1[0] + w2 * p2[0]) / d;
            if x >= test[0] {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

// uniform sample of the unit hemisphere
pub fn sample_z_hemisphere() -> [f64; 3] {
    // i.e. so that z >= 0
    let mut rng = rand::thread_rng();
    let theta_raw: f64 = rng.gen();
    let theta: f64 = std::f64::consts::PI * 2.0 * theta_raw;
    let phi_raw: f64 = rng.gen();
    let phi: f64 = phi_raw.acos();

    [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()] 
}

// returns a unit vector orthogonal to x of length 1
pub fn unit_perp_vec(x: &[f64; 3]) -> [f64; 3] {
    // one of [1,0,0] and [0,1,0] is more than 45 degrees from x 
    let mut v = [1.0, 0.0, 0.0]; 
    let other = [0.0, 1.0, 0.0];
    if dot(&x, &v).abs() > dot(&x, &other).abs() {
        v = other;
    }

    return normalize(&cross(&x, &v));
}

#[cfg(test)]
mod tests {
    use crate::geometry;
    const EPSILON:f64 = 0.00000001;

    #[test]
    fn poly_inclusion_2d_test() {
        let pts = vec![[0.0,0.0],[1.0,0.0],[1.0,1.0]];
        assert!(geometry::poly_inclusion_2d(&pts, &[0.5,0.75]) == false);
        assert!(geometry::poly_inclusion_2d(&pts, &[0.5,0.25]) == true);
    }

    #[test]
    fn mat_mult_3_test(){
        // test a few linear transformations 
        // identity
        let A = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let x1 = [1.0, 0.0, 0.0];
        let x2 = [1.0, 1.0, -1.0];

        let a_x1 = geometry::mat_mult_3(&A, &x1);
        let a_x2 = geometry::mat_mult_3(&A, &x2);
        assert!( (geometry::d(&a_x1, &x1)).abs() < EPSILON );
        assert!( (geometry::d(&a_x2, &x2)).abs() < EPSILON );

        // flip y and z coordinates 
        let A2 = [[1.0,0.0,0.0],[0.0,0.0,1.0],[0.0,1.0,0.0]];
        let a2_x1 = geometry::mat_mult_3(&A2, &x1);
        let a2_x2 = geometry::mat_mult_3(&A2, &x2);
        assert!( (geometry::d(&a2_x1, &x1)).abs() < EPSILON );
        assert!( (geometry::d(&a2_x2, &[1.0,-1.0,1.0])).abs() < EPSILON );

        // random transformation 
        let A3 = [[-1.4, 2.5, -3.4], [0.9, 2.1, 4.4], [1.5,-3.2, 1.2]];
        let x3 = [1.2, -6.7, 0.1];
        let a3_x3 = geometry::mat_mult_3(&A3, &x3);
        assert!( (geometry::d(&a3_x3, 
            &[geometry::dot(&A3[0], &x3), geometry::dot(&A3[1], &x3), geometry::dot(&A3[2], &x3)]
        )).abs() < EPSILON );
    }

    #[test]
    fn dot_test() {
        let v1 = [1.0, 2.0, 3.0];
        let v2 = [2.0, -3.0, 0.0];
        // expected: 1*2 - 3*2 = -4
        assert!( (geometry::dot(&v1, &v2) + 4.0).abs() < EPSILON );
        
        let v3 = [0.0, 0.0, 0.0];
        // expected: 0
        assert!( (geometry::dot(&v1, &v3)).abs() < EPSILON );
    }

    #[test]
    fn sub_test() {
        let v1 = [1.0, 2.0, 3.0];
        let v2 = [2.0, -5.0, 1.2];
        let v3 = [0.0, -2.1, 12.0];

        // note that we usually check if two vectors are equal by 
        // computing mag(sub(x,y)) < EPSILON. Can't exactly do that for the 
        // sub function, now can we?

        // expected: [-1, 7, 1.8]
        let sub_v1_v2 = geometry::sub(&v1, &v2);
        assert!( (sub_v1_v2[0] - -1.0).abs() < EPSILON );
        assert!( (sub_v1_v2[1] - 7.0).abs() < EPSILON );
        assert!( (sub_v1_v2[2] - 1.8).abs() < EPSILON );

        // expected: [1, 4.1, -9]
        let sub_v1_v3 = geometry::sub(&v1, &v3);
        assert!( (sub_v1_v3[0] - 1.0).abs() < EPSILON );
        assert!( (sub_v1_v3[1] - 4.1).abs() < EPSILON );
        assert!( (sub_v1_v3[2] - -9.0).abs() < EPSILON );

        // expected: [2, -2.9, -10.8]
        let sub_v2_v3 = geometry::sub(&v2, &v3);
        assert!( (sub_v2_v3[0] - 2.0).abs() < EPSILON );
        assert!( (sub_v2_v3[1] - -2.9).abs() < EPSILON );
        assert!( (sub_v2_v3[2] - -10.8).abs() < EPSILON );
    }

    #[test]
    fn mag_test() {
        let v1 = [-1.0, 2.0, 3.0];
        let v2 = [1.234, 3.4561, -9.1234];
        let v3 = [0.0, 0.0, 0.0];
        let v4 = [-2.0, -2.0, -2.0];

        assert!( (geometry::mag(&v1) - 3.74165738677).abs() < EPSILON );
        assert!( (geometry::mag(&v2) - 9.83380957564).abs() < EPSILON );
        assert!( (geometry::mag(&v3) - 0.0).abs() < EPSILON );
        assert!( (geometry::mag(&v4) - 3.46410161514).abs() < EPSILON );
    }

    #[test]
    fn cross_test() {
        let v1 = [1.0, 2.2, 3.0];
        let v2 = [0.3, 1.4, 5.3];
        let v3 = [-4.3, 6.3, -0.1];

        let expected_v1v2 = [7.46, -4.4, 0.74];
        let actual_v1v2 = geometry::cross(&v1, &v2);

        let expected_v1v3 = [-19.12, -12.8, 15.76];
        let actual_v1v3 = geometry::cross(&v1, &v3);

        let expected_v2v3 = [-33.53, -22.76, 7.91];
        let actual_v2v3 = geometry::cross(&v2, &v3);

        assert!( geometry::d(&expected_v1v2, &actual_v1v2) < EPSILON );
        assert!( geometry::d(&expected_v1v3, &actual_v1v3) < EPSILON );
        assert!( geometry::d(&expected_v2v3, &actual_v2v3) < EPSILON );
    }

    #[test]
    fn d_test() {
        let v1 = [1.0, 2.0, 3.0];
        let v2 = [0.0, 1.4, 5.3];
        let v3 = [-4.3, 6.3, -0.1];

        assert!( (geometry::d(&v1, &v2) - 2.578759391645525).abs() < EPSILON ); 
        assert!( (geometry::d(&v1, &v3) - 7.495998932764064).abs() < EPSILON );
        assert!( (geometry::d(&v2, &v3) - 8.465222974027323).abs() < EPSILON );
    }

    #[test]
    fn add_test() {
        let v1 = [4.256, 0.236, -1.234];
        let v2 = [1.235, 133.463, 1.23];
        let v3 = [-1.23, -4.46534, -9.91635];
        
        let expected_v1v2 = [5.491, 133.699, -0.004];
        let expected_v1v3 = [3.026, -4.22934, -11.15035];
        let expected_v2v3 = [0.005, 128.99766, -8.68635];

        let actual_v1v2 = geometry::add(&v1, &v2);
        let actual_v1v3 = geometry::add(&v1, &v3);
        let actual_v2v3 = geometry::add(&v2, &v3);

        assert!( geometry::d(&expected_v1v2, &actual_v1v2) < EPSILON );
        assert!( geometry::d(&expected_v1v3, &actual_v1v3) < EPSILON );
        assert!( geometry::d(&expected_v2v3, &actual_v2v3) < EPSILON );
    }
}