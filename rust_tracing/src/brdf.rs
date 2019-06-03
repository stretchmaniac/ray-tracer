use crate::geometry;

#[derive(Clone)]
pub struct Brdf {
    // reflection in r, g and b individually, i.e. reflection_spectrum \in [0,1]^3
    pub reflection_spectrum: [f64; 3],
    // we assume the color of emission is the same as the reflection_spectrum
    pub emissivity: f64,
    pub inside_index_of_refraction: f64
}

impl Brdf {
    // the normal gives the direction to the exterior of the surface
    pub fn sample(&self, in_dir: &[f64; 3], normal: &[f64; 3]) -> [f64; 3] {
        // start with perfectly diffuse material
        // sample hemisphere
        let out_raw = geometry::sample_z_hemisphere();
        // find orthonormal basis with z axis equal to normal
        let new_x = geometry::unit_perp_vec(&normal);
        let new_y = geometry::cross(&normal, &new_x);

        geometry::add(
            &geometry::scale(out_raw[0], &new_x),
            &geometry::add(
                &geometry::scale(out_raw[1], &new_y), 
                &geometry::scale(out_raw[2], &normal)
            )
        )
    }
}