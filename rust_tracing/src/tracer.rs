use crate::geometry;
use crate::brdf::Brdf;
use crate::bvh::BVH;

const MAX_RAY_DEPTH: usize = 5;

// returns a spectral radiance in [r,g,b]
pub fn trace_ray(ray_origin: [f64; 3], ray_dir: [f64; 3], bvh: &BVH) -> [f64; 3] {
    // standard path tracing 
    // what we need to keep track of at each bounce (assuming brdf sampling):
    //     - incoming spectral radiance (can compound)
    //     - (-direction of incoming light) . normal at surface
    //     - reflection spectra at each bounce
    //     - emission spectra at each bounce
    
    // the last entry is for the first bounce
    let mut angle_attenuations: [f64; MAX_RAY_DEPTH] = [0.0; MAX_RAY_DEPTH];
    let mut brdfs: [&Brdf; MAX_RAY_DEPTH] = [&Brdf {emissivity:0.0, reflection_spectrum:[0.0,0.0,0.0],inside_index_of_refraction: 1.0}; MAX_RAY_DEPTH];
    // index for angle_attenuations
    let mut i = 0;
    // for bailing early 
    let mut encountered_emissive_feature = false;

    // for working_ray_origin
    let mut w_ray_origin = ray_origin;
    let mut w_ray_dir = ray_dir;

    // first loop, bouncing
    let mut depth = 0;
    while depth < MAX_RAY_DEPTH { 
        let obj_option = bvh.intersect(&w_ray_origin, &w_ray_dir);
        match obj_option {
            None => {
                // no object was hit. The ray entered empty space, forever alone
                break;
            }
            Some((intersection, normal, brdf)) => {
                // we hit an object. Calculate the direction of the new array
                if brdf.emissivity > 0.0 {
                    encountered_emissive_feature = true;
                }
                w_ray_dir = brdf.sample(&w_ray_dir, &normal);

                w_ray_origin = intersection;

                // propagate ray slightly to avoid intersecting the same surface again 
                // w_ray_origin = geometry::add(&w_ray_origin, &geometry::scale(geometry::EPSILON, &w_ray_dir));

                angle_attenuations[i] = geometry::dot(&w_ray_dir, &normal);
                brdfs[i] = &brdf;
                i += 1; 
            }
        }
        depth += 1;
    }

    // bailing out early if this path has no spectral radiance 
    if encountered_emissive_feature == false {
        return [0.0, 0.0, 0.0];
    }

    // backtrack 
    let mut spectral_radiance = [0.0, 0.0, 0.0];
    while i >= 1 {
        i -= 1;
        // absorption 
        spectral_radiance[0] *= brdfs[i].reflection_spectrum[0];
        spectral_radiance[1] *= brdfs[i].reflection_spectrum[1];
        spectral_radiance[2] *= brdfs[i].reflection_spectrum[2];

        // attenuation
        spectral_radiance = geometry::scale(angle_attenuations[i], &spectral_radiance);

        // emission 
        spectral_radiance = geometry::add(&spectral_radiance, 
            &geometry::scale(brdfs[i].emissivity, &brdfs[i].reflection_spectrum)
        );
    }

    spectral_radiance
}