use crate::geometry;
// main idea of function: Take a set of camera information and return a set of rays that need to 
// be traced. Each ray should have a starting location, direction and be mapped to a pixel 
// Details:
//  Camera details:
//      camera_loc: the location of the center of the sensor of the camera
//      camera_dir: the orientation of a normal (outward) ray emanating from the sensor
//      camera_focal_length: the focal length of the camera lens
//      camera_sensor_to_lens: the shortest distance from the sensor to the lens
//      sensor_diagonal: we assume the sensor is square. The diagonal gives us the size of the sensor
//      aperture_radius: the radius of the lens (aperture). useful for depth of field
//      image_resolution: we assume the sensor is square of size image_resolution x image_resolution
//  Return type: Array of length (image_resolution^2 * 5) * (3 + 3)
//          = (# of rays) * (3 coords for starting position + 3 coords for direction (normalized))
//      The returned array is a flat array. Each 30 entries is one pixel, with pixels reading left to right, top down
pub fn ray_bundle(camera_loc: &[f64; 3], camera_dir: &[[f64; 3]; 3], camera_focal_length: f64, 
        camera_sensor_to_lens: f64, sensor_diagonal: f64, aperture_radius: f64, image_resolution: usize) -> Vec<f64> {
    
    let arr_length = image_resolution * image_resolution * 5 * 6;
    let mut result:Vec<f64> = Vec::with_capacity(arr_length);
    
    let sensor_width = sensor_diagonal / (2.0_f64).sqrt();

    // the location of the center of the lens 
    let lens_center = geometry::add(&camera_loc, 
        &geometry::mat_mult_3(&camera_dir, &[0.0, 0.0, camera_sensor_to_lens])
    );

    // offsets for anti-aliasing
    let offsets: [[f64;2];5] = [
        [0.0,0.0], [-2.0/3.0, -2.0/3.0], [-2.0/3.0, 2.0/3.0], 
        [2.0/3.0, 2.0/3.0], [2.0/3.0, -2.0/3.0]
    ];

    for col in 0..image_resolution {
        for row in 0..image_resolution {
            for j in 0..5 {
                // compute the actual coordinates for the sensor pixel 

                // convert col and row to somewhere within the box [-1, 1] x [-1, 1]
                let offset_x = offsets[j][0];
                let offset_y = offsets[j][1];
                // negative since image is flipped
                let sensor_loc_unit_x = -(2.0 * ((col as f64) + offset_x) / (image_resolution as f64) - 1.0);
                let sensor_loc_unit_y = (2.0 * ((row as f64) + offset_y) / (image_resolution as f64) - 1.0);

                let sensor_loc_raw = [-sensor_loc_unit_x * sensor_width / 2.0, -sensor_loc_unit_y * sensor_width / 2.0, 0.0];

                // now we to rotate the sensor location to match the orientation of the camera
                let sensor_loc: [f64; 3] = geometry::add(&camera_loc,
                    &geometry::mat_mult_3(&camera_dir, &sensor_loc_raw)
                );

                // the idea is that we can calculate the image location (using camera_focal_length), then 
                // sample the lens and return a ray pointing toward the image location
                let naive_ray_dir = geometry::normalize(&geometry::sub(&lens_center, &sensor_loc));
                let image_loc = geometry::add(&lens_center,
                    &geometry::scale(camera_focal_length, &naive_ray_dir)
                );

                // sample the lens (i.e. a circle of radius aperture_radius)
                let disk_sample = geometry::sample_disk(aperture_radius);
                let lens_pt_raw = [disk_sample[0], disk_sample[1], camera_sensor_to_lens];
                let lens_pt = geometry::add(&camera_loc,
                    &geometry::mat_mult_3(&camera_dir, &lens_pt_raw)
                );

                // create the ray 
                let ray_dir = geometry::normalize(&geometry::sub(&image_loc, &lens_pt));

                result.extend_from_slice(&lens_pt);
                result.extend_from_slice(&ray_dir);
            }
        }
    }

    result
}

pub fn spectral_radiance_to_color(spectral_radiance: &[f64; 3]) -> [u8;4] {
    fn f(x:f64) -> f64 {
        2.0 / (1.0 + (std::f64::consts::E).powf(-2.5*x)) - 1.0
    }
    let r = (f(spectral_radiance[0]) * 255.0).floor() as u8;
    let g = (f(spectral_radiance[1]) * 255.0).floor() as u8;
    let b = (f(spectral_radiance[2]) * 255.0).floor() as u8;
    let a = 255u8;

    [r,g,b,a]
}