mod geometry;
mod camera;
mod primitives;
mod brdf;
mod tracer;
mod bvh;

extern crate image;

use std::path::Path;
use std::sync::mpsc;
use std::thread;
use std::time::SystemTime;
use std::fs;

const RESOLUTION:usize = 800;
const NUM_THREADS: usize = 7;

fn main() {
    let start_time = SystemTime::now();

    let rays_original = camera::ray_bundle(
        &[0.0, -1.7, -0.65],
        &[[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        2.0, // focal length
        0.039, // sensor/lens dist
        0.05, // sensor diagonal
        0.000001, // aperture
        RESOLUTION
    );

    let samples_per_pixel = 1950;

    let (tx, rx) = mpsc::channel();

    let mut children = Vec::new();

    // 8 threads total
    for thread_index in 0..NUM_THREADS {
        // copy the necessary data:
        let tx_thread = mpsc::Sender::clone(&tx);
        let rays_pointer = Box::new(rays_original.clone());

        let child = thread::spawn(move || {
            let mut colors = vec![0u8; RESOLUTION*RESOLUTION*4];
            let env = init_environment();

            let rays = *rays_pointer;
            println!("Thread started");
            for y in 0..RESOLUTION {
                if y % 10 == 0 && thread_index == 0 {
                    println!("Completed: {}%", (100*y)/RESOLUTION);
                }
                for x in 0..RESOLUTION {
                    if (y*RESOLUTION + x) % NUM_THREADS == thread_index {
                        let mut total_spectral_radiance = [0.0, 0.0, 0.0];
                        let start_index = (y*RESOLUTION + x) * 30;
                        for sample in 0..samples_per_pixel {
                            for j in 0..5 {
                                let sub_pix_data = &rays[start_index + 6*j..start_index + 6*j + 6];
                                let ray_origin = [sub_pix_data[0], sub_pix_data[1], sub_pix_data[2]];
                                let ray_dir = [sub_pix_data[3], sub_pix_data[4], sub_pix_data[5]];

                                let spectral_radiance = tracer::trace_ray(ray_origin, ray_dir, &env);
                                total_spectral_radiance = geometry::add(&total_spectral_radiance, &spectral_radiance);
                            }
                        }
                        total_spectral_radiance = geometry::scale(1.0/(5.0 * (samples_per_pixel as f64)), &total_spectral_radiance);
                        let color = camera::spectral_radiance_to_color(&total_spectral_radiance);
                        let color_index = (y*RESOLUTION + x) * 4;
                        colors[color_index] = color[0];
                        colors[color_index + 1] = color[1];
                        colors[color_index + 2] = color[2];
                        colors[color_index + 3] = color[3];
                    }
                }
            }

            tx_thread.send(Box::new(colors));
        });

        children.push(child);
    }

    let mut threads_finished = 0;
    let mut master_colors: [u8; RESOLUTION*RESOLUTION*4] = [0u8; RESOLUTION*RESOLUTION*4]; 
    for colors_message_box in rx {
        let colors_message = *colors_message_box;
        println!("message received");
        // collect thread results here
        for i in 0..(RESOLUTION*RESOLUTION*4) {
            if colors_message[i] > 0 {
                master_colors[i] = colors_message[i];
            }
        }
        threads_finished += 1;
        if threads_finished == NUM_THREADS {
            break;
        }
    }

    for child in children {
        child.join();
    }

    write_to_png(&master_colors);
    match start_time.elapsed() {
        Ok(elapsed) => {
            println!("runtime: {}", elapsed.as_millis());
        },
        Err(e) => {
            println!("error");
        }
    }
}

fn init_environment() -> bvh::BVH {
    const z: f64 = -1.0;
    const x: f64 = 2.0;
    const y: f64 = 2.0;
    let floor: primitives::Polygon = primitives::Polygon::new_cached(
        vec![[x,y,z],[-x,y,z],[-x,-y,z],[x,-y,z]],
        brdf::Brdf {
            reflection_spectrum: [0.75, 0.95, 0.5],
            emissivity: 0.0,
            inside_index_of_refraction: 1.0
        },
        [0.0, 0.0, 1.0] // orientation direction
    );

    const h: f64 = 0.5;
    let wall = primitives::Polygon::new_cached(
        vec![[-x-1.0,y,z+h],[-x-1.0,-y,z+h],[-x-1.0,-y,-z+h],[-x-1.0,y,-z+h]],
        brdf::Brdf {
            reflection_spectrum: [1.0,1.0,1.0],
            emissivity: 25.0,
            inside_index_of_refraction: 1.0
        },
        [1.0, 0.0, 0.0]
    );

    let wall2 = primitives::Polygon::new_cached(
        vec![[x,y,z],[x,-y,z],[x,-y,z+3.0],[x,y,z+3.0]],
        brdf::Brdf {
            reflection_spectrum: [0.95,0.95,0.95],
            emissivity: 0.0,
            inside_index_of_refraction: 1.0
        },
        [-1.0, 0.0, 0.0]
    );

    // back wall
    let wall3 = primitives::Polygon::new_cached(
        vec![[x,y,z],[-x,y,z],[-x,y,-z+1.0],[x,y,-z+1.0]],
        brdf::Brdf {
            reflection_spectrum: [0.95, 0.95, 0.95],
            emissivity: 0.0,
            inside_index_of_refraction: 1.0
        },
        [0.0, -1.0, 0.0]
    );

    
    let mut objects: Vec<Box<primitives::Primitive>> = vec![Box::new(floor), Box::new(wall), Box::new(wall2), Box::new(wall3)];
    for obj in polys_from_obj("src/hippo.obj") {
        objects.push(obj);
    }

    bvh::BVH::from_list(objects)
}

fn polys_from_obj(file_name: &str) -> Vec<Box<primitives::Primitive>> {
    let mut objects: Vec<Box<primitives::Primitive>> = Vec::new();

    let contents = fs::read_to_string(file_name).expect("something went wrong with reading the file");

    let mut vertices = Vec::new();
    let mut normals = Vec::new();

    for line in contents.split("\n") {
        if line.starts_with("v ") || line.starts_with("vn ") {
            // a vertex
            let mut coordinates = [0.0; 3];
            let mut count = 0;
            let mut line_better = line.replace("  ", " ");
            line_better = line_better.trim().to_string(); 
            for item in line_better.split(" ") {
                if count > 0{
                    coordinates[count - 1] = item.parse::<f64>().unwrap();
                }
                count += 1;
            }
            coordinates[2] -= 0.717;
            if line.starts_with("v "){
                vertices.push(coordinates);
            }else {
                normals.push(coordinates);
            }
            
        } else if line.starts_with("f "){
            // a face
            let mut vertex_indices = Vec::new();
            let mut normal_indices = Vec::new();
            let mut count = 0;
            let mut line_better = line.replace("  ", " ");
            line_better = line_better.trim().to_string();
            for item in line_better.split(" ") {
                if count > 0 {
                    let sub_items: Vec<&str> = item.split("/").collect();
                    vertex_indices.push(sub_items[0].parse::<usize>().unwrap());
                    normal_indices.push(sub_items[2].parse::<usize>().unwrap());
                }
                count += 1;
            }
            // one indexed
            let mut pts = Vec::new();
            for index in &vertex_indices {
                pts.push(vertices[index - 1].clone());
            }
            
            // vertices are in counter-clockwise order
            let dir = geometry::cross(
                &geometry::sub(&pts[1], &pts[0]),
                &geometry::sub(&pts[2], &pts[1])
            );

            let face = primitives::Polygon::new_cached(
                pts,
                brdf::Brdf {
                    reflection_spectrum: [0.65, 0.65, 0.65],
                    emissivity: 0.0,
                    inside_index_of_refraction: 1.0
                },
                dir
            );
            objects.push(Box::new(face));
        }
    }
    objects
}

fn write_to_png(colors: &[u8; RESOLUTION*RESOLUTION*4]) {
    image::save_buffer(&Path::new("image.png"), colors, RESOLUTION as u32, RESOLUTION as u32, image::RGBA(8));
}
