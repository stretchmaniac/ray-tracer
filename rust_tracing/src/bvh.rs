// bvh - bounding volume hierarchy. Used for faster scene intersection tests in 
// replacement of FeatureCollection

use crate::primitives::Primitive;
use crate::geometry;
use crate::geometry::{min, max};
use crate::brdf::Brdf;

#[derive(Clone)]
pub struct BVH {
    // list of primitives (for one node maybe a few spheres or polygons). Empty if not a leaf node
    base_contents: Vec<Box<Primitive>>,
    // center x, center y, center z, radius
    bounding_shell: [f64; 4],
    // start with binary tree
    left: Option<Box<BVH>>,
    right: Option<Box<BVH>>
}

impl BVH {
    pub fn from_list(objects: Vec<Box<Primitive>>) -> Self {
        // this can be recursive; it is done offline before the tracing begins
        if objects.len() <= 4 { // seems like a decent choice? we'll do some experiments maybe
            return BVH {
                // note that bounding_box must precede base_contents so that objects doesn't get moved
                bounding_shell: BVH::package(&objects),
                base_contents: objects,
                left: None,
                right: None
            };
        }
        
        let bounding_shell = BVH::package(&objects);
        let meanPos = [bounding_shell[0], bounding_shell[1], bounding_shell[2]];

        let mut mean_dir = [0.0, 0.0, 0.0];
        for obj in &objects {
            let prim_ref= &**obj;
            let b = prim_ref.package();
            let middle = [b[0], b[1], b[2]];
            let mut obj_dir = geometry::sub(&middle, &meanPos);
            if geometry::mag(&obj_dir) > geometry::EPSILON {
                obj_dir = geometry::normalize(&obj_dir);
            } else {
                obj_dir = [0.0, 0.0, 0.01];
            }
            if geometry::dot(&obj_dir, &mean_dir) < 0.0 {
                obj_dir = geometry::scale(-1.0, &obj_dir);
            }
            mean_dir = geometry::add(&mean_dir, &obj_dir);
        }

        let mut left_objects = Vec::new();
        let mut right_objects = Vec::new();

        for obj in objects {
            let b = (*obj).package();
            let middle = [b[0], b[1], b[2]];
            if geometry::dot(&geometry::sub(&middle, &meanPos), &mean_dir) > 0.0 {
                left_objects.push(obj);
            }else {
                right_objects.push(obj);
            }
        }

        let left_bvh = BVH::from_list(left_objects);
        let right_bvh = BVH::from_list(right_objects);

        return BVH {
            bounding_shell: bounding_shell,
            base_contents: Vec::new(),
            left: Some(Box::new(left_bvh)),
            right: Some(Box::new(right_bvh))
        }
    }

    pub fn package(objects: &[Box<Primitive>]) -> [f64; 4] {
        if objects.len() == 0 {
            return [0.0, 0.0, 0.0, 0.0];
        }
        let mut mean_pos = [0.0, 0.0, 0.0];
        for box_prim in objects {
            let b = (*box_prim).package();
            mean_pos = geometry::add(&mean_pos, &[b[0], b[1], b[2]]);
        }
        mean_pos = geometry::scale(1.0 / (objects.len() as f64), &mean_pos);

        let mut radius = 0.0;
        for box_prim in objects {
            let b = (*box_prim).package();
            radius = geometry::max(radius, geometry::d(&[b[0], b[1], b[2]], &mean_pos) + b[3]);
        }

        [mean_pos[0], mean_pos[1], mean_pos[2], radius]
    }
    
    pub fn intersect(&self, ray_origin: &[f64; 3], ray_dir: &[f64; 3]) -> Option<([f64; 3], [f64; 3], &Brdf)> {
        // we're looking for speed here, so unfortunately no recursion
        let mut bvhStack = vec![self];
        let mut closest: Option<([f64; 3], [f64; 3], &Brdf)> = None;
        let mut closest_distance = 0.0;

        while bvhStack.len() > 0 {
            let root = bvhStack.pop().unwrap();
            // check if minimum distance to root is farther than closest_distance 
            let b = root.bounding_shell;
            if closest.is_some() && geometry::d(&[b[0],b[1],b[2]], &closest.unwrap().0) - b[3] > closest_distance {
                continue;
            }
            if shell_intersect(&root.bounding_shell, &ray_origin, &ray_dir) {
                for primitive in &root.base_contents {
                    let intersection = primitive.intersect(ray_origin, ray_dir);
                    match intersection {
                        None => {},
                        Some(pos) => {
                            let d = geometry::d(&pos, &ray_origin);
                            let res = Some((pos, primitive.normal_at(&pos), primitive.brdf()));
                            match closest {
                                None => {
                                    closest = res;
                                    closest_distance = d;
                                },
                                Some((closest_pos, normal, brdf)) => {
                                    if d < closest_distance {
                                        closest = res;
                                        closest_distance = d;
                                    }
                                }
                            }
                        }
                    }
                }
                match &root.left {
                    None => {},
                    Some(bvh) => {
                        bvhStack.push(&*bvh);
                    }
                }
                match &root.right {
                    None => {},
                    Some(bvh) => {
                        bvhStack.push(&*bvh);
                    }
                }
            }
        }

        closest
    }
}

fn shell_intersect(shell: &[f64; 4], ray_origin: &[f64; 3], ray_dir: &[f64; 3]) -> bool {
    // see sphere intersection in primitives.rs
    let leg1 = [shell[0] - ray_origin[0], shell[1] - ray_origin[1], shell[2] - ray_origin[2]];
    if geometry::dot(&leg1, &ray_dir) < 0.0 && geometry::mag(&leg1) > shell[3] {
        return false;
    }
    geometry::mag(
        &geometry::cross(
            &ray_dir,
            &leg1
        )
    ) <= shell[3]
}

fn box_intersect(box_coords: &[f64; 6], ray_origin: &[f64; 3], ray_dir: &[f64; 3]) -> bool {
    // must be blazing fast
    let min_x = box_coords[0];
    let max_x = box_coords[1];
    let min_y = box_coords[2];
    let max_y = box_coords[3];
    let min_z = box_coords[4];
    let max_z = box_coords[5];
    if ray_dir[0].abs() >= geometry::EPSILON {
        let mut y1 = ray_origin[1];
        let mut z1 = ray_origin[2];
        if (ray_origin[0] > max_x && ray_dir[0] > 0.0) || (ray_origin[0] < min_x && ray_dir[0] < 0.0) {
            return false;
        }
        if ray_origin[0] < min_x || ray_origin[0] > max_x {
            // we deduce ray_dir[0] > 0 when ray_origin[0] < min_x and ray_dir[0] < 0 in the other case
            let a = if ray_origin[0] < min_x {min_x} else {max_x};
            let t0 = (a - ray_origin[0]) / ray_dir[0];
            y1 += ray_dir[1] * t0;
            z1 += ray_dir[2] * t0;
        }
        let b = if ray_dir[0] < 0.0 {max_x} else {min_x};
        let t = (b - ray_origin[0]) / ray_dir[0];
        let y2 = ray_origin[1] + ray_dir[1] * t;
        let z2 = ray_origin[2] + ray_dir[2] * t;
        return !( (y1 < min_y && y2 < min_y) || (y1 > max_y && y2 > max_y) || (z1 < min_z && z2 < min_z) || (z1 > max_z && z2 > max_z));
    }
    // ray_dir[0] approximately 0
    if ray_origin[0] < min_x || ray_origin[0] > max_x {
        return false;
    }
    // so x component checks out. Do 2D box intersection with y and z coords
    if ray_dir[1].abs() >= geometry::EPSILON {
        let mut z1 = ray_origin[2];
        if (ray_origin[1] > max_y && ray_dir[1] > 0.0) || (ray_origin[1] < min_y && ray_dir[1] < 0.0) {
            return false;
        }
        if ray_origin[1] < min_y || ray_origin[1] > max_y {
            let a = if ray_origin[1] < min_y {min_y} else {max_y};
            let t0 = (a - ray_origin[1]) / ray_dir[1];
            z1 += ray_dir[2] * t0;
        }
        let b = if ray_dir[1] < 0.0 {max_y} else {min_y};
        let t = (b - ray_origin[1]) / ray_dir[1];
        let z2 = ray_origin[2] + ray_dir[2] * t;
        return !( (z1 < min_z && z2 < min_z) || (z1 < max_z && z2 > max_z) );
    }
    // both ray_dir[0] and ray_dir[1] are approximately 0 
    if ray_origin[1] < min_y || ray_origin[1] > max_y {
        return false;
    } 
    return (ray_origin[2] < max_z && ray_dir[2] > 0.0) || (ray_origin[2] > min_z && ray_dir[2] < 0.0);

}