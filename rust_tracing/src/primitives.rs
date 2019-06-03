use crate::geometry;
use crate::geometry::{min, max};
use crate::brdf::Brdf;

pub trait Primitive : PrimitiveClone{
    // a ray consists of a starting position (first 3 slots)
    // and a direction (second 3 slots). This returns the position of the 
    // intersection
    fn intersect(&self, ray_origin: &[f64; 3], ray_dir: &[f64; 3]) -> Option<[f64; 3]>;
    fn normal_at(&self, point: &[f64; 3]) -> [f64; 3];
    // returns a bounding box of the primitive. Ideally the smallest one
    fn package(&self) -> [f64; 4];
    fn brdf(&self) -> &Brdf;
}

// see https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object
trait PrimitiveClone {
    fn clone_box(&self) -> Box<Primitive>;
}

impl<T> PrimitiveClone for T where T: 'static + Primitive + Clone, {
    fn clone_box(&self) -> Box<Primitive> {
        Box::new(self.clone())
    }
}
impl Clone for Box<Primitive> {
    fn clone(&self) -> Box<Primitive> {
        self.clone_box()
    }
}


// a sphere
#[derive(Clone)]
pub struct Sphere {
    pub center: [f64; 3],
    pub radius: f64,
    pub brdf: Brdf
}

impl Primitive for Sphere {
    fn brdf(&self) -> &Brdf {
        &self.brdf
    }
    fn intersect(&self, ray_origin: &[f64; 3], ray_dir: &[f64; 3]) -> Option<[f64; 3]> {
        // line sphere intersection. Speed is critical here
        // check shortest distance to line formed by ray. 9/10 it will be 
        // greater than the radius of the sphere and we can return immediately

        // get to work with that vector arithmetic!
        let leg1 = geometry::sub(&self.center, &ray_origin);

        // heading away from the sphere 
        if geometry::dot(&leg1, &ray_dir) < 0.0 {
            return None;
        }

        let par_area = geometry::mag(
            &geometry::cross(
                &ray_dir,
                &leg1
            )
        );

        // Ray direction is already normalized! (see camera ray_bundle)
        let r = par_area;

        if r > self.radius {
            return None;
        }

        // ok, now we have real work to do
        // sorta. draw it, it's not difficult
        let a = geometry::mag(&leg1);
        let d = (a*a - r*r).sqrt() - (self.radius*self.radius - r*r).sqrt();
        let result = geometry::add(
            &ray_origin,
            &geometry::scale(d, &ray_dir)
        );
        return Some(result);
    }
    fn normal_at(&self, point:&[f64;3]) -> [f64;3] {
        geometry::scale(1.0/self.radius, &geometry::sub(&point, &self.center))
    }
    fn package(&self) -> [f64;4] {
        [
            self.center[0], 
            self.center[1],
            self.center[2],
            self.radius
        ]
    }
}

// a 2D convex polygon embedded in 3D space
// the intersection routine does not require a convex shape, but sampling light sources requires 
// triangularization
#[derive(Clone)]
pub struct Polygon {
    pub points: Vec<[f64; 3]>,
    pub brdf: Brdf,
    // for quick intersection checks
    // doesn't have to be precisely the center
    pub bounding_center: [f64; 3],
    pub bounding_radius: f64,
    // cached info
    // everything must be normalized
    pub normal: [f64; 3],
    // orthonormal basis for plane, using bounding_center as origin
    basis_1: [f64; 3],
    basis_2: [f64; 3],
    // points written in the basis written above
    points_in_basis: Vec<[f64; 2]>
}

impl Polygon {
    pub fn parallelpiped(origin:[f64;3],gen1:[f64;3],gen2:[f64;3],gen3:[f64;3], brdf: Brdf) -> Vec<Polygon> {
        // there are 6 sets of faces
        let face1_pts = vec![origin.clone(),geometry::add(&origin, &gen1), geometry::add(&geometry::add(&origin, &gen1), &gen2), geometry::add(&origin, &gen2)];
        let face1_dir = geometry::scale(-1.0, &gen3);
        let face2_pts = vec![origin.clone(),geometry::add(&origin, &gen2), geometry::add(&geometry::add(&origin, &gen2), &gen3), geometry::add(&origin, &gen3)];
        let face2_dir = geometry::scale(-1.0, &gen1);
        let face3_pts = vec![origin.clone(),geometry::add(&origin, &gen1), geometry::add(&geometry::add(&origin, &gen1), &gen3), geometry::add(&origin, &gen3)];
        let face3_dir = geometry::scale(-1.0, &gen2);

        let face4_pts = face1_pts.iter().map(|&x| geometry::add(&x, &gen3)).collect();
        let face4_dir = gen3.clone();
        let face5_pts = face2_pts.iter().map(|&x| geometry::add(&x, &gen1)).collect();
        let face5_dir = gen1.clone();
        let face6_pts = face3_pts.iter().map(|&x| geometry::add(&x, &gen2)).collect();
        let face6_dir = gen2.clone();

        vec![
            Polygon::new_cached(face1_pts, brdf.clone(), face1_dir),
            Polygon::new_cached(face2_pts, brdf.clone(), face2_dir),
            Polygon::new_cached(face3_pts, brdf.clone(), face3_dir),
            Polygon::new_cached(face4_pts, brdf.clone(), face4_dir),
            Polygon::new_cached(face5_pts, brdf.clone(), face5_dir),
            Polygon::new_cached(face6_pts, brdf.clone(), face6_dir)
        ]
    }

    // direc_vec doesn't have to be normal, just gives an orientation for the surface
    pub fn new_cached(points: Vec<[f64; 3]>, brdf: Brdf, direc_vec: [f64; 3]) -> Polygon {
        // compute all this nonsense

        let mut normal = geometry::cross(
            &geometry::sub(&points[1], &points[0]), 
            &geometry::sub(&points[2], &points[1])
        );

        if geometry::dot(&normal, &direc_vec) < 0.0 {
            normal = geometry::scale(-1.0, &normal);
        }
        if geometry::mag(&normal) < geometry::EPSILON {
            println!("{:?}", points);
        }
        normal = geometry::normalize(&normal);

        let mut basis_1 = geometry::sub(&points[1], &points[0]);
        let mut basis_2 = geometry::cross(&normal, &basis_1);
        basis_1 = geometry::normalize(&basis_1);
        basis_2 = geometry::normalize(&basis_2);

        let mut points_in_basis = Vec::new();
        for pt in points.iter() {
            let offset = geometry::sub(pt, &points[0]);
            points_in_basis.push(
                [
                    geometry::dot(&offset, &basis_1),
                    geometry::dot(&offset, &basis_2)
                ]
            );
        }

        // compute bounding center and bounding radius. 
        // technically we should be solving the smallest circumcircle problem here. 
        // however, most of our polygons will be very small so it probably doesn't matter 
        // that much. If you're looking for performance gains, here is a possible spot

        let mut centroid = [0.0, 0.0, 0.0];
        for p in &points {
            centroid = geometry::add(&centroid, &p);
        }
        centroid = geometry::scale(1.0/(points.len() as f64), &centroid);
        let mut radius = geometry::d(&centroid, &points[0]);
        for p in &points {
            let r = geometry::d(&centroid, &p);
            if r > radius {
                radius = r;
            }
        }

        Polygon {
            points:points, brdf: brdf, bounding_center:centroid,
            bounding_radius: radius, normal: normal, 
            basis_1: basis_1, basis_2: basis_2, points_in_basis: points_in_basis
        }
    }
}

impl Primitive for Polygon {
    fn brdf(&self) -> &Brdf {
        &self.brdf
    }
    fn intersect(&self, ray_origin: &[f64; 3], ray_dir: &[f64; 3]) -> Option<[f64; 3]> {
        // game plan: intersect with plane, project onto suitable 2D coordinates to test for inclusion in polygon

        // we need a quick way to bail like we did with the sphere. Hence bounding_center and bounding_radius
        // but first, line/plane intersection

        // how much the ray moves away from the plane for every unit step forward
        let p = geometry::dot(&self.normal, &ray_dir);

        let test = geometry::dot(&geometry::sub(&ray_origin, &self.points[0]), &self.normal);
        // bail #1: ray is heading away from the plane (or parallel)
        if test <= geometry::EPSILON || p >= -geometry::EPSILON {
            return None;
        }

        let a = geometry::sub(&self.bounding_center, &ray_origin);
        // distance from plane 
        let d = &geometry::dot(&a, &self.normal).abs();

        let plane_intersection = geometry::add(
            &ray_origin,
            &geometry::scale(d/p.abs(), &ray_dir)
        );

        // bail #2: 
        if geometry::d(&plane_intersection, &self.bounding_center) > self.bounding_radius {
            return None;
        }

        // write intersection in terms of basis. We use the first point as the origin
        let intersection_offset = geometry::sub(&plane_intersection, &self.points[0]);
        let intersection_2d = [
            geometry::dot(&self.basis_1, &intersection_offset), 
            geometry::dot(&self.basis_2, &intersection_offset)
        ];

        // now we're down to 2d polygon intersection
        if geometry::poly_inclusion_2d(&self.points_in_basis, &intersection_2d) {
            return Some(plane_intersection);
        }

        return None;
    }

    fn normal_at(&self, point:&[f64;3]) -> [f64;3] {
        self.normal
    }

    // self.points must be non-empty
    fn package(&self) -> [f64;4] {
        let mut meanPt = [0.0, 0.0, 0.0];
        for pt in &self.points {
            meanPt = geometry::add(&meanPt, &pt);
        }
        meanPt = geometry::scale(1.0 / (self.points.len() as f64), &meanPt);

        let mut radius = 0.0;
        for pt in &self.points {
            radius = geometry::max(radius, geometry::d(&pt, &meanPt));
        }
        [meanPt[0], meanPt[1], meanPt[2], radius]
    }
}