

let mathjs = require('mathjs');
let colorConvert = require('color-convert');

class Camera{
    constructor(location, orientation, lensWidth){
        this.location = location;
        this.orientation = orientation;
        this.lensWidth = lensWidth;
    }

    getRays(resolution){
        // divided viewing space into resolution^2 little squares
        // by default the camera points in the -z direction
        let rays = [];
        let d = this.lensWidth / resolution;
        let countX = 0, countY = 0;

        let rotate = (pt, orientation) => {
            return matColMult(orientation, pt);
        }

        for(let i = -this.lensWidth/2 + d/2; i < this.lensWidth/2; i += d){
            countX++;
            countY = 0;
            for(let j = -this.lensWidth/2 + d/2; j < this.lensWidth/2; j += d){
                countY++;
                // this is the center of the pixel
                let rawPt = [i,j,1];
                let childrenPts = [[i+d/3,j+d/3,1], [i-d/3,j+d/3,1], [i+d/3,j-d/3,1], [i-d/3,j-d/3,1]];
                // transform according to the orientation, then translate to the location
                let rotatedPt = rotate(rawPt, this.orientation);
                let rotatedChildren = childrenPts.map(x => rotate(x, this.orientation));

                let newRay = new Ray(plus(this.location,rotatedPt), rotatedPt);
                newRay.children = rotatedChildren.map(x => new Ray(plus(this.location,x), x));
                newRay.pixel = [countX-1, countY-1];
                rays.push(newRay);
            }
        }

        return rays;
    }
}

class Ray{
    // if maxRange < 0, it is infinite
    // start, direction are vectors in R^3 (list)
    constructor(start, direction) {
        this.start = start;
        this.direction = direction;
        if(mag(direction) === 0){
            throw new Error('direction must be nonzero');
        }
    }

    deserialize(){
        if(this.children){
            for(let child of this.children){
                child.__proto__ = Ray.prototype;
                child.deserialize();
            }
        }
    }

    evaluate(t){
        return plus(this.start, scale(t, this.direction));
    }
}

class FeatureCollection{
    constructor(features){
        this.subFeatures = features;
        this.type = 'FeatureCollection'
    }

    deserialize(){
        for(let i = 0; i < this.subFeatures.length; i++){
            let t = this.subFeatures[i].type;
            if(t === 'FeatureCollection'){
                this.subFeatures[i].deserialize();
            }else if(t === 'Sphere'){
                this.subFeatures[i].__proto__ = Sphere.prototype;
            }else if(t === 'PolyPlane'){
                this.subFeatures[i].__proto__ = PolyPlane.prototype;
            }
        }
    }

    intersect(ray){
        let intersections = [];
        for(let f of this.subFeatures){
            let subIntersections = f.intersect(ray);
            for(let subI of subIntersections){
                intersections.push(subI);
            }
        }

        intersections.sort((a,b) => {
            return mag(minus(a.pos, ray.start)) - mag(minus(b.pos, ray.start));
        })

        return intersections;
    }
}

// intersect function returns list of intersections with given ray (t>0), sorted 
// by t
class Sphere{
    constructor(center, radius, properties) {
        this.center = center;
        this.radius = radius;
        this.properties = properties;
        this.type = 'Sphere'
    }

    static uniform(center, radius, opacity, specularity, color, indexOfRefraction){
        let properties = {
            opacity: opacity,
            specularity: specularity,
            color: color,
            indexOfRefraction: indexOfRefraction
        };

        return new Sphere(center,radius,properties);
    }

    // intersection object:
    // {pos, surface:{normal, opacity, specularity, color, indexOfRefraction}}
    intersect(ray){
        // L = ray.start + t ray.direction
        // |x - center| = radius**2
        // transform so that the sphere is at the origin and the ray points in the positive x direction
        let relPos = minus(ray.start, this.center);
        // rotate so that direction points in the positive x direction
        // this is a linear transformation mapping ray.direction to (|ray.direction|,0,0)
        // see https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        
        // a simple way to check if it's going to intersect at all is find the minimum distance 
        // from the point to the ray and see if it's less than radius
        let vec1 = minus(this.center, ray.start);
        let dist = mag(minus(vec1, project(vec1, ray.direction)));
        if(dist > this.radius){
            return [];
        }

        let s = [0,0,0];
        if(ray.direction[1] === 0 && ray.direction[2] === 0){
            if(ray.direction[0] > 0){
                s = relPos;
            }else{
                s = scale(-1, relPos);
            }
        }else{
            let a = normalize(ray.direction);
            let b = [1,0,0];
            let v = cross(a,b);
            let c = dot(a,b);
            let vMat = [
                [0,-v[2],v[1]],
                [v[2],0,-v[0]],
                [-v[1],v[0],0]
            ]
            let vMatSquared = [
                [-1*v[1]**2-v[2]**2, v[0]*v[1], v[0]*v[2]],
                [v[0]*v[1], -1*v[0]**2-v[2]**2, v[1]*v[2]],
                [v[0]*v[2], v[1]*v[2], -1*v[0]**2-v[1]**2]
            ];

            let identity = [[1,0,0],[0,1,0],[0,0,1]];
            let transformationMatrix = matPlus(identity, matPlus(vMat, matScale(1/(1+c), vMatSquared)));
            s = matColMult(transformationMatrix.valueOf(), relPos);
        }

        let discriminant = this.radius**2 - s[1]**2 - s[2]**2;
        if(discriminant < 0){
            return [];
        }

        let magDir = mag(ray.direction);
        let t1 = (-s[0] - Math.sqrt(discriminant))/magDir;
        let t2 = (-s[0] + Math.sqrt(discriminant))/magDir;

        function intersectObject(loc, thisRef){
            let toReturn = {
                pos: loc,
                object:thisRef,
                surface:{}
            }

            if(!thisRef.properties.lightSource){
                toReturn.surface = {
                    normal: normalize(minus(loc, thisRef.center)),
                    opacity: thisRef.properties.opacity,
                    specularity: thisRef.properties.specularity,
                    color: thisRef.properties.color,
                    indexOfRefraction: thisRef.properties.indexOfRefraction
                };
            }

            return toReturn;
        }

        let res = [];
        if(t1 > 0){
            // t1 <= t2
            res.push(intersectObject(ray.evaluate(t1), this));
        }
        // don't want to add a duplicate
        if(t2 > 0 && discriminant !== 0){
            res.push(intersectObject(ray.evaluate(t2), this));
        }

        return res;
    };
}

// please make sure all the points are in a plane and no three 
// adjacent pts are collinear
class PolyPlane{
    constructor(pts, orientationVector, properties){
        this.points = pts;
        this.normal = normalize(cross(minus(pts[2],pts[1]), minus(pts[1],pts[0])));
        this.properties = properties;
        this.type = 'PolyPlane'
        if(dot(this.normal, orientationVector) < 0){
            this.normal = scale(-1, this.normal);
        }
    }

    static uniform(pts, orientationVector, opacity, specularity, color, indexOfRefraction){
        let properties = {
            opacity: opacity,
            specularity: specularity,
            color: color,
            indexOfRefraction: indexOfRefraction
        };

        return new PolyPlane(pts,orientationVector,properties);
    }

    intersect(ray){
        // plane-ray intersection followed by planar polygon intersection
        let den = dot(ray.direction, this.normal);
        if(den === 0){
            // this plane is infinitely thin
            return [];
        }

        let t = dot(this.normal, minus(this.points[0], ray.start)) / den;

        if(t <= 0){
            return [];
        }

        let pt = ray.evaluate(t);

        // project the pts onto a suitable plane 
        let dot1 = Math.abs(dot(this.normal, [1,0,0]));
        let dot2 = Math.abs(dot(this.normal, [0,1,0]));
        let dot3 = Math.abs(dot(this.normal, [0,0,1]));
        let removeIndex = 0;
        if(dot2 > dot1 && dot2 > dot3){
            removeIndex = 1;
        }
        if(dot3 > dot1 && dot3 > dot2){
            removeIndex = 2;
        }

        let polyPoints = [];
        let testPoint = [];
        for(let i = 0; i < this.points.length; i++){
            let newPt = [];
            for(let j = 0; j < 3; j++){
                if(j !== removeIndex){
                    newPt.push(this.points[i][j]);
                }
            }
            polyPoints.push(newPt);
        }
        for(let i = 0; i < 3; i++){
            if(i !== removeIndex){
                testPoint.push(pt[i]);
            }
        }

        let segments = [];
        for(let i = 0; i < polyPoints.length; i++){
            segments.push([polyPoints[i], i === 0 ? polyPoints[polyPoints.length - 1] : polyPoints[i-1]]);
        }

        let count = 0;
        let p = testPoint;
        for(let s of segments){
            let [u,w] = s;
            // x component of intersection, that is
            if((p[1] > Math.min(u[1],w[1]) && p[1] < Math.max(u[1],w[1])) && u[0]+((p[1]-u[1])/(w[1]-u[1]))*(w[0]-u[0]) > p[0]){
                count++;
            }
        }

        function intersectObject(loc, thisRef){
            let toReturn = {
                pos: loc,
                object:thisRef,
                surface:{}
            };

            if(!thisRef.properties.lightSource){
                toReturn.surface = {
                    normal: thisRef.normal,
                    opacity: thisRef.properties.opacity,
                    specularity: thisRef.properties.specularity,
                    color: thisRef.properties.color,
                    indexOfRefraction: thisRef.properties.indexOfRefraction
                }
            };

            return toReturn;
        }

        if(count % 2 === 1){
            // within the polygon
            return [intersectObject(pt, this)];
        }
        return [];

    }
}

function traceRay(featureCollection, environment, lightSources, ray, maxRecursion, currentDepth=0){
    if(currentDepth === maxRecursion){
        // max recursion achieved
        ray.color = ray.medium.color;
        return ray;
    }

    let collisions = featureCollection.intersect(ray);
    if(collisions.length === 0){
        // it hit nothing. Return the inky blackness (or whatever color it is) of the abyss
        ray.color = [0,0,0];
        for(let i = 0; i < 3; i++){
            ray.color[i] = ray.medium.color[i] + (Math.random()-.5)*10;
            ray.color[i] = Math.max(0, ray.color[i]);
            ray.color[i] = Math.min(255, ray.color[i]);
        }
        return ray;
    }
    // the first one
    let collision = collisions[0];

    if(collision.object.properties.lightSource){
        ray.brightness = collision.object.properties.lightSource.surfaceIntensity;
        ray.color = collision.object.properties.lightSource.color;
        return ray;
    }

    // first add ambient lighting
    ray.brightness = environment.ambientConstant;
    ray.color = [0,0,0];

    // next calculate diffuse/specular lighting to all the light sources
    for(let lightFeature of lightSources){
        let light = lightFeature.properties.lightSource;
        let diffuseIntensity = light.totalIntensity / light.diffuseSample.length;
        for(let lightPoint of light.diffuseSample){
            // determine if there is a path to the light source
            let directionToLight = normalize(minus(lightPoint, collision.pos));
            let newRay = new Ray(
                plus(collision.pos,scale(1e-8,directionToLight)), 
                directionToLight
            );
            newRay.medium = ray.medium;

            let testCollisions = featureCollection.intersect(newRay);
            if(testCollisions.length > 0 && testCollisions[0].object === lightFeature){
                // use this light source sample point
                let distToLight = mag(minus(lightPoint, collision.pos));
                let mediumOpacityConst = Math.min(1, ray.medium.opacity * distToLight);
                // it is important that directionToLight and collision.surface.normal are both unit vectors
                // diffuse
                let dBrightness = environment.diffuseConstant*dot(directionToLight, collision.surface.normal) * diffuseIntensity / distToLight**2;
                ray.brightness += dBrightness;
                // specular
                let reflectDir = normalize(plus(
                    scale(2, 
                        minus(
                            project(directionToLight, collision.surface.normal),
                            directionToLight
                        )
                    ), 
                    directionToLight
                ));
                let c1 = Math.abs(dot(reflectDir, normalize(ray.direction)))**environment.specularNarrownessConstant
                ray.brightness += environment.specularStrengthConstant * c1 * diffuseIntensity*collision.surface.specularity / distToLight**2;
            }
        }
    }

    // now do a new ray!
    let outwardDir = scale(-1, ray.direction);
    let bounceDirection = plus(
        scale(2, 
            minus(
                project(outwardDir, collision.surface.normal),
                outwardDir
            )
        ), 
        outwardDir
    );

    // Specular Reflection Ray
    let childRay = new Ray(plus(collision.pos, scale(1e-8,bounceDirection)), bounceDirection);
    // we say the ray is leaving the object if the normal is in the direction of the ray
    let leaving = dot(childRay.direction, collision.surface.normal) > 0;
    if(leaving){
        childRay.medium = environment.medium;
    }else{
        childRay.medium = ray.medium;
    }
    let resultChild = traceRay(featureCollection,environment,lightSources, childRay, maxRecursion, currentDepth + 1);

    // this ray is on the other side of the surface
    // snell's law! 
    let n1 = ray.medium.indexOfRefraction, n2 = 0;
    let refractionIntoSurface = leaving;
    if(!refractionIntoSurface){
        n2 = environment.medium.indexOfRefraction;
    }else{
        n2 = collision.surface.indexOfRefraction;
    }
    let normalWithRay = project(ray.direction, collision.surface.normal);
    let costheta1 = dot(ray.direction, normalWithRay) / (mag(ray.direction)*mag(normalWithRay));
    let sintheta1 = Math.sqrt(1-costheta1**2);
    let sintheta2 = (n1/n2)*sintheta1;
    let totalInternalReflection = false;
    let childRay2 = null;
    if(Math.abs(sintheta2) > 1){
        totalInternalReflection = true;
    }else{
        let theta2 = Math.asin(sintheta2);
        let a = mag(normalWithRay);
        let theta1 = Math.acos(costheta1);
        let f = a*(Math.tan(theta1) - Math.tan(theta2));
        let refractionDirection = plus(
            ray.direction,
            scale(f, normalize(minus(normalWithRay, ray.direction)))
        );
        childRay2 = new Ray(plus(collision.pos, scale(1e-8,normalWithRay)), refractionDirection);
        childRay2.medium = !refractionIntoSurface ? environment.medium : {indexOfRefraction: n2, opacity: collision.surface.opacity, color: collision.surface.color};
        if(refractionIntoSurface && collision.surface.opacity === 1){
            childRay2.color = collision.surface.color;
        }else{
            childRay2 = traceRay(featureCollection,environment,lightSources, childRay2, maxRecursion, currentDepth + 1);
        }
    }

    // combine the colors
    let rayDist = mag(minus(ray.start, collision.pos));
    let opacityConst = Math.min(ray.medium.opacity * rayDist, 1);
    for(let i = 0; i < 3; i++){
        if(ray.color[i] > 255){ray.color[i] = 255;}
        let a1 = collision.surface.color[i]*ray.brightness;
        let a2 = resultChild.color[i];
        let a3 = childRay2.color[i];
        let opacity = !refractionIntoSurface ? environment.medium.opacity : collision.surface.opacity;
        if(totalInternalReflection){
            opacity = 1;
        }
        // this function (z = x + y - xy) satisfies some nice properties
        let normalConst = 1 - (collision.surface.specularity + (1-opacity) - collision.surface.specularity*(1-opacity));
        let specConst = 0, opacConst = 0;
        if((collision.surface.specularity + (1-opacity)) !== 0){
            specConst = (collision.surface.specularity / (collision.surface.specularity + (1-opacity))) * (1 - normalConst);
            opacConst = ((1-opacity) / (collision.surface.specularity + (1-opacity))) * (1 - normalConst);
        }
        ray.color[i] =
            (1-opacityConst)*Math.min(normalConst * a1 + specConst * a2 + opacConst * a3, 255) +
            opacityConst * ray.medium.color[i];
    }

    return ray;

}

function dot(a,b){
    let t = 0;
    for(let i = 0; i < a.length; i++){
        t += a[i]*b[i];
    }
    return t;
}

function cross(a,b){
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}

function plus(...args){
    let res = [];
    for(let i = 0; i < args[0].length; i++){
        let total = 0;
        for(let j = 0; j < args.length; j++){
            total += args[j][i];
        }
        res.push(total);
    }
    return res;
}

// so minus was taking a lot of time, so we're going to ditch map 
// for a less general version (sorry!)
function minus(a,b){
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
}

// see minus for comments on speed
function scale(a, v){
    return [a*v[0], a*v[1], a*v[2]];
}

function mag(a){
    return Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

function normalize(a){
    let m = mag(a);
    if(m === 0){
        throw new Error('attempt to normalize zero vector');
    }
    return scale(1/m, a);
}

// turns out mathjs is pretty slow with matrix multiplication
// mat is a 3x3 matrix, colVec is a 1D column vector (in R^3)
function matColMult(mat, colVec){
    return [dot(mat[0], colVec), dot(mat[1], colVec), dot(mat[2], colVec)];
}

function matScale(a, mat){
    newMat = [];
    for(let i = 0; i < 3; i++){
        newMat.push(scale(a, mat[i]));
    }
    return newMat;
}

function matPlus(mat1, mat2){
    newMat = [];
    for(let i = 0; i < 3; i++){
        newMat.push(plus(mat1[i], mat2[i]));
    }
    return newMat;
}

// 2 3x3 matrices
function matMatMult(mat1, mat2){
    let col1 = [], col2 = [], col3 = [];
    for(let i = 0; i < 3; i++){
        col1.push(mat2[0][i]);
        col2.push(mat2[1][i]);
        col3.push(mat2[2][i]);
    }

    return [
        [dot(mat1[0], col1), dot(mat1[0], col2), dot(mat1[0], col3)],
        [dot(mat1[1], col1), dot(mat1[1], col2), dot(mat1[1], col3)],
        [dot(mat1[2], col1), dot(mat1[1], col2), dot(mat1[2], col3)]
    ];
}

// a onto b, that is
function project(a,b){
    return scale(dot(a,b)/dot(b,b), b);
}

module.exports = {
    Sphere: Sphere,
    PolyPlane: PolyPlane,
    FeatureCollection: FeatureCollection,
    Camera: Camera,
    Ray: Ray,
    dot:dot,
    minus:minus,
    plus:plus,
    project:project,
    normalize:normalize,
    mag:mag,
    cross:cross,
    scale:scale,
    traceRay: traceRay
}
