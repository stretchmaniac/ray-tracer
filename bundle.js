(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
let gaussianRand = require('gauss-random');

class Camera{
    // we model the camera as a sensor behind a convex lens. We need to know a couple things
    // 1. focal length of the lens
    // 2. diagonal of sensor
    // 3. distance sensor --> lens (for zoom)
    // 4. the radius of the lens (aperture)
    // a so-called "normal" camera (i.e. modelled kinda after the eye) has a sensor diagonal
    //  roughly equal to the focal length
    // sensorLocation is the middle of the 2D plane (square) sensor, orientation is a 3 x 3 (orthogonal) rotation
    // matrix that represents rotation from the positive z axis (sensorLocation --> lensLocation)
    constructor(sensorLocation, orientation, focalLength, sensorLensDist, sensorDiagonal, apertureRadius, resolution){
        this.location = sensorLocation;
        this.orientation = orientation;
        this.focalLength = focalLength;
        this.sensorWidth = sensorDiagonal / Math.sqrt(2);
        this.aperture = apertureRadius;
        this.z = sensorLensDist;
        this.resolution = resolution;
    }

    // suppose the detector lies on the interval [-1,1] x [-1,1]
    // to get depth of field, we must consider that a ray converging on a point 
    // on the sensor may come from any point of the lens. Therefore we return a generator 
    // that randomly samples rays that converge on unitLocation. Note that the generated 
    // rays emanate from the plane of the lens, not the sensor
    getRay(unitLocation){
        // there is a question of whether light converging on unitLocation originates uniformly
        // from the lens. It doesn't, and there's some info here http://www.timledlie.org/cs/graphics/finalproj/finalproj.html
        // (specifically one of the papers it links to) about why that is. However, depth of field 
        // is still a reach goal at the moment so we will not worry so much about it
        
        // transform unitLocation to world coordinates (location on sensor)
        // -x and -z to correct for image mirroring (the brain flips the image)
        let sensorPt = [-unitLocation[0]*this.sensorWidth/2,-unitLocation[1]*this.sensorWidth/2,0];
        sensorPt = plus(this.location, matColMult(this.orientation, sensorPt));

        const centerOfLens = plus(this.location, matColMult(this.orientation, [0,0,this.z]));

        //console.log(sensorPt, centerOfLens);
        //console.log(a.b.c);
        
        // lens equation: 1/z - 1/z' = 1/f ==> z' = 1/(1/z - 1/f)
        // the question is what "z" is here. If we simply use the focal length 
        // then we should get a flat plane of "in focus" stuff. If we use the distance 
        // from the sensor location to the middle of the lens we would get a sphere of 
        // focus. In practice I do not think it makes much of a difference
        const zPrime = -1/(1/this.z - 1/this.focalLength);
        const imageLoc = plus(sensorPt, scale(zPrime + this.z, normalize(minus(centerOfLens, sensorPt))));

        // sample from the lens
        let pt = [1,1];
        while(pt[0]**2 + pt[1]**2 > 1){
            pt = [Math.random(), Math.random()];
        }
        pt = [pt[0]*this.aperture, pt[1]*this.aperture];
        let lensPt = [pt[0], pt[1], this.z];
        // rotate, translate 
        lensPt = matColMult(this.orientation, lensPt);
        lensPt = plus(lensPt, this.location);

        return new Ray(lensPt, normalize(minus(imageLoc, lensPt)));
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
        // default medium 
        this.medium = {
            color:[0,0,0]
        };
    }

    evaluate(t){
        return plus(this.start, scale(t, this.direction));
    }
}

class FeatureCollection{
    constructor(features){
        this.subFeatures = features;
        this.type = 'FeatureCollection';

        // for sampling
        this.totalSurfaceArea = 0;
        for(let f of features){
            this.totalSurfaceArea += f.area;
        }
    }

    deserialize(){
        for(let i = 0; i < this.subFeatures.length; i++){
            let t = this.subFeatures[i].type;
            if(t === 'FeatureCollection'){
                this.subFeatures[i].deserialize();
            }else if(t === 'Sphere'){
                this.subFeatures[i].__proto__ = Sphere.prototype;
                this.subFeatures[i].deserialize();
            }else if(t === 'PolyPlane'){
                this.subFeatures[i].__proto__ = PolyPlane.prototype;
                this.subFeatures[i].deserialize();
            }
        }
    }

    sample(){
        // pick a feature and sample from it 
        let r = Math.random();
        let cumulativeArea = 0;
        let feature = undefined;
        for(let f of this.subFeatures){
            cumulativeArea += f.area / this.totalSurfaceArea;
            if(cumulativeArea > r){
                feature = f;
                break;
            }
        }
        if(feature === undefined){
            console.log('ERROR: random feature collection algorithm not working.');
        }
        const s = feature.sample();
        s.object = feature;
        return s;
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

    // test whether the line segment [p1, p2] is uninterrupted
    lineSegmentValid(p1, p1Obj, p2, p2Obj){
        // special case for PolyPlane
        if(p1Obj === p2Obj && p1Obj.type === 'PolyPlane'){
            return false;
        }
        let ray = new Ray(p1, minus(p2, p1));
        const eps = 1e-8;
        ray.start = plus(ray.start, scale(eps, ray.direction));
        const end = minus(p2, scale(eps, ray.direction));
        let intersections = this.intersect(ray);
        if(intersections.length === 0){
            return true;
        }
        return mag(minus(ray.start,end)) < mag(minus(intersections[0].pos, p1));
    }
}

// bidirectional reflectance distribution function
class BRDF{
    // n ranges from 0 to 1, where 0 is perfectly reflective and 
    // 1 is perfectly diffuse
    static variableGlossy(n){
        let brdf = {};
        brdf.method = 'quick_n_dirty';
        if(n === 1){
            brdf.method = 'ideal_diffuse';
        }
        brdf.properties = {
            gloss:1-n,
            lobeSize:(2*n)**2
        };
        brdf.__proto__ = BRDF.prototype;
        return brdf;
    }

    // acts like a coated glass surface
    static variableTransparent(indexOfRefraction, frosty){
        let brdf = {};
        // transparency is pretty simple. We have transmitted and reflected 
        // rays that follow the Fresnel equations 
        // https://en.wikipedia.org/wiki/Fresnel_equations
        // methinks we won't be modelling the polarization of the incoming wave so 
        // we will assume that it is random
        brdf.method = 'transparent';
        brdf.properties = {
            n:indexOfRefraction,
            lobeSize:(2*frosty)**2
        };
        brdf.__proto__ = BRDF.prototype;
        return brdf;
    }

    // reflection/transmission probabilities and 
    // directions given a perfectly transparent (non-frosty) surface
    rt(inVec, normal){
        // we assume a ray with power 1, random polarization in the s and p direction 
        const sPol = Math.random();
        const pPol = 1 - sPol;
        const theta = angleBetweenVecs(normal, inVec);
        // check for total internal reflection
        // n1 -- inside, n2 -- outside
        // n1 sin t = n2 sin t2 ==> t2 = arcsin((n1/n2)sin(t))
        // we assume n2 = 1 for now (and generally assume air layers of negligible thickness around everything)
        let n1 = undefined;
        let n2 = undefined;
        if(dot(inVec, normal) < 0){
            // ray is outside surface
            n1 = 1;
            n2 = this.n;
        }else{
            // ray is inside surface
            n1 = this.n;
            n2 = 1;
        }
        const inner = Math.abs((n1/n2)*Math.sin(theta));
        if(inner > 1){
            //total internal reflection 
            return {
                // reflect, direction, probability
                r:{dir: normalize(reflect(inVec, normal)), p: 1},
                // transmit, direction, probability
                t:{dir: [0,0,1], p:0}
            }
        }

        // calculate refraction angle 
        const refractionAngle = Math.asin(inner);

        //console.log(inVec, normal, theta);

        // calculate reflectance
        const Rs = ((n1*Math.cos(theta)-n2*Math.cos(refractionAngle)) / (n1*Math.cos(theta)+n2*Math.cos(refractionAngle)))**2;
        const Rp = ((n1*Math.cos(refractionAngle)-n2*Math.cos(theta)) / (n1*Math.cos(refractionAngle)+n2*Math.cos(theta)))**2;

        const totalReflectance = sPol*Rs + pPol*Rp;
        return {
            r: {dir:normalize(reflect(inVec, normal)), p:totalReflectance},
            t: {dir:normalize(snell(inVec, normal, refractionAngle)), p:1-totalReflectance}
        };
    }

    // gives a sample of the most likely vectors to 
    // reflect in the direction of outVec
    // should return a list of rays
    inverseReflectionSample(normal, outVec, numToSample){
        if(this.method === 'quick_n_dirty' || this.method === 'transparent' || this.method === 'ideal_diffuse'){
            // we assume that the distribution is symmetric, so an inverse 
            // sample is the same as a normal sample 
            let samples = [];
            let newIn = scale(-1, outVec);
            for(let k = 0; k < numToSample; k++){
                samples.push(new Ray([0,0,0], this.sample(newIn, normal)));
            }
            return samples; 
        }
    }

    sample(inVec, normal){
        if(this.method === 'quick_n_dirty'){
            // a little hacky, but you should have seen what I had before
            const ref = normalize(reflect(inVec, normal)); 
            return BRDF.gaussDir(ref, this.properties.lobeSize, normal);
        }else if(this.method === 'transparent'){
            const probs = this.rt(inVec, normal);
            // console.log(inVec, normal, probs);
            const matchesNormal = dot(normal, inVec) < 0;
            const invNormal = scale(-1, normal);
            if(Math.random() < probs.r.p){
                // reflection 
                return BRDF.gaussDir(probs.r.dir, this.properties.lobeSize, matchesNormal ? normal : invNormal);
            }else{
                // transmission 
                return BRDF.gaussDir(probs.t.dir, this.properties.lobeSize, matchesNormal ? invNormal : normal);
            }
        }else if(this.method === 'ideal_diffuse'){
            // reflection does not depend on inVec
            // instead sample the unit hemisphere
            // theta \in [0, 2\pi], uniform sampling 
            // phi \in [0, pi/2], pdf ~ sin(phi) (with z=cos(phi))
            let theta = Math.random() * Math.PI * 2;
            let phi = Math.acos(Math.random());
            // create local coordinates
            let zhat = normal;
            let xhat = normalize(perpVec(zhat));
            let yhat = cross(zhat, xhat);
            return plus(scale(Math.cos(theta)*Math.sin(phi), xhat), scale(Math.sin(theta)*Math.sin(phi), yhat), scale(Math.cos(phi), zhat));
        }
    }

    static gaussDir(normDir, lobeSize, normal){
        let validDirection = false;
        let finalPt = null;
        while(validDirection === false){
            let randomDirection = [Math.random()-.5, Math.random()-.5, Math.random()-.5];
            while(mag(randomDirection) > .5){
                randomDirection = [Math.random()-.5, Math.random()-.5, Math.random()-.5]
            }
            // this gives us our direction
            let perpComponent = minus(randomDirection, project(randomDirection, normDir));
            perpComponent = normalize(perpComponent);
            // now our angle difference 
            let newAngleOffset = gaussianRand() * lobeSize;
            finalPt = plus(scale(Math.cos(newAngleOffset), normDir), scale(Math.sin(newAngleOffset), perpComponent));

            validDirection = dot(finalPt, normal) >= 0;
        }
        return finalPt;
    }

    static gaussPDF(mean, stdev, test){
        const x = test - mean;
        return (1 / Math.sqrt(2*Math.PI*stdev**2)) * Math.E**(-1*x**2 / (2*stdev**2));
    }

    // gives the associated probability that an input ray reflects in the 
    // direction of outVec
    reflectFunc(inVec, normal, outVec){
        if(this.method === 'quick_n_dirty'){
            const reflected = reflect(inVec, normal);
            const angleDiff = angleBetweenVecs(reflected, outVec);
            return BRDF.gaussPDF(0, this.properties.lobeSize, angleDiff);
        }else if(this.method === 'transparent'){
            const epsilon = 1e-6;
            const prob = this.rt(inVec, normal);
            if(Math.sign(dot(normal, inVec)) === Math.sign(dot(normal, outVec))){
                // we're looking at a transmission here
                const angleDiff = angleBetweenVecs(prob.t.dir, outVec);
                return prob.t.p * BRDF.gaussPDF(0, this.properties.lobeSize, angleDiff);
            }else{
                // looking at reflection
                const angleDiff = angleBetweenVecs(prob.r.dir, outVec);
                return prob.r.p * BRDF.gaussPDF(0, this.properties.lobeSize, angleDiff);
            }
        }else if(this.method === 'ideal_diffuse'){
            return 1/(2*Math.PI);
        }
    }
}

// "Cone" ( --Audrey)
class Frustum{
    constructor(base1Center, base2Center, r1, r2, properties){
        let normedHeight = normalize(minus(base2Center, base1Center));
        let normedV1 = normalize(perpVec(normedHeight));
        let normedV2 = cross(normedHeight, normedHeight, normedV1);
        this.base1 = new Disk(base1Center, scale(r1, normedV1), scale(r1, normedV2), normedHeight, properties);
        this.base2 = new Disk(base2Center, scale(r2, normedV1), scale(r2, normedV2), scale(-1, normedHeight), properties);
        this.base1Center = base1Center;
        this.base2Center = base2Center; 
        this.r1 = r1;
        this.r2 = r2;
        this.properties = properties;
        this.type = 'Frustum';
        let h = mag(minus(base2Center, base1Center));
        this.height = h;
        this.area = Math.PI*(r1**2 + r2**2 + (r1 + r2)*Math.sqrt((r2-r1)**2 + h**2)); 
    }

    static uniform(base1Center, base2Center, r1, r2, color, lightSource){
        const properties = {
            color: color,
            lightSource: lightSource
        };
        return new Frustum(base1Center, base2Center, r1, r2, properties);
    }

    deserialize(){
        if(this.properties.brdf){
            this.properties.brdf.__proto__ = BRDF.prototype;
        }
    }

    sample(){
        let b1Area = Math.PI*this.r1**2;
        let b2Area = Math.PI*this.r2**2;
        let otherArea = this.area - b1Area - b2Area;
        let b1CutOff = b1Area / this.area;
        let b2CutOff = (b1Area + b2Area) / this.area;
        let r = Math.random();
        if(r < b1CutOff){
            // pick from first disk 
            return this.base1.sample();
        }
        if(r < b2CutOff){
            // pick from second disk 
            return this.base2.sample();
        }
        // choose from slanted part. This is 
        // essentially sampling from a sector of 
        // a circle
        let a = this.r2*this.height / (this.r1 - this.r2);
        let bigR = Math.sqrt(this.r1**2 + (this.height + a)**2);
        let littleR = Math.sqrt(this.r2**2 + a**2);
        let sectorAngle = 2*Math.PI*this.r1/bigR;
        // r \in [littleR, bigR]
        // \theta \in [0, sectorAngle]
        while(true){
            let sampledR = Math.random() * (bigR - littleR) + littleR;
            let sampledT = Math.random() * sectorAngle;
            // keep with probability sampledR / bigR
            // see https://mathoverflow.net/questions/9991/how-can-i-sample-uniformly-from-a-surface
            if(Math.random() <= sampledR / bigR){
                let j = (sampledR-littleR)/(bigR-littleR);
                let flatR = j*(this.r2-this.r1)+this.r1; 
                let normedV1 = normalize(this.base1.v1);
                let normedV2 = nromalize(this.base1.v2);
                let actualPt = plus(
                    plus(this.base1Center, scale(j, minus(this.base2Center, this.base1Center))),
                    plus(scale(flatR*Math.cos(sampledT), normedV1), scale(flatR*Math.sin(sampledT), normedV2))
                );
                // find 2 tangent vectors and compute cross product to get normal vector
                // first tangent: vector to bottom
                let bottomPt = plus(this.base1Center, 
                    plus(
                        scale(Math.cos(sampledT), this.base1.v1), 
                        scale(Math.sin(sampledT), this.base1.v2)
                    )
                );
                let tangent1 = minus(actualPt, bottomPt);
                // second tangent: around radius of circle
                let component1 = project(actualPt, this.base1.v1);
                let component2 = project(actualPt, this.base1.v2);
                // negatives are important. The cross product was designed to create an 
                // outward pointing normal
                let newC1 = component2;
                let newC2 = -component1;
                let tangent2 = plus(scale(newC1, this.base1.v1), scale(newC2, this.base1.v2));
                let normal = normalize(cross(tangent1, tangent2));

                return {pt: actualPt, dir: normal, color:this.properties.color}
            }
        }
    }

    intersect(ray){
        // transform to origin 
        
    }
}

class Disk{
    // |v1| = |v2| = r, v1 and v2 span the plane that the disk lies in
    // normalDir has to point to one side of the disk
    constructor(center, v1, v2, normalDir, properties){
        this.center = center;
        this.radius = mag(v1);
        // create orthogonal basis 
        this.v1 = v1;
        this.v2 = cross(v1, cross(v1,v2));
        this.v2 = scale(this.radius, normalize(this.v2));
        this.normal = cross(this.v1, this.v2);
        if(dot(this.normal, normalDir) < 0){
            this.normal = scale(-1, this.normal);
        }
        this.normal = normalize(this.normal);
        this.properties = properties;
        this.type = 'Disk';
        this.area = Math.PI*this.radius**2;
    }

    static uniform(center, radiusVec, color, lightSource){
        const properties = {
            color:color,
            lightSource: lightSource
        };
        return new Disk(center, radiusVec, properties);
    }

    deserialize(){
        if(this.properties.brdf){
            this.properties.brdf.__proto__ = BRDF.prototype;
        }
    }

    sample(){
        // sample from unit circle and extend to disk
        let pt = [2*Math.random() - 1, 2*Math.random() - 1];
        while(pt[0]**2 + pt[1]**2 > 1){
            pt = [2*Math.random() - 1, 2*Math.random() - 1];
        }

        let actualPt = plus(scale(pt[0], this.v1), scale(pt[1], this.v2));
        return {pt: actualPt, dir: BRDF.gaussDir(this.normal, this.properties.brdf.properties.lobeSize, this.normal), color:this.properties.color};
    }

    intersect(ray){
        // transform so that disk lies on xy plane
        // first translate by center 
        let newStart = minus(ray.start, this.center);
        // so this linear transformation maps this.normal to [0,0,1]
        let newDir = matchTransformation(ray.direction, this.normal, [0,0,1]);
        // this was a distance preserving transformation 
        if(newDir[2] === 0){
            // this is more of a convention than anything
            return [];
        }
        let t = -newStart[2] / newDir[2];
        let intersection = plus(newStart, scale(t, newDir));
        if(mag(intersection) <= this.radius){
            // transform back
            let actualIntersection = matchTransformation(intersection, [0,0,1], this.normal);

            let toReturn = {
                pos: actualIntersection,
                object:this,
                surface:{}
            };
            toReturn.surface = {
                normal: this.normal,
                color: this.properties.color
            }
            return [toReturn];
        }
        return [];
    }
}

class Shell{
    constructor(center, outerR, innerR, properties){
        this.center = center;
        this.outerRadius = outerR;
        this.innerRadius = innerR;
        this.properties = properties;
        this.type = 'Shell';
        this.area = 4*Math.PI*(this.outerRadius**2 + this.innerRadius**2);

        // for use in intersection algorithm
        this.innerSphere = new Sphere(center, innerR, properties);
        this.outerSphere = new Sphere(center, outerR, properties);
    }

    static uniform(center, outerR, innerR, color, lightSource){
        let properties = {
            color: color,
            lightSource: lightSource
        };
        return new Shell(center, outerR, innerR, properties);
    }

    deserialize(){
        if(this.properties.brdf){
            this.properties.brdf.__proto__ = BRDF.prototype;
        }
    }

    sample(){
        // strategy very similar to that of a sphere
        let pt = [2*Math.random()-1, 2*Math.random()-1, 2*Math.random()-1];
        while(mag(pt) > 1){
            pt = [2*Math.random()-1, 2*Math.random()-1, 2*Math.random()-1];
        }
        const direction = normalize(pt);

        // choose a radius 
        const probFirstRadius = 4*Math.PI*this.innerRadius**2 / this.area;
        const isInner = Math.random() < probFirstRadius;
        const chosenRadius = isInner ? this.innerRadius : this.outerRadius;

        // map onto surface of sphere 
        pt = plus(this.center, scale(chosenRadius, direction));
        let dir = isInner ? scale(-1, direction) : direction;
        dir = BRDF.gaussDir(dir, this.properties.brdf.properties.lobeSize, dir);
        return {pt: pt, dir: dir, color:this.properties.color};
    }

    intersect(ray){
        // we separate this into 3 cases
        // 1. If the ray starts outside the outer sphere, then compute the intersection 
        //    as if there was only the outer sphere (it would be the first intersection)
        // 2. If the ray starts on the inside of the inner sphere, then compute the intersection
        //    as if there was only the inner sphere (there is guaranteed to be an intersection)
        // 3. If the ray is in between the inner an outer sphere, there is guaranteed to 
        //    be an intersection; use the direction of the ray to decide

        function flipNormals(intersections){
            for(let i of intersections){
                i.surface.normal = scale(-1, i.surface.normal);;
            }
        }

        const rayM = minus(ray.start, this.center);
        const rayRad = mag(rayM);
        if(rayRad >= this.outerRadius){
            return this.outerSphere.intersect(ray);
        }
        if(rayRad <= this.innerRadius){
            let intersections = this.innerSphere.intersect(ray);
            // flip the normals on this intersection 
            flipNormals(intersections);
            return intersections;
        }
        const innerIntersections = this.innerSphere.intersect(ray);
        if(innerIntersections.length === 0){
            return this.outerSphere.intersect(ray);
        }
        // flip the normal 
        flipNormals(innerIntersections);
        return innerIntersections;
    }
}

// intersect function returns list of intersections with given ray (t>0), sorted 
// by t
class Sphere{
    constructor(center, radius, properties) {
        this.center = center;
        this.radius = radius;
        this.properties = properties;
        this.type = 'Sphere';
        this.area = 4*Math.PI*this.radius**2;
    }

    static uniform(center, radius, color, lightSource){
        let properties = {
            color: color,
            lightSource:lightSource,
        };

        return new Sphere(center,radius,properties);
    }

    deserialize(){
        if(this.properties.brdf){
            this.properties.brdf.__proto__ = BRDF.prototype;
        }
    }

    // returns random point on surface of sphere along with outward direction
    sample(){
        // the idea is to sample a unit cube until a point lies inside 
        // the unit sphere. Then project onto the unit sphere and 
        // map onto this sphere
        let pt = [2*Math.random()-1, 2*Math.random()-1, 2*Math.random()-1];
        while(mag(pt) > 1){
            pt = [2*Math.random()-1, 2*Math.random()-1, 2*Math.random()-1];
        }
        const direction = normalize(pt);
        // map onto surface of sphere 
        const projection = plus(this.center, scale(this.radius, direction));
        return {pt: projection, dir: BRDF.gaussDir(direction, this.properties.brdf.properties.lobeSize, direction), color:this.properties.color};
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
            toReturn.surface = {
                normal: normalize(minus(loc, thisRef.center)),
                color: thisRef.properties.color
            };

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
// also please use a convex polygon because I don't want to do difficult triangularization
class PolyPlane{
    constructor(pts, orientationVector, properties){
        this.points = pts;
        this.normal = normalize(cross(minus(pts[2],pts[1]), minus(pts[1],pts[0])));
        this.properties = properties;
        this.type = 'PolyPlane';
        if(dot(this.normal, orientationVector) < 0){
            this.normal = scale(-1, this.normal);
        }

        // for later use (i.e. the sample function) we want to generate a 
        // basis for this polygon and express all the points in that basis, 
        // triangularize, etc.

        // use points[1] as the origin, {points[0] - points[1], points[2] - points[1]} as basis vectors 
        const b1 = normalize(minus(this.points[0], this.points[1]));
        const b2 = normalize(minus(this.points[2], this.points[1]));

        // need to find length of projection onto each basis vector.
        // Lucky for us this is what dot product is good at!
        let points2D = this.points.map(x =>
            [dot(x,b1)/mag(b1), dot(x,b2)/mag(b2)]
        );

        this.basis = [b1, b2];
        this.points2D = points2D;

        function triangleArea(p1, p2, p3){
            // heron's formula for quick use 
            const a = Math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2);
            const b = Math.sqrt((p2[0]-p3[0])**2 + (p2[1]-p3[1])**2);
            const c = Math.sqrt((p1[0]-p3[0])**2 + (p1[1]-p3[1])**2);
            const s = (a+b+c)/2;
            return Math.sqrt(s*(s-a)*(s-b)*(s-c));
        }

        let triangularization2D = [];
        let totalArea = 0;
        for(let i = 2; i < this.points.length; i++){
            // [p1, p2, p3, area]
            let triangle = [points2D[0], points2D[i], points2D[i-1]];
            let area = triangleArea(...triangle);
            totalArea += area;
            triangle.push(area);
            triangularization2D.push(triangle);
        }
        this.area = totalArea;
        this.triangularization2D = triangularization2D;
    }

    static uniform(pts, orientationVector, color, lightSource){
        let properties = {
            color: color,
            lightSource:lightSource
        };

        return new PolyPlane(pts,orientationVector,properties);
    }

    // generate random point from surface of polygon along with direction
    sample(){
        // 1. find coordinates of polygon in 2d coordinates on plane
        // 2. decompose polygon into triangles
        // 3. pick a (weighted) random triangle 
        // 4. sample uniformly from chosen triangle via https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
        // 5. put back in 3d
        
        // 1 and 2 are done in the constructor
        // pick random triangle 
        let r = Math.random();
        let cumulativeArea = 0;
        let chosenTriangle = undefined;
        for(let t of this.triangularization2D){
            let normedArea = t[3] / this.area;
            if(r < cumulativeArea + normedArea){
                chosenTriangle = t;
                break;
            }
            cumulativeArea += normedArea;
        }
        if(chosenTriangle === undefined){
            console.log('ERROR. random triangle algorithm not working.');
        }

        // now pick a random point in our triangle 
        let r1 = Math.random(), r2 = Math.random();
        let f1 = 1 - r1**(1/2);
        let f2 = r1**(1/2)*(1-r2);
        let f3 = r2*r1**(1/2);
        let randomPt = [
            f1*chosenTriangle[0][0] + f2*chosenTriangle[1][0] + f3*chosenTriangle[2][0],
            f1*chosenTriangle[0][1] + f2*chosenTriangle[1][1] + f3*chosenTriangle[2][1]
        ];
        // now express in regular coordinates 
        let realRandomPt = plus(scale(randomPt[0], this.basis[0]), scale(randomPt[1], this.basis[1]));
        return {pt: realRandomPt, dir: BRDF.gaussDir(this.normal, this.properties.brdf.properties.lobeSize, this.normal), color:this.properties.color};
    }

    deserialize(){
        if(this.properties.brdf){
            this.properties.brdf.__proto__ = BRDF.prototype;
        }
    }

    intersect(ray){
        // plane-ray intersection followed by planar polygon intersection
        let den = dot(ray.direction, this.normal);
        const eps = 1e-6;
        if(Math.abs(den) < eps){
            // this plane is infinitely thin
            return [];
        }

        let t = dot(this.normal, minus(this.points[0], ray.start)) / den;

        // NOTE: setting t <= 0 produces rings on flat surfaces
        if(t <= eps){
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
            toReturn.surface = {
                normal: thisRef.normal,
                color: thisRef.properties.color
            }

            return toReturn;
        }

        if(count % 2 === 1){
            // within the polygon
            return [intersectObject(pt, this)];
        }
        return [];
    }
}

// note that featureCollection includes the light sources! (i.e. lightCollection \subset featureCollection)
function traceRay(ray, featureCollection, lightCollection, environment, maxDepth, flags){
    const debug = flags && flags === 'debug';

    let path = [];
    let futureRay = ray;
    while(path.length < maxDepth){
        futureRay.start = plus(futureRay.start, scale(1e-8, futureRay.direction));
        const collisions = featureCollection.intersect(futureRay);
        if(collisions.length === 0){
            path.push({
                color:[0,0,0],
                pos:[0,0,0] // doesn't matter where it is
            });
            break;
        }

        const collision = collisions[0];
        if(debug){
            console.log(collision.pos, collision.surface.normal, futureRay.direction);
        }

        const pathObj = {
            object:collision.object,
            pos: collision.pos,
            normal: collision.surface.normal,
            color: collision.surface.color,
            inDir: futureRay.direction
        };
        path.push(pathObj);

        futureRay = new Ray(pathObj.pos, pathObj.object.properties.brdf.sample(pathObj.inDir, pathObj.normal));
        pathObj.outDir = futureRay.direction;
    }
    let inSpectralRadiance = [0,0,0]; // in [r,g,b], each coordinate independent
    for(let k = path.length - 1; k >= 0; k--){
        // we wish to compute the spectral radiance in the outward direction by the rendering equation
        const workingNode = path[k];
        let emissivity = [0,0,0];
        // i.e. if the node is not an intersection but an endpoint (hit nothing)
        if(workingNode.object === undefined){
            continue;
        }
        if(workingNode.object.properties.lightSource){
            emissivity = scale(workingNode.object.properties.intensity/255, workingNode.color);
        }

        // brdf term cancels (as pdf function)
        let outSpectralRadiance = plus(emissivity, scale(Math.abs(dot(workingNode.outDir, workingNode.normal)),inSpectralRadiance));
        // apply absorbance via color of surface
        for(let i = 0 ; i <= 3; i++){
            outSpectralRadiance[i] *= workingNode.color[i]/255;
        }
        if(k === 0){
            // convert the spectral radiance to a color
            return outSpectralRadiance;
        }
        inSpectralRadiance = outSpectralRadiance;
        incomingSource = workingNode.pos;
    }
    return [0,0,0];
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

function snell(inVec, normal, refractionAngle){
    const theta = angleBetweenVecs(inVec, normal);
    return rotateAboutAxis(inVec, normalize(cross(inVec, normal)), refractionAngle - theta);
}

function rotateAboutAxis(vec, unitAxis, theta){
    // see https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
    let [x,y,z] = vec;
    let [u,v,w] = unitAxis;
    return [
        u*(u*x + v*y + w*z)*(1 - Math.cos(theta)) + x*Math.cos(theta) + (-w*y + v*z)*Math.sin(theta),
        v*(u*x + v*y + w*z)*(1 - Math.cos(theta)) + y*Math.cos(theta) + (w*x - u*z)*Math.sin(theta),
        w*(u*x + v*y + w*z)*(1 - Math.cos(theta)) + z*Math.cos(theta) + (-v*z + u*y)*Math.sin(theta)
    ];
}

function matchTransformation(v, startVec, destVec){
    let axis = normalize(cross(startVec, destVec));
    let acuteAngle = angleBetweenVecs(startVec, destVec);
    return rotateAboutAxis(v, axis, acuteAngle);
}

// correctly accounts for floating point error
function angleBetweenVecs(a,b){
    let inner = dot(a,b) / (mag(a) * mag(b));
    if(inner > 1){
        inner = 1;
    }
    if(inner < 1){
        inner = -1;
    }
    return Math.acos(inner);
}

function perpVec(a){
    let testVecs = [[1,0,0], [0,0,1], [0,1,0]];
    let minIndex = 0;
    let minValue = Math.abs(dot(a, testVecs[0]));
    for(let i of [1,2]){
        let v = Math.abs(dot(a, testVecs[i]));
        if(v < minValue){
            minValue = v;
            minIndex = i;
        }
    }

    let vec = testVecs[minIndex];
    return cross(a, vec);
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

// inVec points into the surface 
function reflect(inVec, normal){
    return minus(inVec, scale(2, project(inVec, normal)));
}

module.exports = {
    Sphere: Sphere,
    PolyPlane: PolyPlane,
    BRDF: BRDF,
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

},{"gauss-random":8}],2:[function(require,module,exports){

const Geometry = require('./geometry.js');
const PoissonDisk = require('poisson-disk-sampling')

module.exports = function(self){

    let featureCollection = null;
    let lightCollection = [];
    let environment = {};
    let basic = false;
    let camera = null;

    self.addEventListener('message',function(e){
        let type = e.data[0];
        if(type === 'overhead'){
            camera = e.data[1];
            featureCollection = e.data[2];
            lightCollection = e.data[3];
            environment = e.data[4];
            basic = e.data[5];
            // reattach methods that were stripped during serialization
            featureCollection.__proto__ = Geometry.FeatureCollection.prototype;
            lightCollection.__proto__ = Geometry.FeatureCollection.prototype;
            featureCollection.deserialize();
            lightCollection.deserialize();
            camera.__proto__ = Geometry.Camera.prototype;
        }else if(type === 'job'){
            // a set of locations to trace
            let locations = e.data[1];
            // work on pixels in batches and update periodically 

            // i.e. how many pixels to do before updating 
            const updateFreq = 200;
            let pixelBuffer = [];
            for(let i = 0; i < updateFreq; i++){
                pixelBuffer.push([0,0,0]);
            }
            let bufferIndex = 0;
            const pixWidth = 2 / camera.resolution;
            
            function sendBuffer(buffer){
                postMessage(['update', pixelBuffer.slice(0, bufferIndex)]);
                bufferIndex = 0;
            }
            for(let pix of locations){
                let [x,y,samples,confidence] = pix;
                // average spectral radiance
                let averageSR = [0,0,0];

                // generate poisson-distributed points 
                // let pds = new PoissonDisk([10,10], 10 / samples**(1/2));
                // let pts = pds.fill();
                // pts = pts.map(x => [pixWidth*x[0]/100, pixWidth*x[1]/100]);

                let pts = [ [ 5.623986229392337, 7.2601422312337505 ],
                    [ 6.644536871719595, 2.5696774186649627 ],
                    [ 0.9570739883176884, 3.081471711298274 ],
                    [ 0.6748887133929911, 8.444477688355523 ] ];
                pts = pts.map(x => [pixWidth*x[0]/10, pixWidth*x[1]/10]);
                let count = 0;
                for(let p of pts){
                    let pos = [x+p[0]-pixWidth/2, y+p[1]-pixWidth/2];
                    let localAverage = [0,0,0];
                    for(let j = 0; j < environment.samplesPerPoint; j++){
                        let ray = camera.getRay(pos);
                        // result spectral radiance
                        const resultSR = Geometry.traceRay(ray, featureCollection, lightCollection, environment, environment.rayDepth);
                        localAverage = Geometry.plus(localAverage, resultSR);
                    }
                    averageSR = Geometry.plus(averageSR, Geometry.scale(1/environment.samplesPerPoint, localAverage));
                    count++;
                }
                averageSR = Geometry.scale(1/count, averageSR);
                // convert to color
                averageSR = averageSR.map(x => 255*(2/(1+2**(-x/5))-1));
                // convert x and y to image coordinates 
                const realX = Math.floor((x+1)/2 * camera.resolution);
                const realY = Math.floor((y+1)/2 * camera.resolution);
                pixelBuffer[bufferIndex][0] = realX;
                pixelBuffer[bufferIndex][1] = realY;
                pixelBuffer[bufferIndex][2] = averageSR;
                bufferIndex++;

                if(bufferIndex >= updateFreq){
                    // update to main thread 
                    sendBuffer(pixelBuffer);
                }
            }
            sendBuffer(pixelBuffer);
            postMessage(['finished']);
        }
    });
}
},{"./geometry.js":1,"poisson-disk-sampling":13}],3:[function(require,module,exports){
const Geometry = require('./geometry.js');
const work = require('webworkify');
let colorConvert = require('color-convert');

/*
    ATTENTION! Hey Alan, because I know you'll forget, here's what you put in the command line 
    to get browserify to auto compile:
        watchify tracer.js -o bundle.js -v
        (arch) ~/node_modules/.bin/watchify tracer.js -o bundle.js -v
    Cheers! (You're welcome)

    And for you other people who don't remember, cd to the ray_tracing folder and type 
        python3 -m http.server
    and go to localhost [whatever it says]
*/

let featureCollection = null;
let lightCollection = null;
let environment = {
    ambientConstant:.4,
    diffuseConstant:.5,
    specularStrengthConstant:2,
    specularNarrownessConstant:80,
    maxRecursion:5,
    samplesPerPoint:5000,
    medium:{
        indexOfRefraction:1,
        opacity:.000,
        color:[54, 110, 169]
    },
    rayDepth:5,
    rayBranching:1
}
let camera = null;

const brightnessTransform = x => x;

// body onload
window.init = function(){
    initFeatures();
    initClickDebug();
    trace(image => {
        // do nothing
    }, false, 'canvas');
}

function initFeatures(){
    const features = [];
    const lightFeatures  = [];
    const whiteColor = [220,220,220];

    // floor
    const floorHeight = .1;
    const floor = Geometry.PolyPlane.uniform(
        [[-5, floorHeight, -1], [-5, floorHeight, 3], [5, floorHeight, 3], [5, floorHeight, -1]], // points
        [0,-1,0], // normal
        whiteColor, // color
        false // is light source
    );
    floor.properties.brdf = Geometry.BRDF.variableGlossy(1);

    // back wall 
    const wallWidth = 1.5;
    const wallHeight = 1;
    const wall1 = Geometry.PolyPlane.uniform(
        [[-wallWidth/2, floorHeight, 2], [wallWidth/2, floorHeight, 2], [wallWidth/2, floorHeight-wallHeight, 2], [-wallWidth/2, floorHeight-wallHeight]],
        [0,0,-1], // normal
        whiteColor, //color 
        false // is light source
    );
    wall1.properties.brdf = Geometry.BRDF.variableGlossy(1);

    // right wall
    const wallX = wallWidth;
    const wall2 = Geometry.PolyPlane.uniform(
        [[wallX, floorHeight, -1], [wallX, floorHeight, 3], [wallX, floorHeight - wallHeight, 3], [wallX, floorHeight-wallHeight, -1]],
        [-1, 0,0], // normal 
        whiteColor, // color
        false // is light source
    );
    wall2.properties.brdf = Geometry.BRDF.variableGlossy(1);

    // ceiling 
    const ceiling = Geometry.PolyPlane.uniform(
        [[wallX, floorHeight-wallHeight, -1],[-wallX, floorHeight-wallHeight,-1],[-wallX,floorHeight-wallHeight,3],[wallX,floorHeight-wallHeight,3]],
        [0,1,0], // normal
        whiteColor, // color 
        false // is light source
    );
    ceiling.properties.brdf = Geometry.BRDF.variableGlossy(1);

    // light source 
    const lightHeight = 1;
    const lightX = -.5
    const p1 = [lightX, floorHeight, 2];
    const p2 = [lightX, floorHeight, .5];
    const p3 = [lightX, floorHeight - lightHeight, .5];
    const p4 = [lightX, floorHeight - lightHeight, 2];
    const normal = Geometry.normalize(Geometry.cross(Geometry.minus(p2, p1), Geometry.minus(p3, p1)));
    let light = Geometry.PolyPlane.uniform(
        [p1, p2, p3, p4], // points
        normal, // normal
        whiteColor, // color
        true // is light source
    );

    //light = Geometry.Sphere.uniform(p1, .1, [255,255,255], true);
    light.properties.brdf = Geometry.BRDF.variableGlossy(1);
    light.properties.intensity = 200;

    // ball on floor 
    const ball1R = .12;
    const ball1 = Geometry.Sphere.uniform([0,floorHeight-ball1R, 1], ball1R, [255,20,20], false);
    ball1.properties.brdf = Geometry.BRDF.variableGlossy(1);

    features.push(floor, light, ball1, wall1, wall2, ceiling);
    lightFeatures.push(light);

    featureCollection = new Geometry.FeatureCollection(features);
    lightCollection = new Geometry.FeatureCollection(lightFeatures);
}

// click anywhere on the picture to send a ray with info showing up on the console
function initClickDebug(){
    const canvas = document.getElementById('canvas');
    let width = canvas.clientWidth;
    let height = canvas.clientHeight;
    canvas.addEventListener('click', function(e){
        let pos = [e.clientX, e.clientY];
        let halfWidth = camera.resolution / 2;
        let normPos = [(pos[0]-halfWidth)/halfWidth, (pos[1]-halfWidth)/halfWidth];
        let ray = camera.getRay(normPos);
        Geometry.traceRay(ray, featureCollection, lightCollection, environment, environment.maxRecursion, 'debug');
    });
}

function trace(traceFinished, basic, imageID){
    const res = 600;
    // this makes ~1.2m focus distance
    camera = new Geometry.Camera([0,0,0], [[1,0,0],[0,1,0],[0,0,1]], 35/1000, 37/1000, 30/1000, .01/1000, res);
    let image = fillArray(res,res,0);

    let workerJobs = []
    // 8 logical cores = 7 workers + main thread
    const threads = 7;
    for(let k = 0; k < threads; k++){
        workerJobs.push([]);
    }
    
    let workerIndex = 0;
    for(let i = -1; i <= 1; i += 2/res){
        for(let j = -1; j <= 1; j += 2/res){
            // assign pixels to workers
            // the third item is samples per pixel
            // the fourth item in the list is the required convergence,
            // essentially the confidence in the pixel necessary to finish
            workerJobs[workerIndex].push([i,j,5,.1]);
            workerIndex = (workerIndex + 1) % threads;
        }
    }

    // for fun shuffle all the lists so we get a pixelated effect
    // for(let i = 0; i < workerJobs.length; i++){
    //   workerJobs[i] = shuffle(workerJobs[i]);
    //}

    // let the workers run wild
    let workers = []

    let onfinish = () => {};

    let finishedWorkers = 0;
    for(let k = 0; k < workerJobs.length; k++){
        let newWorker = work(require('./rayWorker.js'));
        newWorker.postMessage(['overhead', camera, featureCollection, lightCollection, environment, basic]);
        newWorker.postMessage(['job', workerJobs[k]]);
        (function(){
            newWorker.addEventListener('message', function(e){
                let messageType = e.data[0];
                if(messageType === 'update'){
                    let pixelBuffer = e.data[1];
                    drawPixels(pixelBuffer, imageID);
                }else if(messageType === 'finished'){
                    newWorker.terminate();
                    finishedWorkers++;
                    console.log('finished workers: '+finishedWorkers)
                    if(finishedWorkers === threads){
                        // we're done
                        onfinish();
                    }
                }
            });
        })();
    }
}

function drawPixels(pixelBuffer, id){
    let canvas = document.getElementById(id);
    let ctx = canvas.getContext('2d');
    for(let p of pixelBuffer){
        let [x,y,color] = p;
        // we want to make our environment brighter 
        // convert to hsl
        let hsl = colorConvert.rgb.hsl.raw(...color);
        // apply a gamma function to lightness 
        hsl[2] = brightnessTransform(hsl[2]/100) * 100;
        // convert back to rgb
        color = colorConvert.hsl.rgb.raw(...hsl);

        ctx.fillStyle = 'rgb('+color[0]+','+color[1]+','+color[2]+')';
        ctx.fillRect(x,y,1,1);
    }
}

function fillArray(width, height, item){
    let image = [];
    for(let i = 0; i < height; i++){
        let newRow = [];
        for(let j = 0; j < width; j++){
            newRow.push(item);
        }
        image.push(newRow);
    }
    return image;
}

function shuffle(a) {
    for (let i = a.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [a[i], a[j]] = [a[j], a[i]];
    }
    return a;
}
},{"./geometry.js":1,"./rayWorker.js":2,"color-convert":5,"webworkify":15}],4:[function(require,module,exports){
/* MIT license */
var cssKeywords = require('color-name');

// NOTE: conversions should only return primitive values (i.e. arrays, or
//       values that give correct `typeof` results).
//       do not use box values types (i.e. Number(), String(), etc.)

var reverseKeywords = {};
for (var key in cssKeywords) {
	if (cssKeywords.hasOwnProperty(key)) {
		reverseKeywords[cssKeywords[key]] = key;
	}
}

var convert = module.exports = {
	rgb: {channels: 3, labels: 'rgb'},
	hsl: {channels: 3, labels: 'hsl'},
	hsv: {channels: 3, labels: 'hsv'},
	hwb: {channels: 3, labels: 'hwb'},
	cmyk: {channels: 4, labels: 'cmyk'},
	xyz: {channels: 3, labels: 'xyz'},
	lab: {channels: 3, labels: 'lab'},
	lch: {channels: 3, labels: 'lch'},
	hex: {channels: 1, labels: ['hex']},
	keyword: {channels: 1, labels: ['keyword']},
	ansi16: {channels: 1, labels: ['ansi16']},
	ansi256: {channels: 1, labels: ['ansi256']},
	hcg: {channels: 3, labels: ['h', 'c', 'g']},
	apple: {channels: 3, labels: ['r16', 'g16', 'b16']},
	gray: {channels: 1, labels: ['gray']}
};

// hide .channels and .labels properties
for (var model in convert) {
	if (convert.hasOwnProperty(model)) {
		if (!('channels' in convert[model])) {
			throw new Error('missing channels property: ' + model);
		}

		if (!('labels' in convert[model])) {
			throw new Error('missing channel labels property: ' + model);
		}

		if (convert[model].labels.length !== convert[model].channels) {
			throw new Error('channel and label counts mismatch: ' + model);
		}

		var channels = convert[model].channels;
		var labels = convert[model].labels;
		delete convert[model].channels;
		delete convert[model].labels;
		Object.defineProperty(convert[model], 'channels', {value: channels});
		Object.defineProperty(convert[model], 'labels', {value: labels});
	}
}

convert.rgb.hsl = function (rgb) {
	var r = rgb[0] / 255;
	var g = rgb[1] / 255;
	var b = rgb[2] / 255;
	var min = Math.min(r, g, b);
	var max = Math.max(r, g, b);
	var delta = max - min;
	var h;
	var s;
	var l;

	if (max === min) {
		h = 0;
	} else if (r === max) {
		h = (g - b) / delta;
	} else if (g === max) {
		h = 2 + (b - r) / delta;
	} else if (b === max) {
		h = 4 + (r - g) / delta;
	}

	h = Math.min(h * 60, 360);

	if (h < 0) {
		h += 360;
	}

	l = (min + max) / 2;

	if (max === min) {
		s = 0;
	} else if (l <= 0.5) {
		s = delta / (max + min);
	} else {
		s = delta / (2 - max - min);
	}

	return [h, s * 100, l * 100];
};

convert.rgb.hsv = function (rgb) {
	var rdif;
	var gdif;
	var bdif;
	var h;
	var s;

	var r = rgb[0] / 255;
	var g = rgb[1] / 255;
	var b = rgb[2] / 255;
	var v = Math.max(r, g, b);
	var diff = v - Math.min(r, g, b);
	var diffc = function (c) {
		return (v - c) / 6 / diff + 1 / 2;
	};

	if (diff === 0) {
		h = s = 0;
	} else {
		s = diff / v;
		rdif = diffc(r);
		gdif = diffc(g);
		bdif = diffc(b);

		if (r === v) {
			h = bdif - gdif;
		} else if (g === v) {
			h = (1 / 3) + rdif - bdif;
		} else if (b === v) {
			h = (2 / 3) + gdif - rdif;
		}
		if (h < 0) {
			h += 1;
		} else if (h > 1) {
			h -= 1;
		}
	}

	return [
		h * 360,
		s * 100,
		v * 100
	];
};

convert.rgb.hwb = function (rgb) {
	var r = rgb[0];
	var g = rgb[1];
	var b = rgb[2];
	var h = convert.rgb.hsl(rgb)[0];
	var w = 1 / 255 * Math.min(r, Math.min(g, b));

	b = 1 - 1 / 255 * Math.max(r, Math.max(g, b));

	return [h, w * 100, b * 100];
};

convert.rgb.cmyk = function (rgb) {
	var r = rgb[0] / 255;
	var g = rgb[1] / 255;
	var b = rgb[2] / 255;
	var c;
	var m;
	var y;
	var k;

	k = Math.min(1 - r, 1 - g, 1 - b);
	c = (1 - r - k) / (1 - k) || 0;
	m = (1 - g - k) / (1 - k) || 0;
	y = (1 - b - k) / (1 - k) || 0;

	return [c * 100, m * 100, y * 100, k * 100];
};

/**
 * See https://en.m.wikipedia.org/wiki/Euclidean_distance#Squared_Euclidean_distance
 * */
function comparativeDistance(x, y) {
	return (
		Math.pow(x[0] - y[0], 2) +
		Math.pow(x[1] - y[1], 2) +
		Math.pow(x[2] - y[2], 2)
	);
}

convert.rgb.keyword = function (rgb) {
	var reversed = reverseKeywords[rgb];
	if (reversed) {
		return reversed;
	}

	var currentClosestDistance = Infinity;
	var currentClosestKeyword;

	for (var keyword in cssKeywords) {
		if (cssKeywords.hasOwnProperty(keyword)) {
			var value = cssKeywords[keyword];

			// Compute comparative distance
			var distance = comparativeDistance(rgb, value);

			// Check if its less, if so set as closest
			if (distance < currentClosestDistance) {
				currentClosestDistance = distance;
				currentClosestKeyword = keyword;
			}
		}
	}

	return currentClosestKeyword;
};

convert.keyword.rgb = function (keyword) {
	return cssKeywords[keyword];
};

convert.rgb.xyz = function (rgb) {
	var r = rgb[0] / 255;
	var g = rgb[1] / 255;
	var b = rgb[2] / 255;

	// assume sRGB
	r = r > 0.04045 ? Math.pow(((r + 0.055) / 1.055), 2.4) : (r / 12.92);
	g = g > 0.04045 ? Math.pow(((g + 0.055) / 1.055), 2.4) : (g / 12.92);
	b = b > 0.04045 ? Math.pow(((b + 0.055) / 1.055), 2.4) : (b / 12.92);

	var x = (r * 0.4124) + (g * 0.3576) + (b * 0.1805);
	var y = (r * 0.2126) + (g * 0.7152) + (b * 0.0722);
	var z = (r * 0.0193) + (g * 0.1192) + (b * 0.9505);

	return [x * 100, y * 100, z * 100];
};

convert.rgb.lab = function (rgb) {
	var xyz = convert.rgb.xyz(rgb);
	var x = xyz[0];
	var y = xyz[1];
	var z = xyz[2];
	var l;
	var a;
	var b;

	x /= 95.047;
	y /= 100;
	z /= 108.883;

	x = x > 0.008856 ? Math.pow(x, 1 / 3) : (7.787 * x) + (16 / 116);
	y = y > 0.008856 ? Math.pow(y, 1 / 3) : (7.787 * y) + (16 / 116);
	z = z > 0.008856 ? Math.pow(z, 1 / 3) : (7.787 * z) + (16 / 116);

	l = (116 * y) - 16;
	a = 500 * (x - y);
	b = 200 * (y - z);

	return [l, a, b];
};

convert.hsl.rgb = function (hsl) {
	var h = hsl[0] / 360;
	var s = hsl[1] / 100;
	var l = hsl[2] / 100;
	var t1;
	var t2;
	var t3;
	var rgb;
	var val;

	if (s === 0) {
		val = l * 255;
		return [val, val, val];
	}

	if (l < 0.5) {
		t2 = l * (1 + s);
	} else {
		t2 = l + s - l * s;
	}

	t1 = 2 * l - t2;

	rgb = [0, 0, 0];
	for (var i = 0; i < 3; i++) {
		t3 = h + 1 / 3 * -(i - 1);
		if (t3 < 0) {
			t3++;
		}
		if (t3 > 1) {
			t3--;
		}

		if (6 * t3 < 1) {
			val = t1 + (t2 - t1) * 6 * t3;
		} else if (2 * t3 < 1) {
			val = t2;
		} else if (3 * t3 < 2) {
			val = t1 + (t2 - t1) * (2 / 3 - t3) * 6;
		} else {
			val = t1;
		}

		rgb[i] = val * 255;
	}

	return rgb;
};

convert.hsl.hsv = function (hsl) {
	var h = hsl[0];
	var s = hsl[1] / 100;
	var l = hsl[2] / 100;
	var smin = s;
	var lmin = Math.max(l, 0.01);
	var sv;
	var v;

	l *= 2;
	s *= (l <= 1) ? l : 2 - l;
	smin *= lmin <= 1 ? lmin : 2 - lmin;
	v = (l + s) / 2;
	sv = l === 0 ? (2 * smin) / (lmin + smin) : (2 * s) / (l + s);

	return [h, sv * 100, v * 100];
};

convert.hsv.rgb = function (hsv) {
	var h = hsv[0] / 60;
	var s = hsv[1] / 100;
	var v = hsv[2] / 100;
	var hi = Math.floor(h) % 6;

	var f = h - Math.floor(h);
	var p = 255 * v * (1 - s);
	var q = 255 * v * (1 - (s * f));
	var t = 255 * v * (1 - (s * (1 - f)));
	v *= 255;

	switch (hi) {
		case 0:
			return [v, t, p];
		case 1:
			return [q, v, p];
		case 2:
			return [p, v, t];
		case 3:
			return [p, q, v];
		case 4:
			return [t, p, v];
		case 5:
			return [v, p, q];
	}
};

convert.hsv.hsl = function (hsv) {
	var h = hsv[0];
	var s = hsv[1] / 100;
	var v = hsv[2] / 100;
	var vmin = Math.max(v, 0.01);
	var lmin;
	var sl;
	var l;

	l = (2 - s) * v;
	lmin = (2 - s) * vmin;
	sl = s * vmin;
	sl /= (lmin <= 1) ? lmin : 2 - lmin;
	sl = sl || 0;
	l /= 2;

	return [h, sl * 100, l * 100];
};

// http://dev.w3.org/csswg/css-color/#hwb-to-rgb
convert.hwb.rgb = function (hwb) {
	var h = hwb[0] / 360;
	var wh = hwb[1] / 100;
	var bl = hwb[2] / 100;
	var ratio = wh + bl;
	var i;
	var v;
	var f;
	var n;

	// wh + bl cant be > 1
	if (ratio > 1) {
		wh /= ratio;
		bl /= ratio;
	}

	i = Math.floor(6 * h);
	v = 1 - bl;
	f = 6 * h - i;

	if ((i & 0x01) !== 0) {
		f = 1 - f;
	}

	n = wh + f * (v - wh); // linear interpolation

	var r;
	var g;
	var b;
	switch (i) {
		default:
		case 6:
		case 0: r = v; g = n; b = wh; break;
		case 1: r = n; g = v; b = wh; break;
		case 2: r = wh; g = v; b = n; break;
		case 3: r = wh; g = n; b = v; break;
		case 4: r = n; g = wh; b = v; break;
		case 5: r = v; g = wh; b = n; break;
	}

	return [r * 255, g * 255, b * 255];
};

convert.cmyk.rgb = function (cmyk) {
	var c = cmyk[0] / 100;
	var m = cmyk[1] / 100;
	var y = cmyk[2] / 100;
	var k = cmyk[3] / 100;
	var r;
	var g;
	var b;

	r = 1 - Math.min(1, c * (1 - k) + k);
	g = 1 - Math.min(1, m * (1 - k) + k);
	b = 1 - Math.min(1, y * (1 - k) + k);

	return [r * 255, g * 255, b * 255];
};

convert.xyz.rgb = function (xyz) {
	var x = xyz[0] / 100;
	var y = xyz[1] / 100;
	var z = xyz[2] / 100;
	var r;
	var g;
	var b;

	r = (x * 3.2406) + (y * -1.5372) + (z * -0.4986);
	g = (x * -0.9689) + (y * 1.8758) + (z * 0.0415);
	b = (x * 0.0557) + (y * -0.2040) + (z * 1.0570);

	// assume sRGB
	r = r > 0.0031308
		? ((1.055 * Math.pow(r, 1.0 / 2.4)) - 0.055)
		: r * 12.92;

	g = g > 0.0031308
		? ((1.055 * Math.pow(g, 1.0 / 2.4)) - 0.055)
		: g * 12.92;

	b = b > 0.0031308
		? ((1.055 * Math.pow(b, 1.0 / 2.4)) - 0.055)
		: b * 12.92;

	r = Math.min(Math.max(0, r), 1);
	g = Math.min(Math.max(0, g), 1);
	b = Math.min(Math.max(0, b), 1);

	return [r * 255, g * 255, b * 255];
};

convert.xyz.lab = function (xyz) {
	var x = xyz[0];
	var y = xyz[1];
	var z = xyz[2];
	var l;
	var a;
	var b;

	x /= 95.047;
	y /= 100;
	z /= 108.883;

	x = x > 0.008856 ? Math.pow(x, 1 / 3) : (7.787 * x) + (16 / 116);
	y = y > 0.008856 ? Math.pow(y, 1 / 3) : (7.787 * y) + (16 / 116);
	z = z > 0.008856 ? Math.pow(z, 1 / 3) : (7.787 * z) + (16 / 116);

	l = (116 * y) - 16;
	a = 500 * (x - y);
	b = 200 * (y - z);

	return [l, a, b];
};

convert.lab.xyz = function (lab) {
	var l = lab[0];
	var a = lab[1];
	var b = lab[2];
	var x;
	var y;
	var z;

	y = (l + 16) / 116;
	x = a / 500 + y;
	z = y - b / 200;

	var y2 = Math.pow(y, 3);
	var x2 = Math.pow(x, 3);
	var z2 = Math.pow(z, 3);
	y = y2 > 0.008856 ? y2 : (y - 16 / 116) / 7.787;
	x = x2 > 0.008856 ? x2 : (x - 16 / 116) / 7.787;
	z = z2 > 0.008856 ? z2 : (z - 16 / 116) / 7.787;

	x *= 95.047;
	y *= 100;
	z *= 108.883;

	return [x, y, z];
};

convert.lab.lch = function (lab) {
	var l = lab[0];
	var a = lab[1];
	var b = lab[2];
	var hr;
	var h;
	var c;

	hr = Math.atan2(b, a);
	h = hr * 360 / 2 / Math.PI;

	if (h < 0) {
		h += 360;
	}

	c = Math.sqrt(a * a + b * b);

	return [l, c, h];
};

convert.lch.lab = function (lch) {
	var l = lch[0];
	var c = lch[1];
	var h = lch[2];
	var a;
	var b;
	var hr;

	hr = h / 360 * 2 * Math.PI;
	a = c * Math.cos(hr);
	b = c * Math.sin(hr);

	return [l, a, b];
};

convert.rgb.ansi16 = function (args) {
	var r = args[0];
	var g = args[1];
	var b = args[2];
	var value = 1 in arguments ? arguments[1] : convert.rgb.hsv(args)[2]; // hsv -> ansi16 optimization

	value = Math.round(value / 50);

	if (value === 0) {
		return 30;
	}

	var ansi = 30
		+ ((Math.round(b / 255) << 2)
		| (Math.round(g / 255) << 1)
		| Math.round(r / 255));

	if (value === 2) {
		ansi += 60;
	}

	return ansi;
};

convert.hsv.ansi16 = function (args) {
	// optimization here; we already know the value and don't need to get
	// it converted for us.
	return convert.rgb.ansi16(convert.hsv.rgb(args), args[2]);
};

convert.rgb.ansi256 = function (args) {
	var r = args[0];
	var g = args[1];
	var b = args[2];

	// we use the extended greyscale palette here, with the exception of
	// black and white. normal palette only has 4 greyscale shades.
	if (r === g && g === b) {
		if (r < 8) {
			return 16;
		}

		if (r > 248) {
			return 231;
		}

		return Math.round(((r - 8) / 247) * 24) + 232;
	}

	var ansi = 16
		+ (36 * Math.round(r / 255 * 5))
		+ (6 * Math.round(g / 255 * 5))
		+ Math.round(b / 255 * 5);

	return ansi;
};

convert.ansi16.rgb = function (args) {
	var color = args % 10;

	// handle greyscale
	if (color === 0 || color === 7) {
		if (args > 50) {
			color += 3.5;
		}

		color = color / 10.5 * 255;

		return [color, color, color];
	}

	var mult = (~~(args > 50) + 1) * 0.5;
	var r = ((color & 1) * mult) * 255;
	var g = (((color >> 1) & 1) * mult) * 255;
	var b = (((color >> 2) & 1) * mult) * 255;

	return [r, g, b];
};

convert.ansi256.rgb = function (args) {
	// handle greyscale
	if (args >= 232) {
		var c = (args - 232) * 10 + 8;
		return [c, c, c];
	}

	args -= 16;

	var rem;
	var r = Math.floor(args / 36) / 5 * 255;
	var g = Math.floor((rem = args % 36) / 6) / 5 * 255;
	var b = (rem % 6) / 5 * 255;

	return [r, g, b];
};

convert.rgb.hex = function (args) {
	var integer = ((Math.round(args[0]) & 0xFF) << 16)
		+ ((Math.round(args[1]) & 0xFF) << 8)
		+ (Math.round(args[2]) & 0xFF);

	var string = integer.toString(16).toUpperCase();
	return '000000'.substring(string.length) + string;
};

convert.hex.rgb = function (args) {
	var match = args.toString(16).match(/[a-f0-9]{6}|[a-f0-9]{3}/i);
	if (!match) {
		return [0, 0, 0];
	}

	var colorString = match[0];

	if (match[0].length === 3) {
		colorString = colorString.split('').map(function (char) {
			return char + char;
		}).join('');
	}

	var integer = parseInt(colorString, 16);
	var r = (integer >> 16) & 0xFF;
	var g = (integer >> 8) & 0xFF;
	var b = integer & 0xFF;

	return [r, g, b];
};

convert.rgb.hcg = function (rgb) {
	var r = rgb[0] / 255;
	var g = rgb[1] / 255;
	var b = rgb[2] / 255;
	var max = Math.max(Math.max(r, g), b);
	var min = Math.min(Math.min(r, g), b);
	var chroma = (max - min);
	var grayscale;
	var hue;

	if (chroma < 1) {
		grayscale = min / (1 - chroma);
	} else {
		grayscale = 0;
	}

	if (chroma <= 0) {
		hue = 0;
	} else
	if (max === r) {
		hue = ((g - b) / chroma) % 6;
	} else
	if (max === g) {
		hue = 2 + (b - r) / chroma;
	} else {
		hue = 4 + (r - g) / chroma + 4;
	}

	hue /= 6;
	hue %= 1;

	return [hue * 360, chroma * 100, grayscale * 100];
};

convert.hsl.hcg = function (hsl) {
	var s = hsl[1] / 100;
	var l = hsl[2] / 100;
	var c = 1;
	var f = 0;

	if (l < 0.5) {
		c = 2.0 * s * l;
	} else {
		c = 2.0 * s * (1.0 - l);
	}

	if (c < 1.0) {
		f = (l - 0.5 * c) / (1.0 - c);
	}

	return [hsl[0], c * 100, f * 100];
};

convert.hsv.hcg = function (hsv) {
	var s = hsv[1] / 100;
	var v = hsv[2] / 100;

	var c = s * v;
	var f = 0;

	if (c < 1.0) {
		f = (v - c) / (1 - c);
	}

	return [hsv[0], c * 100, f * 100];
};

convert.hcg.rgb = function (hcg) {
	var h = hcg[0] / 360;
	var c = hcg[1] / 100;
	var g = hcg[2] / 100;

	if (c === 0.0) {
		return [g * 255, g * 255, g * 255];
	}

	var pure = [0, 0, 0];
	var hi = (h % 1) * 6;
	var v = hi % 1;
	var w = 1 - v;
	var mg = 0;

	switch (Math.floor(hi)) {
		case 0:
			pure[0] = 1; pure[1] = v; pure[2] = 0; break;
		case 1:
			pure[0] = w; pure[1] = 1; pure[2] = 0; break;
		case 2:
			pure[0] = 0; pure[1] = 1; pure[2] = v; break;
		case 3:
			pure[0] = 0; pure[1] = w; pure[2] = 1; break;
		case 4:
			pure[0] = v; pure[1] = 0; pure[2] = 1; break;
		default:
			pure[0] = 1; pure[1] = 0; pure[2] = w;
	}

	mg = (1.0 - c) * g;

	return [
		(c * pure[0] + mg) * 255,
		(c * pure[1] + mg) * 255,
		(c * pure[2] + mg) * 255
	];
};

convert.hcg.hsv = function (hcg) {
	var c = hcg[1] / 100;
	var g = hcg[2] / 100;

	var v = c + g * (1.0 - c);
	var f = 0;

	if (v > 0.0) {
		f = c / v;
	}

	return [hcg[0], f * 100, v * 100];
};

convert.hcg.hsl = function (hcg) {
	var c = hcg[1] / 100;
	var g = hcg[2] / 100;

	var l = g * (1.0 - c) + 0.5 * c;
	var s = 0;

	if (l > 0.0 && l < 0.5) {
		s = c / (2 * l);
	} else
	if (l >= 0.5 && l < 1.0) {
		s = c / (2 * (1 - l));
	}

	return [hcg[0], s * 100, l * 100];
};

convert.hcg.hwb = function (hcg) {
	var c = hcg[1] / 100;
	var g = hcg[2] / 100;
	var v = c + g * (1.0 - c);
	return [hcg[0], (v - c) * 100, (1 - v) * 100];
};

convert.hwb.hcg = function (hwb) {
	var w = hwb[1] / 100;
	var b = hwb[2] / 100;
	var v = 1 - b;
	var c = v - w;
	var g = 0;

	if (c < 1) {
		g = (v - c) / (1 - c);
	}

	return [hwb[0], c * 100, g * 100];
};

convert.apple.rgb = function (apple) {
	return [(apple[0] / 65535) * 255, (apple[1] / 65535) * 255, (apple[2] / 65535) * 255];
};

convert.rgb.apple = function (rgb) {
	return [(rgb[0] / 255) * 65535, (rgb[1] / 255) * 65535, (rgb[2] / 255) * 65535];
};

convert.gray.rgb = function (args) {
	return [args[0] / 100 * 255, args[0] / 100 * 255, args[0] / 100 * 255];
};

convert.gray.hsl = convert.gray.hsv = function (args) {
	return [0, 0, args[0]];
};

convert.gray.hwb = function (gray) {
	return [0, 100, gray[0]];
};

convert.gray.cmyk = function (gray) {
	return [0, 0, 0, gray[0]];
};

convert.gray.lab = function (gray) {
	return [gray[0], 0, 0];
};

convert.gray.hex = function (gray) {
	var val = Math.round(gray[0] / 100 * 255) & 0xFF;
	var integer = (val << 16) + (val << 8) + val;

	var string = integer.toString(16).toUpperCase();
	return '000000'.substring(string.length) + string;
};

convert.rgb.gray = function (rgb) {
	var val = (rgb[0] + rgb[1] + rgb[2]) / 3;
	return [val / 255 * 100];
};

},{"color-name":7}],5:[function(require,module,exports){
var conversions = require('./conversions');
var route = require('./route');

var convert = {};

var models = Object.keys(conversions);

function wrapRaw(fn) {
	var wrappedFn = function (args) {
		if (args === undefined || args === null) {
			return args;
		}

		if (arguments.length > 1) {
			args = Array.prototype.slice.call(arguments);
		}

		return fn(args);
	};

	// preserve .conversion property if there is one
	if ('conversion' in fn) {
		wrappedFn.conversion = fn.conversion;
	}

	return wrappedFn;
}

function wrapRounded(fn) {
	var wrappedFn = function (args) {
		if (args === undefined || args === null) {
			return args;
		}

		if (arguments.length > 1) {
			args = Array.prototype.slice.call(arguments);
		}

		var result = fn(args);

		// we're assuming the result is an array here.
		// see notice in conversions.js; don't use box types
		// in conversion functions.
		if (typeof result === 'object') {
			for (var len = result.length, i = 0; i < len; i++) {
				result[i] = Math.round(result[i]);
			}
		}

		return result;
	};

	// preserve .conversion property if there is one
	if ('conversion' in fn) {
		wrappedFn.conversion = fn.conversion;
	}

	return wrappedFn;
}

models.forEach(function (fromModel) {
	convert[fromModel] = {};

	Object.defineProperty(convert[fromModel], 'channels', {value: conversions[fromModel].channels});
	Object.defineProperty(convert[fromModel], 'labels', {value: conversions[fromModel].labels});

	var routes = route(fromModel);
	var routeModels = Object.keys(routes);

	routeModels.forEach(function (toModel) {
		var fn = routes[toModel];

		convert[fromModel][toModel] = wrapRounded(fn);
		convert[fromModel][toModel].raw = wrapRaw(fn);
	});
});

module.exports = convert;

},{"./conversions":4,"./route":6}],6:[function(require,module,exports){
var conversions = require('./conversions');

/*
	this function routes a model to all other models.

	all functions that are routed have a property `.conversion` attached
	to the returned synthetic function. This property is an array
	of strings, each with the steps in between the 'from' and 'to'
	color models (inclusive).

	conversions that are not possible simply are not included.
*/

function buildGraph() {
	var graph = {};
	// https://jsperf.com/object-keys-vs-for-in-with-closure/3
	var models = Object.keys(conversions);

	for (var len = models.length, i = 0; i < len; i++) {
		graph[models[i]] = {
			// http://jsperf.com/1-vs-infinity
			// micro-opt, but this is simple.
			distance: -1,
			parent: null
		};
	}

	return graph;
}

// https://en.wikipedia.org/wiki/Breadth-first_search
function deriveBFS(fromModel) {
	var graph = buildGraph();
	var queue = [fromModel]; // unshift -> queue -> pop

	graph[fromModel].distance = 0;

	while (queue.length) {
		var current = queue.pop();
		var adjacents = Object.keys(conversions[current]);

		for (var len = adjacents.length, i = 0; i < len; i++) {
			var adjacent = adjacents[i];
			var node = graph[adjacent];

			if (node.distance === -1) {
				node.distance = graph[current].distance + 1;
				node.parent = current;
				queue.unshift(adjacent);
			}
		}
	}

	return graph;
}

function link(from, to) {
	return function (args) {
		return to(from(args));
	};
}

function wrapConversion(toModel, graph) {
	var path = [graph[toModel].parent, toModel];
	var fn = conversions[graph[toModel].parent][toModel];

	var cur = graph[toModel].parent;
	while (graph[cur].parent) {
		path.unshift(graph[cur].parent);
		fn = link(conversions[graph[cur].parent][cur], fn);
		cur = graph[cur].parent;
	}

	fn.conversion = path;
	return fn;
}

module.exports = function (fromModel) {
	var graph = deriveBFS(fromModel);
	var conversion = {};

	var models = Object.keys(graph);
	for (var len = models.length, i = 0; i < len; i++) {
		var toModel = models[i];
		var node = graph[toModel];

		if (node.parent === null) {
			// no possible conversion, or this node is the source model.
			continue;
		}

		conversion[toModel] = wrapConversion(toModel, graph);
	}

	return conversion;
};


},{"./conversions":4}],7:[function(require,module,exports){
'use strict'

module.exports = {
	"aliceblue": [240, 248, 255],
	"antiquewhite": [250, 235, 215],
	"aqua": [0, 255, 255],
	"aquamarine": [127, 255, 212],
	"azure": [240, 255, 255],
	"beige": [245, 245, 220],
	"bisque": [255, 228, 196],
	"black": [0, 0, 0],
	"blanchedalmond": [255, 235, 205],
	"blue": [0, 0, 255],
	"blueviolet": [138, 43, 226],
	"brown": [165, 42, 42],
	"burlywood": [222, 184, 135],
	"cadetblue": [95, 158, 160],
	"chartreuse": [127, 255, 0],
	"chocolate": [210, 105, 30],
	"coral": [255, 127, 80],
	"cornflowerblue": [100, 149, 237],
	"cornsilk": [255, 248, 220],
	"crimson": [220, 20, 60],
	"cyan": [0, 255, 255],
	"darkblue": [0, 0, 139],
	"darkcyan": [0, 139, 139],
	"darkgoldenrod": [184, 134, 11],
	"darkgray": [169, 169, 169],
	"darkgreen": [0, 100, 0],
	"darkgrey": [169, 169, 169],
	"darkkhaki": [189, 183, 107],
	"darkmagenta": [139, 0, 139],
	"darkolivegreen": [85, 107, 47],
	"darkorange": [255, 140, 0],
	"darkorchid": [153, 50, 204],
	"darkred": [139, 0, 0],
	"darksalmon": [233, 150, 122],
	"darkseagreen": [143, 188, 143],
	"darkslateblue": [72, 61, 139],
	"darkslategray": [47, 79, 79],
	"darkslategrey": [47, 79, 79],
	"darkturquoise": [0, 206, 209],
	"darkviolet": [148, 0, 211],
	"deeppink": [255, 20, 147],
	"deepskyblue": [0, 191, 255],
	"dimgray": [105, 105, 105],
	"dimgrey": [105, 105, 105],
	"dodgerblue": [30, 144, 255],
	"firebrick": [178, 34, 34],
	"floralwhite": [255, 250, 240],
	"forestgreen": [34, 139, 34],
	"fuchsia": [255, 0, 255],
	"gainsboro": [220, 220, 220],
	"ghostwhite": [248, 248, 255],
	"gold": [255, 215, 0],
	"goldenrod": [218, 165, 32],
	"gray": [128, 128, 128],
	"green": [0, 128, 0],
	"greenyellow": [173, 255, 47],
	"grey": [128, 128, 128],
	"honeydew": [240, 255, 240],
	"hotpink": [255, 105, 180],
	"indianred": [205, 92, 92],
	"indigo": [75, 0, 130],
	"ivory": [255, 255, 240],
	"khaki": [240, 230, 140],
	"lavender": [230, 230, 250],
	"lavenderblush": [255, 240, 245],
	"lawngreen": [124, 252, 0],
	"lemonchiffon": [255, 250, 205],
	"lightblue": [173, 216, 230],
	"lightcoral": [240, 128, 128],
	"lightcyan": [224, 255, 255],
	"lightgoldenrodyellow": [250, 250, 210],
	"lightgray": [211, 211, 211],
	"lightgreen": [144, 238, 144],
	"lightgrey": [211, 211, 211],
	"lightpink": [255, 182, 193],
	"lightsalmon": [255, 160, 122],
	"lightseagreen": [32, 178, 170],
	"lightskyblue": [135, 206, 250],
	"lightslategray": [119, 136, 153],
	"lightslategrey": [119, 136, 153],
	"lightsteelblue": [176, 196, 222],
	"lightyellow": [255, 255, 224],
	"lime": [0, 255, 0],
	"limegreen": [50, 205, 50],
	"linen": [250, 240, 230],
	"magenta": [255, 0, 255],
	"maroon": [128, 0, 0],
	"mediumaquamarine": [102, 205, 170],
	"mediumblue": [0, 0, 205],
	"mediumorchid": [186, 85, 211],
	"mediumpurple": [147, 112, 219],
	"mediumseagreen": [60, 179, 113],
	"mediumslateblue": [123, 104, 238],
	"mediumspringgreen": [0, 250, 154],
	"mediumturquoise": [72, 209, 204],
	"mediumvioletred": [199, 21, 133],
	"midnightblue": [25, 25, 112],
	"mintcream": [245, 255, 250],
	"mistyrose": [255, 228, 225],
	"moccasin": [255, 228, 181],
	"navajowhite": [255, 222, 173],
	"navy": [0, 0, 128],
	"oldlace": [253, 245, 230],
	"olive": [128, 128, 0],
	"olivedrab": [107, 142, 35],
	"orange": [255, 165, 0],
	"orangered": [255, 69, 0],
	"orchid": [218, 112, 214],
	"palegoldenrod": [238, 232, 170],
	"palegreen": [152, 251, 152],
	"paleturquoise": [175, 238, 238],
	"palevioletred": [219, 112, 147],
	"papayawhip": [255, 239, 213],
	"peachpuff": [255, 218, 185],
	"peru": [205, 133, 63],
	"pink": [255, 192, 203],
	"plum": [221, 160, 221],
	"powderblue": [176, 224, 230],
	"purple": [128, 0, 128],
	"rebeccapurple": [102, 51, 153],
	"red": [255, 0, 0],
	"rosybrown": [188, 143, 143],
	"royalblue": [65, 105, 225],
	"saddlebrown": [139, 69, 19],
	"salmon": [250, 128, 114],
	"sandybrown": [244, 164, 96],
	"seagreen": [46, 139, 87],
	"seashell": [255, 245, 238],
	"sienna": [160, 82, 45],
	"silver": [192, 192, 192],
	"skyblue": [135, 206, 235],
	"slateblue": [106, 90, 205],
	"slategray": [112, 128, 144],
	"slategrey": [112, 128, 144],
	"snow": [255, 250, 250],
	"springgreen": [0, 255, 127],
	"steelblue": [70, 130, 180],
	"tan": [210, 180, 140],
	"teal": [0, 128, 128],
	"thistle": [216, 191, 216],
	"tomato": [255, 99, 71],
	"turquoise": [64, 224, 208],
	"violet": [238, 130, 238],
	"wheat": [245, 222, 179],
	"white": [255, 255, 255],
	"whitesmoke": [245, 245, 245],
	"yellow": [255, 255, 0],
	"yellowgreen": [154, 205, 50]
};

},{}],8:[function(require,module,exports){
'use strict'

module.exports = gaussRandom

function gaussRandom() {
  return Math.sqrt(-2.0 * Math.log(Math.random())) * Math.cos(2.0 * Math.PI * Math.random())
}
},{}],9:[function(require,module,exports){
"use strict"

function iota(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = i
  }
  return result
}

module.exports = iota
},{}],10:[function(require,module,exports){
/*!
 * Determine if an object is a Buffer
 *
 * @author   Feross Aboukhadijeh <https://feross.org>
 * @license  MIT
 */

// The _isBuffer check is for Safari 5-7 support, because it's missing
// Object.prototype.constructor. Remove this eventually
module.exports = function (obj) {
  return obj != null && (isBuffer(obj) || isSlowBuffer(obj) || !!obj._isBuffer)
}

function isBuffer (obj) {
  return !!obj.constructor && typeof obj.constructor.isBuffer === 'function' && obj.constructor.isBuffer(obj)
}

// For Node v0.10 support. Remove this eventually.
function isSlowBuffer (obj) {
  return typeof obj.readFloatLE === 'function' && typeof obj.slice === 'function' && isBuffer(obj.slice(0, 0))
}

},{}],11:[function(require,module,exports){
module.exports = function moore(range, dimensions) {
  range = range || 1
  dimensions = dimensions || 2

  var size = range * 2 + 1
  var length = Math.pow(size, dimensions) - 1
  var neighbors = new Array(length)

  for (var i = 0; i < length; i++) {
    var neighbor = neighbors[i] = new Array(dimensions)
    var index = i < length / 2 ? i : i + 1
    for (var dimension = 1; dimension <= dimensions; dimension++) {
      var value = index % Math.pow(size, dimension)
      neighbor[dimension - 1] = value / Math.pow(size, dimension - 1) - range
      index -= value
    }
  }

  return neighbors
}

},{}],12:[function(require,module,exports){
var iota = require("iota-array")
var isBuffer = require("is-buffer")

var hasTypedArrays  = ((typeof Float64Array) !== "undefined")

function compare1st(a, b) {
  return a[0] - b[0]
}

function order() {
  var stride = this.stride
  var terms = new Array(stride.length)
  var i
  for(i=0; i<terms.length; ++i) {
    terms[i] = [Math.abs(stride[i]), i]
  }
  terms.sort(compare1st)
  var result = new Array(terms.length)
  for(i=0; i<result.length; ++i) {
    result[i] = terms[i][1]
  }
  return result
}

function compileConstructor(dtype, dimension) {
  var className = ["View", dimension, "d", dtype].join("")
  if(dimension < 0) {
    className = "View_Nil" + dtype
  }
  var useGetters = (dtype === "generic")

  if(dimension === -1) {
    //Special case for trivial arrays
    var code =
      "function "+className+"(a){this.data=a;};\
var proto="+className+".prototype;\
proto.dtype='"+dtype+"';\
proto.index=function(){return -1};\
proto.size=0;\
proto.dimension=-1;\
proto.shape=proto.stride=proto.order=[];\
proto.lo=proto.hi=proto.transpose=proto.step=\
function(){return new "+className+"(this.data);};\
proto.get=proto.set=function(){};\
proto.pick=function(){return null};\
return function construct_"+className+"(a){return new "+className+"(a);}"
    var procedure = new Function(code)
    return procedure()
  } else if(dimension === 0) {
    //Special case for 0d arrays
    var code =
      "function "+className+"(a,d) {\
this.data = a;\
this.offset = d\
};\
var proto="+className+".prototype;\
proto.dtype='"+dtype+"';\
proto.index=function(){return this.offset};\
proto.dimension=0;\
proto.size=1;\
proto.shape=\
proto.stride=\
proto.order=[];\
proto.lo=\
proto.hi=\
proto.transpose=\
proto.step=function "+className+"_copy() {\
return new "+className+"(this.data,this.offset)\
};\
proto.pick=function "+className+"_pick(){\
return TrivialArray(this.data);\
};\
proto.valueOf=proto.get=function "+className+"_get(){\
return "+(useGetters ? "this.data.get(this.offset)" : "this.data[this.offset]")+
"};\
proto.set=function "+className+"_set(v){\
return "+(useGetters ? "this.data.set(this.offset,v)" : "this.data[this.offset]=v")+"\
};\
return function construct_"+className+"(a,b,c,d){return new "+className+"(a,d)}"
    var procedure = new Function("TrivialArray", code)
    return procedure(CACHED_CONSTRUCTORS[dtype][0])
  }

  var code = ["'use strict'"]

  //Create constructor for view
  var indices = iota(dimension)
  var args = indices.map(function(i) { return "i"+i })
  var index_str = "this.offset+" + indices.map(function(i) {
        return "this.stride[" + i + "]*i" + i
      }).join("+")
  var shapeArg = indices.map(function(i) {
      return "b"+i
    }).join(",")
  var strideArg = indices.map(function(i) {
      return "c"+i
    }).join(",")
  code.push(
    "function "+className+"(a," + shapeArg + "," + strideArg + ",d){this.data=a",
      "this.shape=[" + shapeArg + "]",
      "this.stride=[" + strideArg + "]",
      "this.offset=d|0}",
    "var proto="+className+".prototype",
    "proto.dtype='"+dtype+"'",
    "proto.dimension="+dimension)

  //view.size:
  code.push("Object.defineProperty(proto,'size',{get:function "+className+"_size(){\
return "+indices.map(function(i) { return "this.shape["+i+"]" }).join("*"),
"}})")

  //view.order:
  if(dimension === 1) {
    code.push("proto.order=[0]")
  } else {
    code.push("Object.defineProperty(proto,'order',{get:")
    if(dimension < 4) {
      code.push("function "+className+"_order(){")
      if(dimension === 2) {
        code.push("return (Math.abs(this.stride[0])>Math.abs(this.stride[1]))?[1,0]:[0,1]}})")
      } else if(dimension === 3) {
        code.push(
"var s0=Math.abs(this.stride[0]),s1=Math.abs(this.stride[1]),s2=Math.abs(this.stride[2]);\
if(s0>s1){\
if(s1>s2){\
return [2,1,0];\
}else if(s0>s2){\
return [1,2,0];\
}else{\
return [1,0,2];\
}\
}else if(s0>s2){\
return [2,0,1];\
}else if(s2>s1){\
return [0,1,2];\
}else{\
return [0,2,1];\
}}})")
      }
    } else {
      code.push("ORDER})")
    }
  }

  //view.set(i0, ..., v):
  code.push(
"proto.set=function "+className+"_set("+args.join(",")+",v){")
  if(useGetters) {
    code.push("return this.data.set("+index_str+",v)}")
  } else {
    code.push("return this.data["+index_str+"]=v}")
  }

  //view.get(i0, ...):
  code.push("proto.get=function "+className+"_get("+args.join(",")+"){")
  if(useGetters) {
    code.push("return this.data.get("+index_str+")}")
  } else {
    code.push("return this.data["+index_str+"]}")
  }

  //view.index:
  code.push(
    "proto.index=function "+className+"_index(", args.join(), "){return "+index_str+"}")

  //view.hi():
  code.push("proto.hi=function "+className+"_hi("+args.join(",")+"){return new "+className+"(this.data,"+
    indices.map(function(i) {
      return ["(typeof i",i,"!=='number'||i",i,"<0)?this.shape[", i, "]:i", i,"|0"].join("")
    }).join(",")+","+
    indices.map(function(i) {
      return "this.stride["+i + "]"
    }).join(",")+",this.offset)}")

  //view.lo():
  var a_vars = indices.map(function(i) { return "a"+i+"=this.shape["+i+"]" })
  var c_vars = indices.map(function(i) { return "c"+i+"=this.stride["+i+"]" })
  code.push("proto.lo=function "+className+"_lo("+args.join(",")+"){var b=this.offset,d=0,"+a_vars.join(",")+","+c_vars.join(","))
  for(var i=0; i<dimension; ++i) {
    code.push(
"if(typeof i"+i+"==='number'&&i"+i+">=0){\
d=i"+i+"|0;\
b+=c"+i+"*d;\
a"+i+"-=d}")
  }
  code.push("return new "+className+"(this.data,"+
    indices.map(function(i) {
      return "a"+i
    }).join(",")+","+
    indices.map(function(i) {
      return "c"+i
    }).join(",")+",b)}")

  //view.step():
  code.push("proto.step=function "+className+"_step("+args.join(",")+"){var "+
    indices.map(function(i) {
      return "a"+i+"=this.shape["+i+"]"
    }).join(",")+","+
    indices.map(function(i) {
      return "b"+i+"=this.stride["+i+"]"
    }).join(",")+",c=this.offset,d=0,ceil=Math.ceil")
  for(var i=0; i<dimension; ++i) {
    code.push(
"if(typeof i"+i+"==='number'){\
d=i"+i+"|0;\
if(d<0){\
c+=b"+i+"*(a"+i+"-1);\
a"+i+"=ceil(-a"+i+"/d)\
}else{\
a"+i+"=ceil(a"+i+"/d)\
}\
b"+i+"*=d\
}")
  }
  code.push("return new "+className+"(this.data,"+
    indices.map(function(i) {
      return "a" + i
    }).join(",")+","+
    indices.map(function(i) {
      return "b" + i
    }).join(",")+",c)}")

  //view.transpose():
  var tShape = new Array(dimension)
  var tStride = new Array(dimension)
  for(var i=0; i<dimension; ++i) {
    tShape[i] = "a[i"+i+"]"
    tStride[i] = "b[i"+i+"]"
  }
  code.push("proto.transpose=function "+className+"_transpose("+args+"){"+
    args.map(function(n,idx) { return n + "=(" + n + "===undefined?" + idx + ":" + n + "|0)"}).join(";"),
    "var a=this.shape,b=this.stride;return new "+className+"(this.data,"+tShape.join(",")+","+tStride.join(",")+",this.offset)}")

  //view.pick():
  code.push("proto.pick=function "+className+"_pick("+args+"){var a=[],b=[],c=this.offset")
  for(var i=0; i<dimension; ++i) {
    code.push("if(typeof i"+i+"==='number'&&i"+i+">=0){c=(c+this.stride["+i+"]*i"+i+")|0}else{a.push(this.shape["+i+"]);b.push(this.stride["+i+"])}")
  }
  code.push("var ctor=CTOR_LIST[a.length+1];return ctor(this.data,a,b,c)}")

  //Add return statement
  code.push("return function construct_"+className+"(data,shape,stride,offset){return new "+className+"(data,"+
    indices.map(function(i) {
      return "shape["+i+"]"
    }).join(",")+","+
    indices.map(function(i) {
      return "stride["+i+"]"
    }).join(",")+",offset)}")

  //Compile procedure
  var procedure = new Function("CTOR_LIST", "ORDER", code.join("\n"))
  return procedure(CACHED_CONSTRUCTORS[dtype], order)
}

function arrayDType(data) {
  if(isBuffer(data)) {
    return "buffer"
  }
  if(hasTypedArrays) {
    switch(Object.prototype.toString.call(data)) {
      case "[object Float64Array]":
        return "float64"
      case "[object Float32Array]":
        return "float32"
      case "[object Int8Array]":
        return "int8"
      case "[object Int16Array]":
        return "int16"
      case "[object Int32Array]":
        return "int32"
      case "[object Uint8Array]":
        return "uint8"
      case "[object Uint16Array]":
        return "uint16"
      case "[object Uint32Array]":
        return "uint32"
      case "[object Uint8ClampedArray]":
        return "uint8_clamped"
    }
  }
  if(Array.isArray(data)) {
    return "array"
  }
  return "generic"
}

var CACHED_CONSTRUCTORS = {
  "float32":[],
  "float64":[],
  "int8":[],
  "int16":[],
  "int32":[],
  "uint8":[],
  "uint16":[],
  "uint32":[],
  "array":[],
  "uint8_clamped":[],
  "buffer":[],
  "generic":[]
}

;(function() {
  for(var id in CACHED_CONSTRUCTORS) {
    CACHED_CONSTRUCTORS[id].push(compileConstructor(id, -1))
  }
});

function wrappedNDArrayCtor(data, shape, stride, offset) {
  if(data === undefined) {
    var ctor = CACHED_CONSTRUCTORS.array[0]
    return ctor([])
  } else if(typeof data === "number") {
    data = [data]
  }
  if(shape === undefined) {
    shape = [ data.length ]
  }
  var d = shape.length
  if(stride === undefined) {
    stride = new Array(d)
    for(var i=d-1, sz=1; i>=0; --i) {
      stride[i] = sz
      sz *= shape[i]
    }
  }
  if(offset === undefined) {
    offset = 0
    for(var i=0; i<d; ++i) {
      if(stride[i] < 0) {
        offset -= (shape[i]-1)*stride[i]
      }
    }
  }
  var dtype = arrayDType(data)
  var ctor_list = CACHED_CONSTRUCTORS[dtype]
  while(ctor_list.length <= d+1) {
    ctor_list.push(compileConstructor(dtype, ctor_list.length-1))
  }
  var ctor = ctor_list[d+1]
  return ctor(data, shape, stride, offset)
}

module.exports = wrappedNDArrayCtor

},{"iota-array":9,"is-buffer":10}],13:[function(require,module,exports){
"use strict";

var zeros = require('zeros'),
    moore = require('moore'),
    sphereRandom = require('./sphere-random');

/**
 * Get the squared euclidean distance from two points of arbitrary, but equal, dimensions
 * @param {Array} point1
 * @param {Array} point2
 * @returns {number} Squared euclidean distance
 */
var squaredEuclideanDistance = function squaredEuclideanDistance (point1, point2) {
    var result = 0,
        i = 0;

    for (; i < point1.length; i++) {
        result += Math.pow(point1[i] - point2[i], 2);
    }

    return result;
};

/**
 * Get the neighbourhood ordered by distance, including the origin point
 * @param {int} dimensionNumber Number of dimensions
 * @returns {Array} Neighbourhood
 */
var getNeighbourhood = function getNeighbourhood (dimensionNumber) {
    var neighbourhood = moore(2, dimensionNumber),
        origin = [],
        dimension;

    for (dimension = 0; dimension < dimensionNumber; dimension++) {
        origin.push(0);
    }

    neighbourhood.push(origin);

    // sort by ascending distance to optimize proximity checks
    // see point 5.1 in Parallel Poisson Disk Sampling by Li-Yi Wei, 2008
    // http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.460.3061&rank=1
    neighbourhood.sort(function (n1, n2) {
        var squareDist1 = 0,
            squareDist2 = 0;

        for (var dimension = 0; dimension < dimensionNumber; dimension++) {
            squareDist1 += Math.pow(n1[dimension], 2);
            squareDist2 += Math.pow(n2[dimension], 2);
        }

        if (squareDist1 < squareDist2) {
            return -1;
        } else if(squareDist1 > squareDist2) {
            return 1;
        } else {
            return 0;
        }
    });

    return neighbourhood;
};


/**
 * PoissonDiskSampling constructor
 * @param {Array} shape Shape of the space
 * @param {float} minDistance Minimum distance between each points
 * @param {float} [maxDistance] Maximum distance between each points, defaults to minDistance * 2
 * @param {int} [maxTries] Number of times the algorithm has to try to place a point in the neighbourhood of another points, defaults to 30
 * @param {function|null} [rng] RNG function, defaults to Math.random
 * @constructor
 */
var PoissonDiskSampling = function PoissonDiskSampling (shape, minDistance, maxDistance, maxTries, rng) {
    maxDistance = maxDistance || minDistance * 2;

    this.shape = shape;
    this.dimension = this.shape.length;
    this.minDistance = minDistance;
    this.squaredMinDistance = minDistance * minDistance;
    this.deltaDistance = maxDistance - minDistance;
    this.cellSize = minDistance / Math.sqrt(this.dimension);
    this.maxTries = maxTries || 30;
    this.rng = rng || Math.random;

    this.neighbourhood = getNeighbourhood(this.dimension);

    this.currentPoint = null;
    this.processList = [];
    this.samplePoints = [];

    // cache grid

    this.gridShape = [];

    for (var i = 0; i < this.dimension; i++) {
        this.gridShape.push(Math.ceil(shape[i] / this.cellSize));
    }

    this.grid = zeros(this.gridShape, 'uint32'); //will store references to samplePoints
};

PoissonDiskSampling.prototype.shape = null;
PoissonDiskSampling.prototype.dimension = null;
PoissonDiskSampling.prototype.minDistance = null;
PoissonDiskSampling.prototype.squaredMinDistance = null;
PoissonDiskSampling.prototype.deltaDistance = null;
PoissonDiskSampling.prototype.cellSize = null;
PoissonDiskSampling.prototype.maxTries = null;
PoissonDiskSampling.prototype.rng = null;
PoissonDiskSampling.prototype.neighbourhood = null;

PoissonDiskSampling.prototype.currentPoint = null;
PoissonDiskSampling.prototype.processList = null;
PoissonDiskSampling.prototype.samplePoints = null;
PoissonDiskSampling.prototype.gridShape = null;
PoissonDiskSampling.prototype.grid = null;

/**
 * Add a totally random point in the grid
 * @returns {Array} The point added to the grid
 */
PoissonDiskSampling.prototype.addRandomPoint = function () {
    var point = new Array(this.dimension);

    for (var i = 0; i < this.dimension; i++) {
        point[i] = this.rng() * this.shape[i];
    }

    return this.directAddPoint(point);
};

/**
 * Add a given point to the grid
 * @param {Array} point Point
 * @returns {Array|null} The point added to the grid, null if the point is out of the bound or not of the correct dimension
 */
PoissonDiskSampling.prototype.addPoint = function (point) {
    var dimension,
        valid = true;

    if (point.length === this.dimension) {
        for (dimension = 0; dimension < this.dimension && valid; dimension++) {
            valid = (point[dimension] >= 0 && point[dimension] <= this.shape[dimension]);
        }
    } else {
        valid = false;
    }

    return valid ? this.directAddPoint(point) : null;
};

/**
 * Add a given point to the grid, without any check
 * @param {Array} point Point
 * @returns {Array} The point added to the grid
 * @protected
 */
PoissonDiskSampling.prototype.directAddPoint = function (point) {
    var internalArrayIndex = 0,
        stride = this.grid.stride,
        dimension;

    this.processList.push(point);
    this.samplePoints.push(point);

    for (dimension = 0; dimension < this.dimension; dimension++) {
        internalArrayIndex += ((point[dimension] / this.cellSize) | 0) * stride[dimension];
    }

    this.grid.data[internalArrayIndex] = this.samplePoints.length; // store the point reference

    return point;
};

/**
 * Check whether a given point is in the neighbourhood of existing points
 * @param {Array} point Point
 * @returns {boolean} Whether the point is in the neighbourhood of another point
 * @protected
 */
PoissonDiskSampling.prototype.inNeighbourhood = function (point) {
    var dimensionNumber = this.dimension,
        stride = this.grid.stride,
        neighbourIndex,
        internalArrayIndex,
        dimension,
        currentDimensionValue,
        existingPoint;

    for (neighbourIndex = 0; neighbourIndex < this.neighbourhood.length; neighbourIndex++) {
        internalArrayIndex = 0;

        for (dimension = 0; dimension < dimensionNumber; dimension++) {
            currentDimensionValue = ((point[dimension] / this.cellSize) | 0) + this.neighbourhood[neighbourIndex][dimension];

            if (currentDimensionValue >= 0 && currentDimensionValue < this.gridShape[dimension]) {
                internalArrayIndex += currentDimensionValue * stride[dimension];
            }
        }

        if (this.grid.data[internalArrayIndex] !== 0) {
            existingPoint = this.samplePoints[this.grid.data[internalArrayIndex] - 1];

            if (squaredEuclideanDistance(point, existingPoint) < this.squaredMinDistance) {
                return true;
            }
        }
    }

    return false;
};

/**
 * Try to generate a new point in the grid, returns null if it wasn't possible
 * @returns {Array|null} The added point or null
 */
PoissonDiskSampling.prototype.next = function () {
    var tries,
        angle,
        distance,
        currentPoint,
        newPoint,
        inShape,
        i;

    while (this.processList.length > 0) {
        if (this.currentPoint === null) {
            this.currentPoint = this.processList.shift();
        }

        currentPoint = this.currentPoint;

        for (tries = 0; tries < this.maxTries; tries++) {
            inShape = true;
            distance = this.minDistance + this.deltaDistance * this.rng();

            if (this.dimension === 2) {
                angle = this.rng() * Math.PI * 2;
                newPoint = [
                    Math.cos(angle),
                    Math.sin(angle)
                ];
            } else {
                newPoint = sphereRandom(this.dimension, this.rng);
            }

            for (i = 0; inShape && i < this.dimension; i++) {
                newPoint[i] = currentPoint[i] + newPoint[i] * distance;
                inShape = (newPoint[i] >= 0 && newPoint[i] <= this.shape[i] - 1)
            }

            if (inShape && !this.inNeighbourhood(newPoint)) {
                return this.directAddPoint(newPoint);
            }
        }

        if (tries === this.maxTries) {
            this.currentPoint = null;
        }
    }

    return null;
};

/**
 * Automatically fill the grid, adding a random point to start the process if needed.
 * Will block the thread, probably best to use it in a web worker or child process.
 * @returns {Array[]} Sample points
 */
PoissonDiskSampling.prototype.fill = function () {
    if (this.samplePoints.length === 0) {
        this.addRandomPoint();
    }

    while(this.next()) {}

    return this.samplePoints;
};

/**
 * Get all the points in the grid.
 * @returns {Array[]} Sample points
 */
PoissonDiskSampling.prototype.getAllPoints = function () {
    return this.samplePoints;
};

/**
 * Reinitialize the grid as well as the internal state
 */
PoissonDiskSampling.prototype.reset = function () {
    var gridData = this.grid.data,
        i = 0;

    // reset the cache grid
    for (i = 0; i < gridData.length; i++) {
        gridData[i] = 0;
    }

    // new array for the samplePoints as it is passed by reference to the outside
    this.samplePoints = [];

    // reset the internal state
    this.currentPoint = null;
    this.processList.length = 0;
};

module.exports = PoissonDiskSampling;

},{"./sphere-random":14,"moore":11,"zeros":16}],14:[function(require,module,exports){
"use strict";

// sphere-random module by Mikola Lysenko under the MIT License
// waiting for https://github.com/scijs/sphere-random/pull/1 to be merged

module.exports = sampleSphere;

/**
 * @param {int} d Dimensions
 * @param {Function} rng
 * @returns {Array}
 */
function sampleSphere(d, rng) {
    var v = new Array(d),
        d2 = Math.floor(d/2) << 1,
        r2 = 0.0,
        rr,
        r,
        theta,
        h,
        i;

    for (i = 0; i < d2; i += 2) {
        rr = -2.0 * Math.log(rng());
        r =  Math.sqrt(rr);
        theta = 2.0 * Math.PI * rng();

        r2+= rr;
        v[i] = r * Math.cos(theta);
        v[i+1] = r * Math.sin(theta);
    }

    if (d % 2) {
        var x = Math.sqrt(-2.0 * Math.log(rng())) * Math.cos(2.0 * Math.PI * rng());
        v[d - 1] = x;
        r2+= Math.pow(x, 2);
    }

    h = 1.0 / Math.sqrt(r2);

    for (i = 0; i < d; ++i) {
        v[i] *= h;
    }

    return v;
}

},{}],15:[function(require,module,exports){
var bundleFn = arguments[3];
var sources = arguments[4];
var cache = arguments[5];

var stringify = JSON.stringify;

module.exports = function (fn, options) {
    var wkey;
    var cacheKeys = Object.keys(cache);

    for (var i = 0, l = cacheKeys.length; i < l; i++) {
        var key = cacheKeys[i];
        var exp = cache[key].exports;
        // Using babel as a transpiler to use esmodule, the export will always
        // be an object with the default export as a property of it. To ensure
        // the existing api and babel esmodule exports are both supported we
        // check for both
        if (exp === fn || exp && exp.default === fn) {
            wkey = key;
            break;
        }
    }

    if (!wkey) {
        wkey = Math.floor(Math.pow(16, 8) * Math.random()).toString(16);
        var wcache = {};
        for (var i = 0, l = cacheKeys.length; i < l; i++) {
            var key = cacheKeys[i];
            wcache[key] = key;
        }
        sources[wkey] = [
            'function(require,module,exports){' + fn + '(self); }',
            wcache
        ];
    }
    var skey = Math.floor(Math.pow(16, 8) * Math.random()).toString(16);

    var scache = {}; scache[wkey] = wkey;
    sources[skey] = [
        'function(require,module,exports){' +
            // try to call default if defined to also support babel esmodule exports
            'var f = require(' + stringify(wkey) + ');' +
            '(f.default ? f.default : f)(self);' +
        '}',
        scache
    ];

    var workerSources = {};
    resolveSources(skey);

    function resolveSources(key) {
        workerSources[key] = true;

        for (var depPath in sources[key][1]) {
            var depKey = sources[key][1][depPath];
            if (!workerSources[depKey]) {
                resolveSources(depKey);
            }
        }
    }

    var src = '(' + bundleFn + ')({'
        + Object.keys(workerSources).map(function (key) {
            return stringify(key) + ':['
                + sources[key][0]
                + ',' + stringify(sources[key][1]) + ']'
            ;
        }).join(',')
        + '},{},[' + stringify(skey) + '])'
    ;

    var URL = window.URL || window.webkitURL || window.mozURL || window.msURL;

    var blob = new Blob([src], { type: 'text/javascript' });
    if (options && options.bare) { return blob; }
    var workerUrl = URL.createObjectURL(blob);
    var worker = new Worker(workerUrl);
    worker.objectURL = workerUrl;
    return worker;
};

},{}],16:[function(require,module,exports){
"use strict"

var ndarray = require("ndarray")

function dtypeToType(dtype) {
  switch(dtype) {
    case 'uint8':
      return Uint8Array;
    case 'uint16':
      return Uint16Array;
    case 'uint32':
      return Uint32Array;
    case 'int8':
      return Int8Array;
    case 'int16':
      return Int16Array;
    case 'int32':
      return Int32Array;
    case 'float':
    case 'float32':
      return Float32Array;
    case 'double':
    case 'float64':
      return Float64Array;
    case 'uint8_clamped':
      return Uint8ClampedArray;
    case 'generic':
    case 'buffer':
    case 'data':
    case 'dataview':
      return ArrayBuffer;
    case 'array':
      return Array;
  }
}

module.exports = function zeros(shape, dtype) {
  dtype = dtype || 'float64';
  var sz = 1;
  for(var i=0; i<shape.length; ++i) {
    sz *= shape[i];
  }
  return ndarray(new (dtypeToType(dtype))(sz), shape);
}

},{"ndarray":12}]},{},[3]);
