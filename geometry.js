let gaussianRand = require('gauss-random');

class Camera{
    constructor(location, orientation, lensWidth, resolution){
        this.location = location;
        this.orientation = orientation;
        this.lensWidth = lensWidth;
        this.resolution = resolution;
    }

    // suppose the detector lies on the interval [-1,1] x [-1,1]
    getRay(unitLocation){
        let rawPt = [unitLocation[0]*this.lensWidth/2, -unitLocation[1]*this.lensWidth/2, 1];
        let rotatedPt = normalize(matColMult(this.orientation, rawPt));
        let newRay = new Ray(plus(this.location, rotatedPt), rotatedPt);

        return newRay;
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
        return feature.sample();
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
    lineSegmentValid(p1, p2){
        let ray = new Ray(p1, minus(p2, p1));
        let intersections = this.intersect(ray);
        if(intersections.length === 0){
            return true;
        }
        return mag(minus(p2,p1)) < mag(minus(intersections[0].pos, p1));
    }
}

// bidirectional reflectance distribution function
class BRDF{
    // n ranges from 0 to 1, where 0 is perfectly reflective and 
    // 1 is perfectly diffuse
    static variableGlossy(n){
        let brdf = {};
        brdf.method = 'quick_n_dirty';
        brdf.properties = {
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
        // reflection/transmission probabilities and 
        // directions given a perfectly transparent (non-frosty) surface
        brdf.rt = function(inVec, normal){
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

            // calculate reflectance
            const Rs = ((n1*Math.cos(theta)-n2*Math.cos(refractionAngle)) / (n1*Math.cos(theta)+n2*Math.cos(refractionAngle)))**2;
            const Rp = ((n1*Math.cos(refractionAngle)-n2*Math.cos(theta)) / (n1*Math.cos(refractionAngle)+n2*Math.cos(theta)))**2;

            const totalReflectance = sPol*Rs + pPol*Rp;
            return {
                r: {dir:normalize(reflect(inVec, normal)), p:totalReflectance},
                t: {dir:normalize(snell(inVec, normal, refractionAngle)), p:1-totalReflectance}
            };
        }
        brdf.__proto__ = BRDF.prototype;
        return brdf;
    }

    // gives a sample of the most likely vectors to 
    // reflect in the direction of outVec
    // should return a list of rays
    inverseReflectionSample(normal, outVec, numToSample){
        if(this.method === 'quick_n_dirty' || this.method === 'transparent'){
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
            if(Math.random() < probs.r.p){
                // reflection 
                return BRDF.gaussDir(probs.r.dir, this.properties.lobeSize, normal);
            }else{
                // transmission 
                return BRDF.gaussDir(probs.t.dir, this.properties.lobeSize, scale(-1, normal));
            }
        }
    }

    static gaussDir(normDir, lobeSize, normal){
        let validDirection = false;
        let finalPt = null;
        while(validDirection === false){
            let randomDirection = [Math.random()-.5, Math.random()-.5, Math.random()-.5];
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
        }
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
        let projection = normalize(pt);
        let direction = projection;
        // map onto surface of sphere 
        projection = plus(this.center, scale(this.radius, projection));
        return {pt: projection, dir: direction, color:this.properties.color};
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
        return {pt: realRandomPt, dir: this.normal, color:this.properties.color};
    }

    deserialize(){
        if(this.properties.brdf){
            this.properties.brdf.__proto__ = BRDF.prototype;
        }
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
function traceRay(ray, featureCollection, lightCollection, environment, maxDepth){
    // 1. construct bounce path
    // 2. construct random light path
    // 3. connect ends for bidirectional path tracing
    let latestRay = ray;
    let path = [];
    while(path.length < maxDepth){
        // trace latestRay to new intersection
        // move the ray the slightest bit forward
        latestRay.start = plus(latestRay.start, scale(1e-8, latestRay.direction));
        // find first collision 
        let collisions = featureCollection.intersect(latestRay);
        if(collisions.length === 0){
            // this ray has reached the inky blackness of the utter void
            path.push({color:[0,0,0], rayColor: [0,0,0], pos:plus(latestRay.start, latestRay.direction)});
            break;
        }
        let collision = collisions[0];

        let collisionObj = collision.object;
        let collisionPos = collision.pos;
        let normal = collision.surface.normal;
        let color = collision.surface.color;

        // add to path 
        let node = {
            pos: collisionPos,
            object: collisionObj,
            normal: normal,
            color: color
        };
        path.push(node);

        if(node.object.properties.lightSource){
            node.rayColor = node.color;
            break;
        }

        if(path.length === maxDepth){
            node.rayColor = [0,0,0];
        }

        // now prepare to continue (i.e. reflection)
        let newDir = collisionObj.properties.brdf.sample(latestRay.direction, normal);
        latestRay = new Ray(collisionPos, newDir);
    }

    let lightSample = lightCollection.sample();

    latestRay = new Ray(lightSample.pt, lightSample.dir);
    lightPath = [
        {pos:lightSample.pt, rayColor:lightSample.color}
    ];
    while(lightPath.length < maxDepth){
        // same kinda deal as before
        latestRay.start = plus(latestRay.start, scale(1e-8, latestRay.direction));
        // find first collision 
        let collisions = featureCollection.intersect(latestRay);
        if(collisions.length === 0){
            // this ray has reached the inky blackness of the utter void
            break;
        }
        let collision = collisions[0];
        let collisionObj = collision.object;
        let collisionPos = collision.pos;
        let normal = collision.surface.normal;
        let color = collision.surface.color;

        // add to path 
        lightPath.push({
            pos: collisionPos,
            object: collisionObj,
            normal: normal,
            color: color
        });

        // now prepare to continue (i.e. reflection)
        let newDir = collisionObj.properties.brdf.sample(latestRay.direction, normal);
        latestRay = new Ray(collisionPos, newDir);
    }

    // now connect the dots
    // we're looking for contiguous paths from a light source 
    // to the camera

    // calculate colors for lightRay, starting from the light source
    for(let k = 1; k < lightPath.length; k++){
        let prevNode = lightPath[k-1];
        let currentNode = lightPath[k];

        // combine color of surface with color from prevNode
        let newColor = [0,0,0];
        for(let j = 0; j < 3; j++){
            newColor[j] = prevNode.rayColor[j] * currentNode.color[j] / 255;
        }
        currentNode.rayColor = newColor;
    }

    // now follow the normal ray up
    for(let k = 1; k < path.length; k++){
        let prevNode = path[path.length - k];
        let currentNode = path[path.length - k - 1];
        let outDir = undefined;
        if(k < path.length - 1){
            let nextNode = path[path.length - k - 2];
            outDir = normalize(minus(nextNode.pos, currentNode.pos));
        }else{
            let nextPos = ray.start;
            outDir = normalize(minus(nextPos, currentNode.pos));
        }

        // integrate over all light path options (+ prevNode)
        let inDir = normalize(minus(currentNode.pos, prevNode.pos));
        let = currentNode.object.properties.brdf.reflectFunc(inDir, currentNode.normal, outDir);
        // apply absorption
        let absorb = c => {
            let nC = []
            for(let i = 0; i < 3; i++){
                nC.push(c[i] * currentNode.color[i] / 255);
            }
            return nC;
        }
        let newColor = scale(totalWeight, absorb(prevNode.rayColor));
        for(let j = 0; j < lightPath.length; j++){
            let lightNode = lightPath[j];
            // test whether this path is uninterrupted
            if(featureCollection.lineSegmentValid(currentNode.pos, lightNode.pos)){
                inDir = normalize(minus(currentNode.pos, lightNode.pos));
                let weight = currentNode.object.properties.brdf.reflectFunc(inDir, currentNode.normal, outDir);
                newColor = plus(newColor, scale(weight, absorb(lightNode.rayColor)));
                totalWeight += weight;
            }
        }
        let finalColor = scale(1/totalWeight, newColor);
        currentNode.rayColor = finalColor;
    }

    return path[0].rayColor;
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
