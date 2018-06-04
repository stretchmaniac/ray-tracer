
const Geometry = require('./geometry.js');

let featureCollection = null;
let lightSources = [];
let environment = {
    ambientConstant:.4,
    diffuseConstant:.5,
    specularStrengthConstant:2,
    specularNarrownessConstant:80,
    maxRecursion:5,
    medium:{
        indexOfRefraction:1,
        opacity:.000,
        color:[54, 110, 169]
    }
}

// body onload
window.init = function(){
    initFeatures();
    let image = trace();
    show(image);
}

function initFeatures(){

    let features = [];

    // ball sitting on a large plane
    sphere = Geometry.Sphere.uniform([0, 3, 10], 3, .1, 0, [0, 30, 0], 1.5);
    sphere2 = Geometry.Sphere.uniform([2.5,1,7], 1, 1, 1, [235, 0, 150], 1);
    let w = 30;
    let plane = Geometry.PolyPlane.uniform(
        [[w,0,-2*w], [w,0,2*w], [-w,0,2*w], [-w,0,-2*w]],
        [0,1,0],
        1, .1, [1,1,1], 1
    );
    // make the plane checkered
    plane.properties.colorFunction = function(loc){

        function mod(n, m) {
            return ((n % m) + m) % m;
        }

        // checker width
        let width = 4;
        let color1 = [150,150,150];
        let color2 = [230,230,230];
        let [x,y,z] = loc;
        let value = mod(Math.round(x/2) + Math.round(z/2), 2);
        return value === 0 ? color1 : color2;
    };

    if(Geometry.dot(plane.normal, [0,1,0]) < 0){
        plane.normal = Geometry.scale(-1, plane.normal);
    }

    features.push(sphere);
    features.push(sphere2);
    features.push(plane);

    // light sources consist of a feature used for specular reflection
    // and a sampling of points on the feature used for diffuse lighting
    // for now we'll use spheres with a point source at the middle
    let sources = [[5,25,10], [15,25,15]]
    let radius = 3;
    for(let center of sources){
        let newFeature = new Geometry.Sphere(center, radius, 
            {
                lightSource: {
                    surfaceIntensity: .8,
                    totalIntensity: 500,
                    color:[255, 255, 255],
                    diffuseSample: [center]
                }
            }
        );
        features.push(newFeature);
        lightSources.push(newFeature);
    }

    featureCollection = new Geometry.FeatureCollection(features);
}

function trace(){
    let camera = new Geometry.Camera([0,5,-5], [[1,0,0],[0,1,0],[0,0,5]], 5);
    const res = 800;
    let rays = camera.getRays(res);
    let image = fillArray(res,res,0);
    let childRays = fillArray(res,res,0);
    
    for(let ray of rays){
        ray.medium = environment.medium;
        for(let child of ray.children){
            child.medium = ray.medium;
        }
        let resultRay = traceRay(ray, environment.maxRecursion);
        image[ray.pixel[0]][ray.pixel[1]] = resultRay.color;
        childRays[ray.pixel[0]][ray.pixel[1]] = ray.children;
    }

    // now for some edge cleanup
    const edgeColorDifferential = 10;

    colorDiff = (a,b) => Geometry.mag(Geometry.minus(a,b));
    let newImage = fillArray(res, res, 0);

    for(let i = 0; i < res; i++){
        for(let j = 0; j < res; j++){
            newImage[i][j] = image[i][j];
            let castMore = 
                (i > 0 && colorDiff(image[i][j], image[i-1][j]) > edgeColorDifferential) ||
                (j > 0 && colorDiff(image[i][j], image[i][j-1]) > edgeColorDifferential) ||
                (i < res - 1 && colorDiff(image[i][j],image[i+1][j]) > edgeColorDifferential) ||
                (j < res - 1 && colorDiff(image[i][j], image[i][j+1]) > edgeColorDifferential);
            if(castMore){
                auxRays = childRays[i][j];
                resultColors = auxRays.map(x => traceRay(x, environment.maxRecursion).color);
                // take the weighted average of all the colors
                let total1 = Geometry.scale(1/4, Geometry.plus(...resultColors));
                newImage[i][j] = Geometry.plus(Geometry.scale(.5,image[i][j]), Geometry.scale(.5,total1));
            }
        }
    }

    image = newImage;

    return image;
}

function show(image){
    let canvas = document.getElementById('canvas');
    let ctx = canvas.getContext('2d');
    for(let y = 0; y < image.length; y++){
        for(let x = 0; x < image[y].length; x++){
            let color = image[x][image.length - y - 1];
            ctx.fillStyle = 'rgb('+color[0]+','+color[1]+','+color[2]+')';
            ctx.fillRect(x,y,1,1);
        }
    }
}

function traceRay(ray, maxRecursion, currentDepth=0){
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
            let directionToLight = Geometry.normalize(Geometry.minus(lightPoint, collision.pos));
            let newRay = new Geometry.Ray(
                Geometry.plus(collision.pos,Geometry.scale(1e-8,directionToLight)), 
                directionToLight
            );
            newRay.medium = ray.medium;

            let testCollisions = featureCollection.intersect(newRay);
            if(testCollisions.length > 0 && testCollisions[0].object === lightFeature){
                // use this light source sample point
                let distToLight = Geometry.mag(Geometry.minus(lightPoint, collision.pos));
                let mediumOpacityConst = Math.min(1, ray.medium.opacity * distToLight);
                // it is important that directionToLight and collision.surface.normal are both unit vectors
                // diffuse
                let dBrightness = environment.diffuseConstant*Geometry.dot(directionToLight, collision.surface.normal) * diffuseIntensity / distToLight**2;
                ray.brightness += dBrightness;
                // specular
                let reflectDir = Geometry.normalize(Geometry.plus(
                    Geometry.scale(2, 
                        Geometry.minus(
                            Geometry.project(directionToLight, collision.surface.normal),
                            directionToLight
                        )
                    ), 
                    directionToLight
                ));
                let c1 = Math.abs(Geometry.dot(reflectDir, Geometry.normalize(ray.direction)))**environment.specularNarrownessConstant
                ray.brightness += environment.specularStrengthConstant * c1 * diffuseIntensity*collision.surface.specularity / distToLight**2;
            }
        }
    }

    // now do a new ray!
    let outwardDir = Geometry.scale(-1, ray.direction);
    let bounceDirection = Geometry.plus(
        Geometry.scale(2, 
            Geometry.minus(
                Geometry.project(outwardDir, collision.surface.normal),
                outwardDir
            )
        ), 
        outwardDir
    );

    // Specular Reflection Ray
    let childRay = new Geometry.Ray(Geometry.plus(collision.pos, Geometry.scale(1e-8,bounceDirection)), bounceDirection);
    // we say the ray is leaving the object if the normal is in the direction of the ray
    let leaving = Geometry.dot(childRay.direction, collision.surface.normal) > 0;
    if(leaving){
        childRay.medium = environment.medium;
    }else{
        childRay.medium = ray.medium;
    }
    let resultChild = traceRay(childRay, maxRecursion, currentDepth + 1);

    // this ray is on the other side of the surface
    // snell's law! 
    let n1 = ray.medium.indexOfRefraction, n2 = 0;
    let refractionIntoSurface = leaving;
    if(!refractionIntoSurface){
        n2 = environment.medium.indexOfRefraction;
    }else{
        n2 = collision.surface.indexOfRefraction;
    }
    let normalWithRay = Geometry.project(ray.direction, collision.surface.normal);
    let costheta1 = Geometry.dot(ray.direction, normalWithRay) / (Geometry.mag(ray.direction)*Geometry.mag(normalWithRay));
    let sintheta1 = Math.sqrt(1-costheta1**2);
    let sintheta2 = (n1/n2)*sintheta1;
    let totalInternalReflection = false;
    let childRay2 = null;
    if(Math.abs(sintheta2) > 1){
        totalInternalReflection = true;
    }else{
        let theta2 = Math.asin(sintheta2);
        let a = Geometry.mag(normalWithRay);
        let theta1 = Math.acos(costheta1);
        let f = a*(Math.tan(theta1) - Math.tan(theta2));
        let refractionDirection = Geometry.plus(
            ray.direction,
            Geometry.scale(f, Geometry.normalize(Geometry.minus(normalWithRay, ray.direction)))
        );
        childRay2 = new Geometry.Ray(Geometry.plus(collision.pos, Geometry.scale(1e-8,normalWithRay)), refractionDirection);
        childRay2.medium = !refractionIntoSurface ? environment.medium : {indexOfRefraction: n2, opacity: collision.surface.opacity, color: collision.surface.color};
        if(refractionIntoSurface && collision.surface.opacity === 1){
            childRay2.color = collision.surface.color;
        }else{
            childRay2 = traceRay(childRay2, maxRecursion, currentDepth + 1);
        }
    }

    // combine the colors
    let rayDist = Geometry.mag(Geometry.minus(ray.start, collision.pos));
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