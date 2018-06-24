
const Geometry = require('./geometry.js');
const work = require('webworkify');

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
    trace(image => {
        show(image);
    });
}

function initFeatures(){

    let features = [];

    let y = -3;
    x = 20;
    z = 20;
    let plane1 = Geometry.PolyPlane.uniform([[x,y,z],[0,y,z],[0,y,0],[x,y,0]],[0,1,0],1,0,[100,100,100],1);
    let plane2 = Geometry.PolyPlane.uniform([[-x,y,z],[0,y,z],[0,y,0],[-x,y,0]],[0,1,0],1,0,[255,255,255],1);
    let plane3 = Geometry.PolyPlane.uniform([[x,y,-z],[0,y,-z],[0,y,0],[x,y,0]],[0,1,0],1,0,[255,255,255],1);
    let plane4 = Geometry.PolyPlane.uniform([[-x,y,-z],[0,y,-z],[0,y,0],[-x,y,0]],[0,1,0],1,0,[100,100,100],1);

    features.push(plane1,plane2,plane3,plane4);

    // ball on the ground
    let ball = Geometry.Sphere.uniform([2,0,5], 3, 1, 1, [100,255,0], 1);

    features.push(ball);

    // light sources consist of a feature used for specular reflection
    // and a sampling of points on the feature used for diffuse lighting
    // for now we'll use spheres with a point source at the middle
    let sources = [[4,15,0]]
    let radius = 1;
    for(let center of sources){
        let newFeature = new Geometry.Sphere(center, radius, 
            {
                lightSource: {
                    surfaceIntensity: 1,
                    totalIntensity: 100,
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

function trace(traceFinished){
    let camera = new Geometry.Camera([0,0,-5], [[1,0,0],[0,1,0],[0,0,3]], 5);
    const res = 500;
    let rays = camera.getRays(res);
    let image = fillArray(res,res,0);
    let childRays = fillArray(res,res,0);

    let workerJobs = []
    // 8 logical cores = 7 workers + main thread
    const threads = 6;
    for(let k = 0; k < threads; k++){
        workerJobs.push([]);
    }
    
    let workerIndex = 0;
    for(let ray of rays){
        ray.medium = environment.medium;
        for(let child of ray.children){
            child.medium = ray.medium;
        }

        workerJobs[workerIndex].push(ray);
        workerIndex = (workerIndex + 1) % threads;
    }

    // let the workers run wild
    let workers = []
    let finishedRays = []

    let onfinish = () => {
        for(let ray of finishedRays){
            image[ray.pixel[0]][ray.pixel[1]] = ray.color;
            childRays[ray.pixel[0]][ray.pixel[1]] = ray.children;
        }

        let finalImage = edgeCleanup(image, res, childRays);
        traceFinished(finalImage);
    }

    let finishedWorkers = 0;
    for(let k = 0; k < workerJobs.length; k++){
        let newWorker = work(require('./rayWorker.js'));
        newWorker.postMessage(['overhead', featureCollection, lightSources, environment]);
        newWorker.postMessage(['job', workerJobs[k]]);
        (function(){
            newWorker.addEventListener('message', function(e){
                for(let r of e.data[0]){
                    // deserialize
                    r.__proto__ = Geometry.Ray.prototype;
                    r.deserialize();
                    finishedRays.push(r);
                }
                finishedWorkers++;
                console.log('finished workers: '+finishedWorkers)
                if(finishedWorkers === threads){
                    // we're done
                    onfinish();
                }
            });
        })();
    }
}

function edgeCleanup(image, res, childRays){
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
                resultColors = auxRays.map(x => Geometry.traceRay(featureCollection,environment,lightSources, x, environment.maxRecursion).color);
                // take the weighted average of all the colors
                let total1 = Geometry.scale(1/4, Geometry.plus(...resultColors));
                newImage[i][j] = Geometry.plus(Geometry.scale(.5,image[i][j]), Geometry.scale(.5,total1));
            }
        }
    }

    return newImage;
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