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
    samplesPerPoint:200,
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

    // floor
    const floorHeight = .1;
    const floor = Geometry.PolyPlane.uniform(
        [[-5, floorHeight, -1], [-5, floorHeight, 3], [5, floorHeight, 3], [5, floorHeight, -1]], // points
        [0,-1,0], // normal
        [255,255,255], // color
        false // is light source
    );
    floor.properties.brdf = Geometry.BRDF.variableGlossy(1);

    // back wall 
    const wallWidth = 1.5;
    const wallHeight = 1;
    const wall1 = Geometry.PolyPlane.uniform(
        [[-wallWidth/2, floorHeight, 2], [wallWidth/2, floorHeight, 2], [wallWidth/2, floorHeight-wallHeight, 2], [-wallWidth/2, floorHeight-wallHeight]],
        [0,0,-1], // normal
        [255,255,255], //color 
        false // is light source
    );
    wall1.properties.brdf = Geometry.BRDF.variableGlossy(1);

    // light source 
    const p1 = [-.25, -.3, .9];
    const p2 = [-.25, -.3, 1.1];
    const p3 = [-.3,-.1,1];
    const normal = Geometry.normalize(Geometry.cross(Geometry.minus(p2, p1), Geometry.minus(p3, p1)));
    let light = Geometry.PolyPlane.uniform(
        [p1, p2, p3], // points
        normal, // normal
        [255,255,255], // color
        true // is light source
    );

    light = Geometry.Sphere.uniform(p1, .1, [255,255,255], true);
    light.properties.brdf = Geometry.BRDF.variableGlossy(1);
    light.properties.intensity = 2000;

    // ball on floor 
    const ball1R = .07;
    const ball1 = Geometry.Sphere.uniform([0,floorHeight-ball1R, 1], ball1R, [255,20,20], false);
    ball1.properties.brdf = Geometry.BRDF.variableGlossy(1);

    features.push(floor, light, ball1, wall1);
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
    camera = new Geometry.Camera([0,0,0], [[1,0,0],[0,1,0],[0,0,1]], 35/1000, 37/1000, 30/1000, .01/1000, 600);
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