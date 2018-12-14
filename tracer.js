const Geometry = require('./geometry.js');
const work = require('webworkify');
let colorConvert = require('color-convert');

/*
    ATTENTION! Hey Alan, because I know you'll forget, here's what you put in the command line 
    to get browserify to auto compile:
        watchify tracer.js -o bundle.js -v
    Cheers! (You're welcome)

    And for you other people who don't remember, cd to the ray_tracing folder and type 
        python3 -m http.server
    and go to localhost [whatever it says]
*/

let featureCollection = null;
let lightSources = [];
let environment = {
    ambientConstant:.4,
    diffuseConstant:.5,
    specularStrengthConstant:2,
    specularNarrownessConstant:80,
    maxRecursion:5,
    samplesPerPoint:50,
    medium:{
        indexOfRefraction:1,
        opacity:.000,
        color:[54, 110, 169]
    },
    rayDepth:5,
    rayBranching:1
}

const brightnessTransform = x => x**(1/2);

// body onload
window.init = function(){
    initFeatures();
    trace(image => {
        // do nothing
    }, false, 'canvas');
}

function initFeatures(){

    let features = [];

    let k = 6;
    let planes = [
        [[k,k,k],[k,-k,k],[k,-k,-k],[k,k,-k]],
        [[-k,k,k],[-k,-k,k],[-k,-k,-k],[-k,k,-k]],
        [[k,k,k],[-k,k,k],[-k,k,-k],[k,k,-k]],
        [[k,-k,k],[-k,-k,k],[-k,-k,-k],[k,-k,-k]],
        [[k,k,k],[-k,k,k],[-k,-k,k],[k,-k,k]],
        [[k,k,-k],[-k,k,-k],[-k,-k,-k],[k,-k,-k]]
    ];
    let normals = [[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]];
    for(let p=0; p<planes.length; p++){
        let plane = Geometry.PolyPlane.uniform(planes[p],normals[p],[255,255,255],p === 0 ? true : false);
        plane.properties.brdf = Geometry.BRDF.variableGlossy(p === 3 ? .4 : .8);
        features.push(plane);
    }

    // cube on the ground 
    function transform(x){
        return Geometry.plus(Geometry.scale(.05,x), [-1,-k+k/20,-2]);
    }
    let transformedNormals = normals.map(x => Geometry.scale(-1,x));
    for(let p = 0; p < planes.length; p++){
        let transformedPlane = planes[p].map(x => transform(x));
        let newPlane = Geometry.PolyPlane.uniform(transformedPlane, transformedNormals[p], [249, 229, 118], false);
        newPlane.properties.brdf = Geometry.BRDF.variableGlossy(1);
        features.push(newPlane);
    }

    let r = 1.5;
    let pos = [-2,-k+r,0];

    // ball on the ground
    let ball = Geometry.Sphere.uniform(pos, r, [100,255,255], false);
    ball.properties.brdf = Geometry.BRDF.variableGlossy(0);

    features.push(ball);

    featureCollection = new Geometry.FeatureCollection(features);
}

function trace(traceFinished, basic, imageID){
    const res = 600;
    let camera = new Geometry.Camera([-3,-5,-5.5], [[1,0,0],[0,1,0],[0,0,3]], 5, res);
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
    for(let i = 0; i < workerJobs.length; i++){
        workerJobs[i] = shuffle(workerJobs[i]);
    }

    // let the workers run wild
    let workers = []

    let onfinish = () => {};

    let finishedWorkers = 0;
    for(let k = 0; k < workerJobs.length; k++){
        let newWorker = work(require('./rayWorker.js'));
        newWorker.postMessage(['overhead', camera, featureCollection, lightSources, environment, basic]);
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