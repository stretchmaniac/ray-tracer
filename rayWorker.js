
const Geometry = require('./geometry.js');
const PoissonDisk = require('poisson-disk-sampling')

module.exports = function(self){

    let featureCollection = null;
    let lightSources = [];
    let environment = {};
    let basic = false;
    let camera = null;

    self.addEventListener('message',function(e){
        let type = e.data[0];
        if(type === 'overhead'){
            camera = e.data[1];
            featureCollection = e.data[2];
            lightSources = e.data[3];
            environment = e.data[4];
            basic = e.data[5];
            // reattach methods that were stripped during serialization
            featureCollection.__proto__ = Geometry.FeatureCollection.prototype;
            featureCollection.deserialize();
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
                let averageColor = [0,0,0];

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
                    let ray = camera.getRay(pos);
                    let localAverage = [0,0,0];
                    for(let j = 0; j < environment.samplesPerPoint; j++){
                        const result = Geometry.traceRay(featureCollection, lightSources, environment, ray, environment.rayDepth, environment.rayBranching);
                        localAverage = Geometry.plus(localAverage, result.color);
                    }
                    averageColor = Geometry.plus(averageColor, Geometry.scale(1/environment.samplesPerPoint, localAverage));
                    count++;
                }
                averageColor = Geometry.scale(1/count, averageColor);
                // convert x and y to image coordinates 
                const realX = Math.floor((x+1)/2 * camera.resolution);
                const realY = Math.floor((y+1)/2 * camera.resolution);
                pixelBuffer[bufferIndex][0] = realX;
                pixelBuffer[bufferIndex][1] = realY;
                pixelBuffer[bufferIndex][2] = averageColor;
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