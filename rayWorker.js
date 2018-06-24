
const Geometry = require('./geometry.js');

module.exports = function(self){

    let featureCollection = null;
    let lightSources = [];
    let environment = {};

    self.addEventListener('message',function(e){
        let type = e.data[0];
        if(type === 'overhead'){
            featureCollection = e.data[1];
            lightSources = e.data[2];
            environment = e.data[3];
            // reattach methods that were stripped during serialization
            featureCollection.__proto__ = Geometry.FeatureCollection.prototype
            featureCollection.deserialize()
        }else{
            let rays = e.data[1];
            console.log(rays[0])
            for(let k = 0; k < rays.length; k++){
                // deserialize the rays 
                rays[k].__proto__ = Geometry.Ray.prototype;
                rays[k].deserialize();
                rays[k] = Geometry.traceRay(featureCollection, environment, lightSources, rays[k], environment.maxRecursion);
            }

            self.postMessage([rays]);
        }
    });
}