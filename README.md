# ray-tracer
Building a custom path tracer in Javascript for educational purposes.

Current Features:
- 3D Primitives: Spheres, Planar Polygons
- Materials: variable glossy
- CPU multithreading using Web Workers
- Unidirectional path tracing

To Do:
- Bidirectional path tracing with multiple importance sampling
- transparent and frosty materials, fog
- GPU acceleration, if possible

Unidirectional path tracing output:
![Path traced sphere](https://i.imgur.com/PXCq98s.png)

Update: rewritten with Rust. Now with Bounding Volume Hierarchy for faster intersection!
![Path traced hippo](https://i.imgur.com/LVHLGRO.png)
