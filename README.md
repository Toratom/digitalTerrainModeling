# digitalTerrainModeling

## The Project

This project was realized at Télécom Paris in the framework of the IGR205 course. The goal was to implement several methods to generate and simulate digital terrains.
The terrain is modeled by several heightmaps, one per material layer.
During this project we also had to implement a simple water simulation model, we used the pipe model.

We implemented the following methods:
 * Terrain generation from noise.
 * Terrain generation by the fault algorithm
 * Simulation of thermal errosion
 * Simulation of hydraulic errosion

For more information, please refer to the file: ```report.pdf```

This project was done in C++ and OpenGL


## Organization of the Repository

The repository is organised in two branches : ```main``` and ```GPU```. Indeed, we implemented some of the previous algorithms on GPU with compute shaders. The GPU version also offers a mode to edit the heightmap.
