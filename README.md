# Digital Terrain Modeling

## Overview

This project was realized at Télécom Paris in the framework of the IGR205 course. The goal was to implement several methods to generate and simulate digital terrains.
The terrain is modeled by several heightmaps, one per material layer.
During this project we also had to implement a simple water simulation model, we used the pipe model.

We implemented the following methods:
 * Terrain generation from noise
 * Terrain generation by the fault algorithm
 * Simulation of thermal erosion
 * Simulation of hydraulic erosion

This project was done in C++ and OpenGL. The two simulations were implemented on CPU and GPU.
For more information, please refer to the file: ```report.pdf```

https://github.com/Toratom/digitalTerrainModeling/assets/74663696/49fe69ce-9c8d-48c3-bf4b-a55384df2670

## Setup

The repository is organized into two branches: ```main``` and ```GPU```. Indeed, we implemented some of the previous algorithms on GPU with compute shaders.
The GPU version also offers a mode to edit the heightmaps, which can be started by clicking on the key "E".
For both branches, the installation process is the same. The installation process is described here for Windows, it requires [CMake](https://cmake.org/download/) and [Visual Studio](https://visualstudio.microsoft.com/fr/free-developer-offers/)
  * Clone the repository.
  * Add a folder ```build``` at the root of the project.
  * Open CMake:
    * Select as "where is the source code", the root folder of the project i.e. the one containing the file ```CMakeLists.txt```.
    * Select as "where to build the binaries", the build folder previously created.
    * Click a first-time "Configure", then a second time. Then click "Generate". And finally, "Open Project".
  * If everything goes well, Visual Studio automatically opens. In Visual Studio:
    * In the "Solution Explorer", select the "Solution 'projetIGR205'".
    * Below in the "Properties" panel, click on the wrench to open the property page.
    * In the "Single startup project" drop-down list, select "projetIGR205".
  * Now, you should be able to compile the project.

## Contributors

[Lucas Teissier](https://github.com/LucasTsr) - [Octave Le Tullier](https://github.com/OctaveLT) - [Thomas Poyet](https://github.com/Toratom)
