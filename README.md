# PSA: Pre-Stress Algorithm
This repository contains the core data corresponding to our accepted manuscript, titled "Large-Scale Finite Element Modeling of Pre-Stress in Articular Cartilage".

Pre-stressing in multi-phasic models of articular cartilage affects its soft structure, subsequently altering the material due to fluid diffusion within the solid structure. Given that the initial state of the cartilage prior to the pre-stressed state (where finite element analysis typically begins) is unknown, we have developed a unified optimizer to find that state. This optimizer has proven suitable for large-scale pre-stressing analysis. The value of such an algorithm lies in its efficiency; it is not overly time-consuming and allows us to include the crucial osmotic pressure phase (which causes pre-stressing) alongside the other phases at a large scale.

## Dependencies
- Visual Studio 2019 or later (Community edition or higher)
- Intel® Parallel Studio XE 2020 (Update 4 or later)
- Abaqus 2021 (full version)

## Installation
Firstly, you need to install all the dependencies, which have been successfully tested on Windows 10, and then link them. The process is described in [this tutorial](http://dx.doi.org/10.13140/RG.2.2.33539.32800), which uses a similar version. Afterward, download this repository into your system. Given that the address of the root directory of the local repository on your system may differ from the one set in the code, you will need to change it to your local address. You can do this by modifying the values passed to the `os.chdir` function in the Python file and the `FilLoc` variable in the Fortran file. In the latter case, if the length of the new address is changed, you will also need to correct the length defined after `CHARACTER FilLoc` that is the length + 8.

## How to Run
Execute the Python scripts within Abaqus. This process first initializes the state variables for each point, which are then saved by data sharding to speed up access in Fortran. Finally, it implements the pre-stressing algorithm and generates the results.

Enjoy!
