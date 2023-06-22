# Large-scale pre-stress in multi-phasic cartilage
This repository contains the core data of our under-review manuscript, entitled "Large-Scale Finite Element Modeling of Pre-Stress in Articular Cartilage".

Pre-stressing in multiphasic models of articular cartilage causes its soft structure to change, consequently changing the material (due to the fluid diffusion inside the solid structure). Since the initial state of cartilage before pre-stressing (where finite element analysis starts) is unknown, a unified optimizer is developed, which proved to be suitable in large-scale pre-stressing analysis. The importance of such a model is that it is not very time-consuming and allows us to include the important osmotic pressure phase (that causes pre-stressing) next to the other phases.

## Dependency
- Visual Studio = 2019 (at least the Community edition)
- Intel® Parallel Studio XE = 2020 with Update 4
- Abaqus/CAE = 2021 (full version)

## Installation
You should install all the dependencies (tested successfully on Windows 10) and link them. This is described (using another similar version) in [this tutorial](http://dx.doi.org/10.13140/RG.2.2.33539.32800). Then, download this repository to your system. As the address of the root directory of the local repository in your system may not be the same as the one set in the code, you should change it to your local address by correcting the values passed to: 1) the `os.chdir` function inside Python file; 2) the `FilLoc` variable in the Fortran file (in this case, if the length of the new address is changed, you should also correct the length assigned to `LenFil` defined after `CHARACTER FilLoc`).

## How to run
Run the Python scripts inside Abaqus. It first initializes state variables, in the values in shraded text files, and then implementes pre-stressing. Finally, it generates the convergence plot.

Enjoy!
