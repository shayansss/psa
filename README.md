# PSA: Pre-Stress Algorithm
This repository contains the core data corresponding to our accepted manuscript, titled "[Large-Scale Finite Element Modeling of Pre-stress in Articular Cartilage](http://dx.doi.org/10.1007/978-3-031-55315-8_12)".

Pre-stressing in multi-phasic models of articular cartilage affects its soft structure, subsequently altering the material due to fluid diffusion within the solid structure. Given that the initial state of the cartilage prior to the pre-stressed state (where finite element analysis typically begins) is unknown, we have developed a unified optimizer to find that state. This optimizer has proven suitable for large-scale pre-stressing analysis. The value of such an algorithm lies in its efficiency; it is not overly time-consuming and allows us to include the crucial osmotic pressure phase (which causes pre-stressing) alongside the other phases at a large scale.

## Citation
If this research data is useful for your work, kindly please consider citing our work ([DOI](http://dx.doi.org/10.1007/978-3-031-55315-8_12)):

```
@InProceedings{sajjadinia2024a,
author="Sajjadinia, Seyed Shayan and Carpentieri, Bruno and Holzapfel, Gerhard A.",
editor="Skalli, Wafa and Laporte, S{\'e}bastien and Benoit, Aur{\'e}lie",
title="Large-Scale Finite Element Modeling of Pre-stress in Articular Cartilage",
booktitle="Computer Methods in Biomechanics and Biomedical Engineering II",
year="2024",
publisher="Springer Nature Switzerland",
address="Cham",
pages="105--112",
isbn="978-3-031-55315-8",
doi="10.1007/978-3-031-55315-8\_12"
}
```

This study is also based on our previous work ([DOI](https://doi.org/10.1016/j.jmbbm.2020.104203)), which can cited as follows:

```
@article{sajjadinia2021a,
author="Sajjadinia, Seyed Shayan and Carpentieri, Bruno and Holzapfel, Gerhard A.",
title="A Backward Pre-stressing Algorithm for Efficient Finite Element Implementation of In Vivo Material and Geometrical Parameters into Fibril-reinforced Mixture Models of Articular Cartilage",
journal="Journal of the Mechanical Behavior of Biomedical Materials",
volume="114",
pages="104203",
year="2021",
issn="1751--6161",
doi="10.1016/j.jmbbm.2020.104203"
}
```

## Dependencies
- Visual Studio 2019 or later (Community edition or higher)
- Intel® Parallel Studio XE 2020 (Update 4 or later)
- Abaqus 2021 (full version)

## Installation
Firstly, you need to install all the dependencies, which have been successfully tested on Windows 10, and then link them. The process is described in [this tutorial](http://dx.doi.org/10.13140/RG.2.2.33539.32800), which uses a similar version. Afterward, download this repository into your system. Given that the address of the root directory of the local repository on your system may differ from the one set in the code, you will need to change it to your local address. You can do this by modifying the values passed to the `os.chdir` function in the Python file and the `FilLoc` variable in the Fortran file. In the latter case, if the length of the new address is changed, you will also need to correct the length defined after `CHARACTER FilLoc` that is the length + 8.

## How to Run
Execute the Python scripts within Abaqus. This process first initializes the state variables for each point, which are then saved by data sharding to speed up access in Fortran. Finally, it implements the pre-stressing algorithm and generates the results.

Enjoy!
