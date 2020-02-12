# 2DFluidSim

[![youtube](youtube_image.png)](https://www.youtube.com/watch?v=BLFr2GfCn5c)

## Introduction

This code accompanied our submission to Proceedings of the ACM (PACM) in May 2019 (paper included above). This repo is divided into two components: a C++ solver and a Houdini plugin.

The C++ solver solves the incompressible Euler equations in spherical coordinates on a staggered grid defined on the sphere. The goal of this method is to accurately model fluid behavior regardless of the underlying grid definition (i.e, the result should be grid-agnostic). The spherical coordinate system is simple and convenient, but it introduces polar singularities that our method mitigates with an appropriate choice of boundary conditions and approximation methods. More technical details can be found in the video posted above and in the paper included in this repository.

The Houdini plugin enables the authoring of fluid simulations on a sphere. Perhaps the most natural use case is to use this tool to generate detailed renderings of soap bubbles, but as you can see in the render below, creative users will find many applications for its use.

![lollipop](lolipop_real_4k.jpg)

## Installation & Usage

### Hardware Requirements

- Memory

  64 bit operating system
  8 GB or higher and 64 GB strongly recommended.

- Processor

  64-bit AMD or Intel processor (requires SSE or higher)

- Graphics Card

  NVidia GeForce GPUs: 390.77 or higher
  NVidia Quadro GPUs: 390.77 or higher
  AMD: 17.Q1 or higher

### Software Requirements

Windows 8 or 10 (64 bit) and the latest version of Houdini.

### Setup 

Environment settings: CUDA 9.1 and Visual Studio 2015 are required.

Have these environment variables set:
H16_VERSION				The version of your Houdini distribution. e.g. 16.5.405
H16_PATH					The Path to the Houdini installation. e.g. C:/Program Files/Side Effects Software/Houdini %H16_VERSION%
CUSTOM_DSO_PATH		The location you intend to install the plugin. e.g. D:/houdiniPlugin
HOUDINI_DSO_PATH	Set this value to %CUSTOM_DSO_PATH%;&

Open Kamino.sln, set solution configuration and solution platform to Release64+x64. Build the solution.

When it’s done, copy the executable KaminoCoreSolver.exe in x64/Release to your $CUSTOM_DSO_PATH. If you don’t have opencv2.4 installed on your computer (and in your $PATH variable), you also need to copy all the opencv*.dll (there are 6 of them) libraries to the $CUSTOM_DSO_PATH.

Make two folders in the $CUSTOM_DSO_PATH, one called “output”, and the other called “particles”.

### User Workflow

In your Houdini scene, create an instance node in the network editor. Double click the node, press tab, and type “Kamino” to search for the Kamino node.

Create a Kamino node. You should now have access to the following parameters:

- radius: The radius of the sphere.
- nTheta: The resolution of the simulation grid. The total amount of grids would be 2*nTheta*nTheta.
- particleDensity: The maximum amount of particles per grid.
- dt: The time step.
- DT: The frame rate. It should be assigned to 1.0/FPS.
- frames: The amount of frames to be generated.
- spin: The spinning speed component along the equator direction.
- uphi, utheta, vphi, vtheta: Parameters used to initialize the velocity field distribution. The higher the parameter, the more turbulent the velocity will be.
- densityImage: The full path to an image that defines the initial density distribution of the particles. The image would be resized to fit the resolution of the simulation grid.
- solidImage: The full path to an image that specifies whether a cell should be solid or fluid. White cells (greyscale value > 128) would be solid grids, and black ones would be fluid grids.
- colorImage: The full path to an image that would assign colors to the particles. Extremely useful for creating various patterns flowing on the sphere, for example, Latte art patterns.

Press “Run” button to run the simulation. Wait for the pop up window to finish its job.
Load the velocity field and particles from $HIP/output and $HIP/particles and do whatever you want to them. Typical usage would be to feed the output of the particle node to OpenVDB nodes to generate a volumetric fluid.

### Content Creation

In true Houdini fashion, users can leverage the physically-based solver we’ve developed to produce visually compelling results.

Users can scrub through frames to see the flow patterns developing on the sphere’s surface. We typically recommend that users first output a few frames on a lower grid resolution to get a feel for the parameters they’ve specified. After a desirable look is achieved, users should increase the grid resolution, increase the particle density (if so desired), and output the required number of frames.
