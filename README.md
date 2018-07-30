# 2DFluidSim

----- Hardware Requirements -----

Memory
64 bit operating system
8 GB or higher and 64 GB strongly recommended.

Processor
64-bit AMD or Intel processor (requires SSE or higher)
Graphics Card
NVidia GeForce GPUs: 390.77 or higher
NVidia Quadro GPUs: 390.77 or higher
AMD: 17.Q1 or higher

----- Software Requirements -----

Windows 8 or 10 (64 bit) and the latest version of Houdini.

----- Setup -----

Note: We will try and release a bootstrapping shell code to do these steps for you in the future.

Environment settings: CUDA 9.1 and Visual Studio 2015 are required.

Have these environment variables set:
H16_VERSION				The version of your Houdini distribution. e.g. 16.5.405
H16_PATH					The Path to the Houdini installation. e.g. C:/Program Files/Side Effects Software/Houdini %H16_VERSION%
CUSTOM_DSO_PATH		The location you intend to install the plugin. e.g. D:/houdiniPlugin
HOUDINI_DSO_PATH	Set this value to %CUSTOM_DSO_PATH%;&

Clone the repository from GitHub: ProjectKamino

Open Kamino.sln, set solution configuration and solution platform to Release64+x64. Build the solution.

When it’s done, copy the executable KaminoCoreSolver.exe in x64/Release to your $CUSTOM_DSO_PATH. If you don’t have opencv2.4 installed on your computer in your $PATH variable, you also need to copy all the opencv*.dll (there are 6 of them) to the $CUSTOM_DSO_PATH.
Make two folders in the $CUSTOM_DSO_PATH, one called “output”, and the other called “particles”.

----- User Workflow -----

In your Houdini scene, create an instance node in the network editor. Double click the node, press tab, and type “Kamino” to search for the Kamino node.

Create a Kamino node. You should now have access to the following parameters:

*radius*: The radius of the sphere.

*nTheta*: The resolution of the simulation grid. The total amount of grids would be 2*nTheta*nTheta.

*particleDensity*: The maximum amount of particles per grid.

*dt*: The time step.

*DT*: The frame rate. It should be assigned to 1.0/FPS.

*frames*: The amount of frames to be generated.

*spin*: The spinning speed component along the equator direction.

*uphi, utheta, vphi, vtheta*: Parameters used to initialize the velocity field distribution. The higher the parameter, the more turbulent the velocity would be.

*densityImage*: The full path to an image that defines the initial density distribution of the particles. The image would be resized to fit the resolution of the simulation grid.

*solidImage*: The full path to an image that specifies whether a cell should be solid or fluid. White cells (greyscale value > 128) would be solid grids, and black ones would be fluid grids.

*colorImage*: The full path to an image that would assign colors to the particles. Extremely useful for creating various patterns flowing on the sphere, for example, Latte art patterns.

Press “Run” button to run the simulation. Wait for the pop up window to finish its job.
Load the velocity field and particles from $HIP/output and $HIP/particles and do whatever you want to them. Typical usage would be feeding the output of the particles to OpenVDB nodes and generate fluid.

----- Content Creation -----

In true Houdini fashion, users can leverage the physically-based solver we’ve developed to produce visually compelling results.

Users can scrub through frames to see the flow patterns developing on the sphere’s surface. We typically recommend that users first output a few frames on a lower grid resolution to get a feel for the parameters they’ve specified. After a desirable look is achieved, users should increase the grid resolution, increase the particle density (if so desired), and output the required number of frames.
