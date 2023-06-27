Installation of VM2D
====================

.. Brief instructions for installing the **VM2D** code are given
.. at the project GitHub-page <https://github.com/vortexmethods/VM2D>`_

This section describes the information required to install **VM2D** on your computer.

Installation is possible both on computers under *Linux* and *Windows*.


Necessary software 
------------------

Compiling **VM2D** requires the following software to be installed:

* `cmake <https://cmake.org/>`_ automation system for building software from the source code,
* C++ compiler supporting *OpenMP* technologie and the *С++11* standard,
* *MPI* libraries, for example, `Open MPI <https://www.open-mpi.org/>`_ or `Microsoft MPI <https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi>`_,
* `Eigen <http://eigen.tuxfamily.org>`_ libraries (not necessarily, there are source codes in the "include" folder),
* if there is a *Nvidia* GPU supporting `CUDA <https://ru.wikipedia.org/wiki/CUDA>`_ technologie, `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`_ is required for its usage.

`Paraview <https://www.paraview.org/>`_ is useful for viewing the simulation results.


Downloading source codes
------------------------

To install the **VM2D** code on your computer, you need to download the source codes.
If *Git* is installed on your computer, just execute the command

      ``git clone https://github.com/vortexmethods/VM2D.git VM2D``

As a result, a subfolder ``VM2D`` will be created in the current folder and all files from the repository will be loaded into it, including source codes, examples of tasks to be solved, scripts, etc.	  
	  
If *Git* is not installed on your computer, the source codes can be downloaded directly from the `GitHub-page of the VM2D project <https://github.com/vortexmethods/VM2D>`_


Preparing for installation
--------------------------
  
Preparing for compiling the source codes involves creating a folder ``build`` in the directory with the downloaded source codes, going to this folder and executing the command

      ``cmake ..``
	  
If necessary, you should set the necessary keys for configuring the compilers used, specifying compilation options, etc. It may also require some modification of the ``CMakeLists.txt`` file containing the *CMake* settings and located in the project root directory.

In particular, to prepare source codes for compiling them on *Windows* using *MS Visual Studio*, you should, depending on the version, use one of the following commands (the ``Win64`` option is required to use the ability to perform calculations on graphics cards *Nvidia CUDA*, for *Visual Studio 2019* it is enabled by default for 64-bit systems)

      ``cmake -G"Visual Studio 14 2015 Win64" ..``
	  
      ``cmake -G"Visual Studio 15 2017 Win64" ..``
	  
      ``cmake -G"Visual Studio 16 2019" ..``
	  
	  
If you use a compiler on *Windows* that differs from the default one (as a rule, the *MSVC* compiler built into *MS Visual Studio*), for example, the *Intel* compiler, when preparing the source codes for compilation, you must specify, depending on version, key (note that Intel C++ Compiler 19 is integrated into Visual Studio 2019 only from Upd.4 version)
	  
      ``-T"Intel C++ Compiler <ver>"``
	  
On * Linux *, if you want to use some other compiler instead of *gcc/g++*, for example, *Intel* compiler (*icpc*), you should execute the command

	  ``CXX=icpc camke ..``	  
  

Source code compilation
-----------------------	  
	  
Compilation procedure depends on operating system used.

On *Windows*, as a rule, a project is created using *CMake* for its further opening and compilation in *MS Visual Studio* (see above). 

On ``Linux``, you can just execute the command

      ``make``

If compilation (on *Linux*) fails with an error related to the lack of necessary *MPI* libraries, you can try to clear the ``build`` folder, go to it, and instead of the ``cmake..`` command, execute one of the following commands:

      ``CXX=mpiCC cmake ..``
	  
	  ``CXX=mpiicpc cmake ..``
	  	  
and after its successful execution, repeat the compilation using the ``make`` command. 

	  
	 
Computations using GPUs (Nvidia CUDA)
-------------------------------------
	 
If the `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>` _ is installed, during the ``cmake`` execution the program will be automatically configured to  calculations using GPUs and *NVidia CUDA* technology.


When performing calculations on a node with several GPUs, by default all calculations will be performed on the device with the index ``#0``. To avoid this, you need to open the source code file ``VM2D/src/VM2D/Gpu2D/Gpu2D.cpp``, read the comments written in it and comment out the required line.
	 
If the *CUDA Toolkit* is installed, but it is planned to perform calculations on CPUs without using the GPU capabilities, then in the CMakeLists.txt file you need to comment out the line

      ``add_definitions(-DUSE_CUDA)``
	  
by adding symbol ``#`` to get
	  
	  ``#add_definitions(-DUSE_CUDA)``	

