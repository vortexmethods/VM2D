Working directory content
=========================

To get started with **VM2D** you need to prepare *working directory*.
The working directory must contain 4 text files:

*     ``problems``
*     ``defaults``
*     ``switchers``
*     ``mechanics``

Examples of such files are placed in the *run* folder of the source repository.

It is possible to use \texttt{VM2D} for solution of one particular problem as well as for solution of the set of similar (or not similar) problems. Every problem to be solved is indicated by some *label* (name), which is a text string, and for every problem in the working directory separate subdirectory should be created with the same name as problem's *label*.


File ``problems``
-----------------

Labels of all problems that are supposed to be solved should be listed in the file ``problems`` placed in the working directory. Next to the label of each task, the parameters of its start are specified in parentheses, and some parameters may also be listed that will be passed to the *passport* of the task. The typical structure of the ``problems`` file is shown below.

.. code-block:: c
   :caption: problems
   :name: problems
	
   /*--------------------------------*- VM2D -*-----------------*---------------*\
   | ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
   | ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
   | ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
   |  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
   |   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
   |                                                                             |
   | Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
   *-----------------------------------------------------------------------------*
   | File name: problems                                                         |
   | Info: Problems to be solved by using VM2D                                   |
   \*---------------------------------------------------------------------------*/
   problems = {
       wing00deg(np = 1, angle =  0, tau = 1.5e-2),
       wing05deg(np = 2, angle =  5, tau = 1.5e-2),
       wing10deg(np = 2, angle = 10, tau = 2.0e-2)
   };

In this example, three problems are listed: ``wing00deg``, ``wing05deg`` and ``wing10deg``, for each problem starting parameters are specified in parentheses. 

There is no mandatory parameters, thus, in the simplest case, parentheses can be empty. There are two parameters that influence the starting procedure:

+-----------+-------------+------------------------------------+
| Parameter | Default     | Description                        |
| name      | value       |                                    |
+===========+=============+====================================+
| pspfile   | passport    | name of the file where             |
|           |             | *passport* of the solved problem   |
|           |             | is contained                       |
+-----------+-------------+------------------------------------+
| np        | 1           | number of MPI-processes that       |
|           |             | will be created for problem solving|
+-----------+-------------+------------------------------------+

When starting calculations on a multi-node computing cluster, when each node is a system with shared memory, it is advisable to launch only one *MPI*-process per node, since all cores in the node will be automatically involved using the *OpenMP* technology.

All other parameters (in the above example these include *angle* and *tau*) can be specified by user arbitrarily and have any type - integer, floating point, boolean, string, and also a list. Their setting can be convenient for further use inside *passports* of problems.




File ``defaults``
-----------------

The **VM2D** code implements three priority levels for the values of various parameters.

1. **Low-level priority.** For some parameters, default values are specified directly in the source code of the program. These include, in particular, the above parameters *pspfile* and *np*.

2. **Mid-level priority.** User can arbitrarily set the default values for some parameters in the ``defaults`` file, so as not to explicitly indicate them in the passports of specific problems. The file, like the others, is organized as a dictionary containing ``key = value`` pairs. An example of the structure of this file is shown below.

.. code-block:: c
   :caption: defaults
   :name: defaults
	
   /*--------------------------------*- VM2D -*-----------------*---------------*\
   | ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
   | ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
   | ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
   |  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
   |   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
   |                                                                             |
   | Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
   *-----------------------------------------------------------------------------*
   | File name: defaults                                                         |
   | Info: Default values for various parameters                                 |
   \*---------------------------------------------------------------------------*/
   
   angle = 0.0;
   scale = 1.0;
   basePoint = {0.0, 0.0};

   linearSystemSolver = linearSystemGauss;
   velocityComputation = velocityBiotSavart;
   panelsType = panelsRectilinear; 
   boundaryConditionSatisfaction = boundaryConstantLayerAverage;
   fileType = text;

   mechanicalSystem = mechanicsRigidImmovable;
   np = 1;

   airfoilsDir = "../settings/airfoils/";
   wakesDir = "../settings/wakes/";


3. **High-level priority.** Direct indication of parameter values in the problem's passport has the highest priority, while the default values for these parameters are ignored.




File ``switchers``
------------------

To simplify the readability of the values of some keys, they can be specified in the corresponding files using verbal expressions. The value of the key *mechanicalSystem* given, for example, in the previous section, specified by the word *mechanicsRigidImmovable*, obviously means that a fixed and non-deformable airfoil is used as a mechanical system.

At the same time, the **VM2D** source code assumes that such keys correspond to integer values. The correspondence of word expressions to integer values is specified in the file ``switchers``, an example of its structure is given below.

.. code-block:: c
   :caption: switchers
   :name: switchers

   /*--------------------------------*- VM2D -*-----------------*---------------*\
   | ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
   | ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
   | ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
   |  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
   |   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
   |                                                                             |
   | Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
   *-----------------------------------------------------------------------------*
   | File name: switchers                                                        |
   | Info: Verbal notation for various parameters                                |
   \*---------------------------------------------------------------------------*/

   Impulse = 0;  //Impulsively started flow
   RampLin = 1;  //Flow acceleration according to linear law
   RampCos = 2;  //Flow acceleration according to cosine law

   panelsRectilinear = 0;   //Rectilinear panels
   panelsCurvilinear = 1;   //Curvilinear panels
   
   text = 0;                //vtk-files storage in text format
   binary = 1;              //vtk-files storage in binary format

   linearSystemGauss = 0;      //Linear system solving by Gaussian elimination
   linearSystemGMRES = 1;      //Linear system solving by GMRES
   linearSystemBiCGStab = 2;   //Linear system solving by BiCGStab
   linearSystemFMM = 3;        //Linear system solving by fast multipole method

   mechanicsRigidImmovable = 0;   //Immovable non-deformable body
   mechanicsRigidGivenLaw = 1;    //Non-deformable body moving according to given law
   mechanicsRigidOscillPart = 2;  //Non-deformable body with elastic constrain, partitioned approach
   
   
Отметим, что приведенные параметры даны исключительно в качестве примера, не весь функционал может быть реализован в текущей версии программы (и наоборот, реальный файл может содержать намного большее количество различных параметров).
Note that the given parameters are listed as an example, not all functionality can be implemented in the current version of the program (and vice versa, a real file can contain a much larger number of different parameters).

   

File ``mechanics``
------------------

Since **VM2D** allows solving a wide class of problems, including FSI problems, the specific type of the problem being solved should be indicated by specifying the type of *mechanical system* used. The available (software implemented) mechanical systems are listed in the file ``mechanics``, where they are assigned short label names.

An example of the structure of such a file is given below.


.. code-block:: c
   :caption: switchers
   :name: switchers

   /*--------------------------------*- VM2D -*-----------------*---------------*\
   | ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
   | ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
   | ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
   |  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
   |   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
   |                                                                             |
   | Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
   *-----------------------------------------------------------------------------*
   | File name: mechanics                                                        |
   | Info: Dictionary with mechanical systems description                        |
   \*---------------------------------------------------------------------------*/

   mech0 = mechanicsRigidImmovable();
   mech1 = mechanicsRigidGivenLaw();
   mech2 = mechanicsRigidOscillPart(sh={0, $sh, 0}, m=$m);
   
There are two types of mechanical systems: *mechanicsRigidImmovable* and *mechanicsRigidGivenLaw* have no parameters (the law of motion of a body in the second case is set directly inside the **VM2D** code), while the third mechanical system is determined by two parameters --- dimensionless oscillation frequency with elastic constraints, which is denoted by *sh*, and body mass *m*. Both of them are set implicitly, but through the **$** symbol, which means dereferencing the parameter. Accordingly, when choosing this type of mechanical system for a specific problem being solved, the user must specify numerical values ​​for the *sh* and *m* parameters.





