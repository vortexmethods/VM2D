Problem configuration
=====================


Task passport ``passport``
--------------------------


In the previous section, example of file ``problems`` was considered, where it was supposed to solve three similar problems, which, as the names suggest, consist in modeling the flow around a wing installed at some angle of attack to the incoming flow. We will present this file here again.

.. code-block:: c
   :caption: problems
   :name: problemsRepeat
	
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
   
To perform these three simulations, subdirectories must be created in the *working directory* with names corresponding to the task labels: ``wing00deg``, ``wing00deg`` and ``wing00deg``.

Inside each of these subdirectories, there must be a file with *task passport*, which by default is named ``passport``, but if necessary, it can be renamed (see the previous section).

The file *task passport*, which could be used for these tasks, has the following form.

.. code-block:: c
   :caption: passport
   :name: passport
	
   /*--------------------------------*- VM2D -*-----------------*---------------*\
   | ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
   | ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
   | ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
   |  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
   |   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
   |                                                                             |
   | Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
   *-----------------------------------------------------------------------------*
   | File name: passport                                                         |
   | Info: Parameters of the problem to be solved                                |
   \*---------------------------------------------------------------------------*/

   //Physical Properties
   rho = 1.0;          // flow density
   vInf = {1.0, 0.0};  // incident flow velocity
   //vRef = 1.0;       // reference velocity magnitude, isn't used here
   nu = 0.001;         // kinematic viscosity coefficient

   //Time Discretization Properties
   timeStart = 0.0;    // *physical time at which the simulation starts
   timeStop = 50.0;    //  physical time at which the simulation stops
   dt = $tau;          //  the specific value is set in a batch file, e.g., 1.5e-2
   accelVel = Impulse; // *RampLin(T) and RampCos(T) modes are supported

   fileType = text;      // *text or binary format for output vtk-files
   saveVtx = 100;        // *frequency (in steps) of vortices storing
   saveVP = 100;         // *the same velocity & pressure
   nameLength = 5;       // *number of digits in file names

   //Wake Discretization Properties
   eps = 0.0015;       //  vortex core smoothing radius (should be small)
   epscol = 0.0040;    //  vortex particles merging distance
   distFar = 20.0;     // *the distance of the vortex wake cropping
   delta = 5e-6;       // *distance from airfoil to the generating vortices
   vortexPerPanel = 1; // *minimal number of vortices generating over each panel
   maxGamma = 1.0e-4;  // *maximal vortex particle strength

   //Numerical Schemes
   linearSystemSolver = linearSystemGauss;  // *fast methods are under constr.
   velocityComputation = velocityBiotSavart;// *Barnes-Hut method is supported
   panelsType = panelsRectilinear;          // *curvilin. panels are under constr.
   boundaryConditionSatisfaction = boundaryLinearLayerAverage; // *T^1 scheme

   //Files and parameters
   airfoilsDir = "../settings/airfoils/"; // *path to discretized airfoils
   wakesDir    = "../settings/wakes/";    // *path to vortex wakes

   airfoil = {
       "naca0012"(                                          // *file name
           basePoint = {1.0, 0.0}, scale = 1.0,             // *geometry
           inverse = false,                                 // *external flow
           angle = $angle,                                  // *AoA
           mechanicalSystem = mechanicsRigidImmovable()) }; // *mechanics

   fileWake = { };   // *previously stored vortex wake can be loaded
   fileSource = { }; // *positions and intensities of point sources/sinks
   

The meaning of most of the parameters should be clear from their names and the short comments given. An asterisk next to a comment denotes those options with *low-level priority* defaults, and, therefore, they can be omitted if the default values are appropriate.
  
The *low-level priority* default values that are set directly in the **VM2D** code are listed in the following table.

+-------------------+-----------------+--------------------------------------+
| Parameter         | Default         | Description                          |
| name              | value           |                                      |
+===================+=================+======================================+
| timeStart         | 0.0             | physical time at which               |
|                   |                 | the simulation starts                |
+-------------------+-----------------+--------------------------------------+
| accelVel          | RampLin(1.0)    | linearly flow acceleration from      |
|                   |                 | zero in 1 second                     |
+-------------------+-----------------+--------------------------------------+
| fileType          | text            | text format of output vtk files      |
+-------------------+-----------------+--------------------------------------+
| saveVtx           | 100             | saving of vortex particles positions |
|                   |                 | at every 100 time steps              |
+-------------------+-----------------+--------------------------------------+
| saveVP            | 0               | do not calculate velocity and        |
|                   |                 | pressure in the specified points     |
|                   |                 | (step should be set for saving)      |
+-------------------+-----------------+--------------------------------------+
| nameLength        | 5               | length of file name                  |
+-------------------+-----------------+--------------------------------------+
| distFar           | 10.0            | distance at which vortex particles   |
|                   |                 | are removed form the simulation      |
+-------------------+-----------------+--------------------------------------+
| delta             | 1.0e-5          | small distance at which vortex       |
|                   |                 | particles are placed over            |
|                   |                 | the airfoil surface line             |
|                   |                 | after being generated                |
+-------------------+-----------------+--------------------------------------+
| vortexPerPanel    | 1               | minimal number of vortex particles   | 
|                   |                 | generated from one panel of the      |
|                   |                 | airfoil for every time step          |
+-------------------+-----------------+--------------------------------------+
| maxGamma          | 0.0             | maximal value of vortex particle     |
|                   |                 | circulation is not bounded from above|
|                   |                 | (to bound, the value of acceptable   | 
|                   |                 | circulation should me set)           |
+-------------------+-----------------+--------------------------------------+
| panelsType        |panelsRectilinear| rectilinear panels at the airfoil    |
+-------------------+-----------------+--------------------------------------+
| boundaryCondition | constLayerAver  | the vorticity generated at the       |
|                   |                 | airfoil is represented as a vortex   |
|                   |                 | sheet with piecewise-constant        |
|                   |                 | vorticity (across panels)            |
+-------------------+-----------------+--------------------------------------+
| linearSystemSolver|linearSystemGauss| solving of SLAE (discrete analogue of|
|                   |                 | the boundary integral equation) by   |
|                   |                 | Gaussian elimination (LU-decompos.)  |
+-------------------+-----------------+--------------------------------------+
|velocityComputation|biotSavart       | direct computation of vortex elements|
|                   |                 | velosities (Biot-Savart law)         |
+-------------------+-----------------+--------------------------------------+
| airfoilsDir       | ../settings/    | path to files with airfoils          |
|                   | airfoils/       | geometries                           |
+-------------------+-----------------+--------------------------------------+
| wakesDir          | ../settings/    | path to files with vortex wakes      |
|                   | wakes/          |                                      |
+-------------------+-----------------+--------------------------------------+
| fileWake          | ""              | vortex wake is not loaded            |
|                   | (empty string)  |                                      |
+-------------------+-----------------+--------------------------------------+
| fileSource        | ""              | sources are not loaded               |
|                   | (empty string)  |                                      |
+-------------------+-----------------+--------------------------------------+
| airfoil           | {}              | airfoils are not loaded              |
|                   | (empty list)    |                                      |
+-------------------+-----------------+--------------------------------------+
| basePoint         | {0.0, 0.0}      | coordinates of the base point        |
|                   |                 | (center) of the downloading airfoil  |
+-------------------+-----------------+--------------------------------------+
| scale             | 1.0             | scaling factor of the airfoil        |
+-------------------+-----------------+--------------------------------------+
| angle             | 0.0             | angle of incidence                   |
+-------------------+-----------------+--------------------------------------+
| inverse           | false           | means modeling of the flow external  |
|                   |                 | to the given airfoil                 |
+-------------------+-----------------+--------------------------------------+
| mechanicalSystem  | rigidImmovable  | immovable non-deformable airfoil     |
+-------------------+-----------------+--------------------------------------+

These parameters can be overridden, if necessary, in the ``defaults`` file, and default values for other parameters can also be set there.

Here are additional explanations for some parameters:

*  *vRef* sets the speed scale, must be specified in tasks where the incident flow is absent, i.e. *vInf = {0.0, 0.0}*, since it acts as a velocity scale when calculating dimensionless parameters; if there is an incident flow, the vector module *vInf* is automatically used as *vRef*;
 
*  *accelVel* defines the way to "accelerate" the incoming stream from zero to the value *vInf*, can take the following values:
    
    *  *Impulse* - the stream starts instantly,
    *  *RampLin(T)* - linear increase of the flow velocity from zero to *vinf* in *T* seconds,
    *  *RampCos(T)* - increase of the flow velocity from zero to *vinf* in *T* seconds according to cosine law;
    
*   *saveVtx* and *saveVP* determine the frequency of saving the vortex wake to files, as well as the velocities and pressures calculated at specified points; setting zero values means that there is no need to save the corresponding data (velocities and pressure are not calculated in this case);
 
*   *maxGamma* determines the maximum value of the circulation of vortex particles generated on the airfoil, as well as resulting from the restructuring of the vortex wake (in the latter case, when moving away from the airfoil, this value automatically increases slightly); setting a zero value means that the circulation of vortex particles is not boubded from above.


The *airfoil* parameter specifies the name of the file containing information about the airfoil geometry. Files with airfoils are stored in a directory whose name is determined by the *airfoilsDir* parameter (in the repository this is the directory `` settings/airfoil``); the contents of these files are the positions of the vertices of the airfoil vertexes counterclockwise. Their Cartesian coordinates are given as the value of the *r* key and are a list of pairs of numbers. Examples can be found in the repository files.
If the file name does not contain spaces, then it can be specified without quotes, otherwise double quotes are required.

The following parameters are indicated in brackets after the airfoil file name:

+-----------------+-----------------------------------------+
| Parameter       | Description                             |
| name            |                                         |
+=================+=========================================+
| basePoint       | coordinates of the point where airfoil  |
|                 | base point (center) should be placed    |
+-----------------+-----------------------------------------+
| scale           | scaling factor for the airfoil          |
+-----------------+-----------------------------------------+
| inverse         | boolean parameter which takes value     |
|                 | *true* if internal flow is simulated    |
+-----------------+-----------------------------------------+
| angle           | angle of incidence (clockwise)          |
+-----------------+-----------------------------------------+
| mechanicalSystem| the label (type) of the mechanical      |
|                 | system                                  |
+-----------------+-----------------------------------------+


If the flow around the airfoil system is being simulated, they should be sequentially listed in the list corresponding to the *airfoil* key:

.. code-block:: c
   :caption: passport
   :name: passport2

    airfoil = {
        square_160points(
            basePoint = {0.0, 0.0},
            angle = 45.0,
            scale = 1.0
            ),
        circle_200points(
            basePoint = {1.2, -0.2},
            angle = 0.0,
            scale = 0.5
            )

In this example, the interference of two airfoils is considered: a small (reduced by 2 times compared to the base shape) circular cylinder, the shape of which is specified in the file ``circle_200points``, is located behind the square airfoil from the file ``square_160points`` (set by the ``diamond`` at an angle of attack of 45 degrees) in its vortex wake. All other parameters are taken from the default values.

It is also possible to load the previously saved distribution of vortex particles (parameter *fileWake*) and positions of stationary point sources/sinks in the flow domain (parameter *fileSource*). The directory where these files are located is specified by the *wakesDir* parameter. Sources and sinks are typically used to simulate internal flows.


In the considered example of the *task passport* two parameters are not defined explicitly: the value of the calculation time-step *dt* and the angle of incidence *angle*. The corresponding template parameters are marked with *$*, which means that these parameters should be taken from the values ​​defined by the user in the ``problems`` file for this task (in parenthesis after the task label). This allows you to make "universal" passport  and automate the procedure for solving problems: in order to use this opportunity, you must specify the *copyPath* key in the file ``problems`` by name of the directory where the "universal" passport is stored, for example like this:


.. code-block:: c
   :caption: problems
   :name: problems2

   problems = {
       wing00deg(np = 1, angle =  0, tau = 1.5e-2, copyPath = "./wingBase"),
       wing05deg(np = 2, angle =  5, tau = 1.5e-2, copyPath = "./wingBase"),
       wing10deg(np = 2, angle = 10, tau = 2.0e-2, copyPath = "./wingBase"),
   };
   
As a result, all the necessary directories corresponding to the labels of the tasks to be solved will be created automatically, and all the files contained in the folder specified in the *copyPath* key will be copied into them (at least there should be a file with *task passport*, which must contain template parameters, otherwise all calculations will be identical!).



Points of velocity and pressure calculation ``pointsVP``
--------------------------------------------------------

File ``pointsVP`` lists the points where during the simulation the velocities and pressure will be periodically calculated in order to save the corresponding fields to files.

In this file, you can specify two keys: 

*    *points* - points where velocities and pressure will be calculated and saved in vtk-files
*    *history* - points, for each of which a file will additionally be created containing the history of corresponding parameters at a given point.

File ``pointsVP`` may be absent if the parameter *saveVP* is set to zero in the simulation passport.
