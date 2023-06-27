VM2D run and processing the results
===================================

VM2D run
--------
		  
The program is launched from the *working directory*, i.e. a directory containing the files ``problems``, ``defaults``, ``switchers`` and ``mechanics``, as well as folders for specific tasks to be solved. To start simulation, you need to call the command

      ``VM2D``

(if necessary, specifying the full path to the compiled executable file; this path can be specified in the system variable ``PATH``)


Saving results
--------------

Vortex wakes are saved in vtk files in the ``snapshots`` subdirectory, the velocity and pressure fields calculated at the given points --- in the ``velPress`` subdirectory. Both of these folders are located in the task subdirectory and are created automatically.

Also, information about the loads acting on the airfoil or on each airfoil from the system of airfoils is stored in the task subdirectory. This information is saved in files named ``forces-airfoil-n'', where ``n'' corresponds to the airfoil number in the list from the *airfoil* key in the *task passport* file. These files are saved in two formats: text (without extension) and *csv*. 

If the airfoil is movable, then information about its position is saved to the file ``position-airfoil-n`` also in two formats.

To view vtk and csv files, *paraview* free open source postprocessor is useful.
