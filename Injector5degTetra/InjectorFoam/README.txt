interPhaseChangeFoam case - uses water flowing through 5deg injector "wedge" under a pressure differential. No turbulence model applied. Mesh = 3D Tetrahedral created by SALOME NETGEN-1D-2D-3D.

To view sim results: open injector.foam on ParaView

*** NOTE ***
to run the sim, you will have to delete all numbered folders (e.g. '0.01' --- EXCEPT 0 ), and all folders with prefix 'processor'

decomposeParDict - currently configured to run with 4 cores, uses scotch decomposition (OpenFOAM automatic algorithm)

to run parallel: cd to working directory, run:
	decomposePar
then run:
	mpirun -np 4 interPhaseChangeFoam -parallel
once complete, run:
	reconstructPar 
to reconstruct parallelised fields. 
to see new results type:	
	touch <FILENAME>.foam

**optional**
mpirun -np 4 interPhaseChangeFoam -parallel > log.interPhaseChangeFoam &  (this reduces the difference between ExecutionTime and ClockTime)
tail -f log.interPhaseChangeFoam (this allows you to follow it as usual)