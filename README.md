# Elias Lopez @2016 msn.eliaslopez@gmail.com
carpeta para testear codigo mpi en python para correr autodock
importante:
	*	tener instalado autodock4, mpi4py, mpirun, etc. 
	*	la funcion trivialmente paralelizable esta embedida en la funcion f(i)
	*	La funcion esta implementada para python 2.7 NO esta probada en python3
	*	La funcion implementada supone que las carpetas para docking contienen todos
		los archivo y el particular el archivo dock.dpf.  La salida en cada carpeta 
		la llama "nombrecarpeta"_dock_out.pdb


MODO DE USO:

	mpirun -np 4  python mpi4py_trivial_autodock.py -f 1bmk* 
	donde
		np: Numero de procesadores.
		1bmk* son todas carpetas replicas del mismo sistema (runs de docking = 2).
