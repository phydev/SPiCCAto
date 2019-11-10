# Makefile

# compilers/flags
# compiler = ifort 
# debug = -g -traceback
# flags = -O3 -fast -r8 $(debug) 

compiler = gfortran
tests = -fprofile-arcs -ftest-coverage
flags = -03 -fdefault-real-8 $(tests)

# the program

#dependencies
run-SPiCCAto: global.o routines.o class_mesh.o main.o
	$(compiler) $(flags) global.o routines.o class_mesh.o main.o -o run-SPiCCAto 
global.mod: global.o global.F90
	$(compiler) $(flags) -c global.F90
global.o: global.F90
	$(compiler) $(flags) -c global.F90
routines.mod: global.mod routines.o routines.F90 
	$(compiler) $(flags) -c global.mod routines.F90
routines.o: global.mod routines.F90
	$(compiler) $(flags) -c global.mod routines.F90
class_mesh.mod: global.mod routines.mod class_mesh.F90
	$(compiler) $(flags) -c global.mod routines.mod class_mesh.F90
class_mesh.o: global.mod routines.mod class_mesh.F90
	$(compiler) $(flags) -c global.mod routines.mod class_mesh.F90
main.o: global.mod routines.mod class_mesh.mod main.F90
	$(compiler) $(flags) -c main.F90
clean:
	rm -f  run-SPiCCAto global.mod global.o routines.mod routines.o class_mesh.mod class_mesh.o main.o



