### Makefile du programme

# Compilateur
CF := gfortran

# Compilation
all : exe

exe : main.f90 tools.o functions.o method.o
	@echo 'Compilation de main.f90 ...'
	$(CF) -o exe main.f90 tools.o functions.o method.o

functions.o : functions.f90 tools.o
	@echo 'Compilation de functions.f90 ...'
	$(CF) -o functions.o -c functions.f90

method.o : method.f90 tools.o
	@echo 'Compilation de method.f90 ...'
	$(CF) -o method.o -c method.f90

tools.o : tools.f90
	@echo 'Compilation de tools.f90 ...'
	$(CF) -o tools.o -c tools.f90

# Commandes supplémentaires
clean :
	rm *.o
	rm *.mod
