#==============================================================================
#	VOID FINDER
#==============================================================================
# Sebastian Bustamante (Universidad de Antioquia), sebastian.bustamante@gmail.com


#Calculating density field from trace of eigenvalues (Tweb method)
OPT += -DDENSITY_TRACE


CC	=	gcc
CFLAGS	=	-g -I. -c $(OPT)

#Void Finder - Program
Void_Finder:inout.o finder.o tools.o void_finder.o
	gcc inout.o finder.o void_finder.o tools.o -lm -o Void_Finder.out
	rm -r *.o
	
#Density Profile of voids - Program
Density:inout.o finder.o tools.o density.o
	gcc inout.o finder.o density.o tools.o -lm -o Density_Profile.out
	rm -r *.o

edit:
	kate *.c *.h &

clean:
	rm -r *.o *.out *.png *.tmp