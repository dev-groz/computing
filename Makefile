output: task1.o lagrange.o newton.o printing.o
	gcc task1.o lagrange.o newton.o printing.o -o task1
 
task1.o: task1.c
	gcc -c task1.c

lagrange.o: lagrange.c
	gcc -c lagrange.c
	
newton.o: newton.c
	gcc -c newton.c

printing.o: printing.c
	gcc -c printing.c

clean:
	del *.o
	del task1.exe