OBJS = GISradar.o

make: $(OBJS)
	gcc -o GISradartest.exe GISradartest.c GISradar.o -lm
GISradar.o: 
	gcc -c GISradar.c -lm
	
run:
	./GISradartest.exe
	
clean:
	rm *.exe *.o
