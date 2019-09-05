CC = g++
EXEC = Heff
LIBS = 
FLAGS = -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11


all: main.o
	$(CC) *.o -o $(EXEC) $(LIBS)

main.o : main.cpp
	$(CC) main.cpp -c $(FLAGS)


clear :
	rm -f *.o

mr_proper :
	rm -f *.o $(EXEC)
