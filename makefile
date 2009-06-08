geounweld : geounweld.cpp mpcomplex.o
	g++ -o geounweld geounweld.cpp -L /usr/local/lib/*.a -L mpcomplex.o

mpcomplex.o : mpcomplex.cpp mpcomplex.h
	g++ -c -o mpcomplex.o mpcomplex.cpp
