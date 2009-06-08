mpcomplex.o : mpcomplex.cpp mpcomplex.h
	gcc -c -o mpcomplex.o mpcomplex.cpp

geounweld : geounweld.cpp mpcomplex.o
	gcc -o geounweld geounweld.cpp -L /usr/local/lib/*.a -L mpcomplex.o