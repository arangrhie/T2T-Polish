
all: samToAlignment samToErrorRate

samToAlignment: bamcat.o samToAlignment.o
	g++ -g -o samToAlignment -lz -lpthread bamcat.o samToAlignment.o

samToErrorRate: bamcat.o samToErrorRate.o
	g++ -g -o samToErrorRate -lz -lpthread bamcat.o samToErrorRate.o

bamcat.o: bamcat.c
	gcc -g -c bamcat.c

samToAlignment.o: samToAlignment.C
	g++ -g -c samToAlignment.C

samToErrorRate.o: samToErrorRate.C
	g++ -g -c samToErrorRate.C

objclean:
	rm -f bamcat.o
	rm -f samToAlignment.o
	rm -f samToErrorRate.o

clean: objclean
	rm -f samToAlignment
	rm -f samToErrorRate

.PHONY: all objclean clean

