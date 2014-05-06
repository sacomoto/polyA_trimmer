all:kseq.h polyA_trimmer.cpp
	g++ -O3 polyA_trimmer.cpp -o polyA_trimmer -lz

clean:
	rm -f *.o polyA_trimmer