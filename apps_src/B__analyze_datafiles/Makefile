.PHONY: clean

../../bin/analyze_datafiles: analyze_datafiles.cpp
	g++ -Wall -march=native -O3 -g -fPIC analyze_datafiles.cpp -o ../../bin/analyze_datafiles

clean:
	@- rm ../../bin/analyze_datafiles
