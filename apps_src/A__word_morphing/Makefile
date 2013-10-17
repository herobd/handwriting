.PHONY: clean

../../bin/word_morphing: word_morphing.cpp
	g++ -Wall -march=native -O3 -g -rdynamic -fPIC word_morphing.cpp -I../../src -o ../../bin/word_morphing -L../../lib/ -ldocumentproj_2013.08.30 -ljpeg -ltiff -lpng -pthread -lm


clean:
	@- rm ../../bin/word_morphing
