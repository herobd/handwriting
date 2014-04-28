.PHONY: clean

../../bin/word_morphing: word_morphing.cpp
	g++ -Wall -march=native -O0 -g -rdynamic -fPIC word_morphing.cpp -I../../src -o ../../bin/word_morphing -L../../lib/ -ldocumentproj_2013.08.30 -ljpeg -ltiff -lpng -pthread -lm #Brian edit, should be O3

test: test.cpp ../../obj/dmorphink.o
	g++ -Wall -march=native -O0 -g -rdynamic -fPIC test.cpp -I../../src -o ../../bin/test_word_morphing -L../../lib/ -ldocumentproj_2013.08.30 -ljpeg -ltiff -lpng -pthread -lm

testdif: ../../bin/test_word_morphing
	gdb --args ../../bin/test_word_morphing 3 7

testsame: ../../bin/test_word_morphing
	gdb --args ../../bin/test_word_morphing 14 21

testhard: ../../bin/test_word_morphing
	gdb --args ../../bin/test_word_morphing 1421 774

clean:
	@- rm ../../bin/word_morphing
	@- rm ../../bin/test_word_morphing
