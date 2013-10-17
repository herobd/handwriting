.PHONY: clean all

all:	chinese_word_morphing extract_gnt_image

chinese_word_morphing: chinese_word_morphing.cpp
	g++ -Wall -O3 -g -fPIC chinese_word_morphing.cpp -I../../src -o chinese_word_morphing -L../../lib/ -ldocumentproj_2013.08.30 -ljpeg -ltiff -lpng -pthread -lm

extract_gnt_image: extract_gnt_image.cpp
	g++ -Wall -O3 -g -fPIC extract_gnt_image.cpp -I../../src -o extract_gnt_image -L../../lib/ -ldocumentproj_2013.08.30 -ljpeg -ltiff -lpng -pthread -lm

clean:
	@- rm chinese_word_morphing extract_gnt_image
