#C++ compiler flags (use += so we don't clobber what user passes to make)
CXXFLAGS+= -Wall  -march=native -O3 -g -fPIC
#linker flags (use += so we don't clobber what user passes to make)
LDFLAGS+= -L../../lib
LFLAGS+= -ldocumentproj_2013.08.30 -ljpeg -ltiff -lpng -lm -pthread
#directories for headers, separated by spaces with -I in front of each
INC+=-I. -I../../src
BINPATH = ../../bin


.PHONY: clean all

all: $(BINPATH)/word_clustering $(BINPATH)/NxNtrainMatrixChunk $(BINPATH)/combineNxNChunks

$(BINPATH)/word_clustering: word_clustering.cpp
	g++ word_clustering.cpp -o $(BINPATH)/word_clustering $(CXXFLAGS) $(LDFLAGS) $(LFLAGS) $(INC)

$(BINPATH)/NxNtrainMatrixChunk: NxNtrainMatrixChunk.cpp
	g++ NxNtrainMatrixChunk.cpp -o $(BINPATH)/NxNtrainMatrixChunk $(CXXFLAGS) $(LDFLAGS) $(LFLAGS) $(INC)

$(BINPATH)/combineNxNChunks: combineNxNChunks.cpp
	g++ combineNxNChunks.cpp -o $(BINPATH)/combineNxNChunks $(CXXFLAGS) $(LDFLAGS) $(LFLAGS) $(INC)


clean:
	rm -f $(BINPATH)/word_clustering $(BINPATH)/NxNtrainMatrixChunk $(BINPATH)/combineNxNChunks

