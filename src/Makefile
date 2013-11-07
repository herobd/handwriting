#make has predefined values for $(CC) $(CPP) $(CXX) $(LD) and $(RM).
#CC is the C compiler, CPP is the C/C++ preprocessor, CXX is C++ compiler,
#LD is the linker, and RM is the remove command

#directories for headers, separated by spaces with -I in front of each
INC+=-I.
#where to put the library file when it is created
LIBDIR=../lib
#where to dump object files
OBJDIR=../obj
#C compiler flags (use += so user can pass stuff in as options to make)
CFLAGS+= -Wall -O0 -g -fPIC -pthread #Brian edit, used to be O3
#C/C++ PreProcessor flags (use += so we don't clobber what user passes to make)
CPPFLAGS+=
#C++ compiler flags (use += so we don't clobber what user passes to make)
#CXXFLAGS+= -Wall -march=native -O3 -pg -rdynamic -fPIC
CXXFLAGS+= -Wall -march=native -O0 -g -fPIC -pthread #Brian edit, used to be O3
#linker flags (use += so we don't clobber what user passes to make)
LDFLAGS+=
#name of the library to create
TARGET=$(LIBDIR)/libdocumentproj_2013.08.30.a

#the wildcard does globbing, sets CXX_SRCS to all files in this dir that match
CXX_SRCS= $(wildcard *.cpp)
#replace the .cpp suffix with .o
CXX_OBJS_TMP:= ${CXX_SRCS:.cpp=.o}
#now substitute the ./ with the OBJDIR
#CXX_OBJS:= $(subst ./,$(OBJDIR),$(CXX_OBJS_TMP))#this and next line equivalent
CXX_OBJS:= $(addprefix  $(OBJDIR)/,$(CXX_OBJS_TMP))

#header dependency lists for the sources
#I was auto-generating them, but now they are just pasted at bottom of Makefile
#DEPS= ${CXX_OBJS:.o=.depends}



#the .PHONY means there is no file associated with these targets
.PHONY: clean all distclean

all: $(TARGET)

$(TARGET): $(CXX_OBJS)
	ar -rs $(TARGET) $(CXX_OBJS)

# the $< represents the source file and the $@ represents the output file
# the -MM flag to gcc generates dependencies
#$(OBJDIR)/%.depends: $(notdir %.cpp)
#	@ echo -n "$(OBJDIR)/" > $@
#	$(CXX) -MM $(CXXFLAGS) $< >> $@

$(OBJDIR)/%.o: $(notdir %.cpp)
	$(CXX) -c $(CXXFLAGS) $< -o $@

# @ means don't print the command (run it silently) and - means to ignore errors
clean:
	@- $(RM) $(CXX_OBJS) $(TARGET)

distclean: clean

#include the list of header dependencies for each source file
#-include $(DEPS)


#########################################################################
# dependencies are listed below.  To generate them, use g++ -MM
# and then paste them below
#########################################################################


../obj/dbackgroundremover.o: dbackgroundremover.cpp dbackgroundremover.h \
 dimage.h ddefs.h dinttypes.h dsize.h dmedianfilter.h dmaxfilter.h \
 dkernel2d.h dconvolver.h dedgedetector.h dtimer.h dvariancefilter.h

../obj/dbilateralfilter.o: dbilateralfilter.cpp dbilateralfilter.h dimage.h \
 ddefs.h dinttypes.h dsize.h dprogress.h dinstancecounter.h \
 paris_sylvain_bilateral_headers/fast_color_bf.h \
 paris_sylvain_bilateral_headers/mixed_vector.h \
 paris_sylvain_bilateral_headers/array.h \
 paris_sylvain_bilateral_headers/array_n.h \
 paris_sylvain_bilateral_headers/geom.h \
 paris_sylvain_bilateral_headers/math_tools.h \
 paris_sylvain_bilateral_headers/msg_stream.h \
 paris_sylvain_bilateral_headers/../dprogress.h \
 paris_sylvain_bilateral_headers/fast_lbf.h dtimer.h

../obj/dcolorspace.o: dcolorspace.cpp dcolorspace.h dimage.h ddefs.h dinttypes.h \
 dsize.h

../obj/dconnectedcompdistance.o: dconnectedcompdistance.cpp \
 dconnectedcompdistance.h dimage.h ddefs.h dinttypes.h dsize.h

../obj/dconnectedcompinfo.o: dconnectedcompinfo.cpp dconnectedcompinfo.h ddefs.h \
 dinttypes.h dinstancecounter.h

../obj/dconnectedcomplabeler.o: dconnectedcomplabeler.cpp \
 dconnectedcomplabeler.h dimage.h ddefs.h dinttypes.h dsize.h \
 dconnectedcompinfo.h

../obj/dconvolver.o: dconvolver.cpp dconvolver.h dimage.h ddefs.h dinttypes.h \
 dsize.h dkernel2d.h dprogress.h

../obj/ddistancemap.o: ddistancemap.cpp ddistancemap.h dimage.h ddefs.h \
 dinttypes.h dsize.h

../obj/ddynamicprogramming.o: ddynamicprogramming.cpp ddynamicprogramming.h \
 dfeaturevector.h ddefs.h dinttypes.h dinstancecounter.h dimage.h dsize.h \
 dmath.h dcolorspace.h

../obj/dedgedetector.o: dedgedetector.cpp dedgedetector.h dimage.h ddefs.h \
 dinttypes.h dsize.h dkernel2d.h dconvolver.h

../obj/dfeaturevector.o: dfeaturevector.cpp dfeaturevector.h ddefs.h dinttypes.h \
 dinstancecounter.h

../obj/dglobalskew.o: dglobalskew.cpp dglobalskew.h dimage.h ddefs.h dinttypes.h \
 dsize.h dprofile.h dedgedetector.h dthresholder.h dmath.h

../obj/dhistogram.o: dhistogram.cpp

../obj/dhough.o: dhough.cpp dhough.h dimage.h ddefs.h dinttypes.h dsize.h \
 dmath.h dprogress.h dinstancecounter.h dedgedetector.h dconvolver.h \
 dkernel2d.h

../obj/dimage.o: dimage.cpp dmemalign.h dimage.h ddefs.h dinttypes.h dsize.h \
 dinstancecounter.h dimageio.h drect.h dpoint.h

../obj/dimageio.o: dimageio.cpp dimageio.h ddefs.h dinttypes.h dimage.h dsize.h

../obj/dinstancecounter.o: dinstancecounter.cpp dinstancecounter.h

../obj/dkernel2d.o: dkernel2d.cpp dkernel2d.h dimage.h ddefs.h dinttypes.h \
 dsize.h dmemalign.h dinstancecounter.h

../obj/dmaxfilter.o: dmaxfilter.cpp dmaxfilter.h dimage.h ddefs.h dinttypes.h \
 dsize.h dprogress.h dinstancecounter.h dthreads.h

../obj/dmeanfilter.o: dmeanfilter.cpp dmeanfilter.h dimage.h ddefs.h dinttypes.h \
 dsize.h dprogress.h dinstancecounter.h

../obj/dmedialaxis.o: dmedialaxis.cpp dmedialaxis.h dimage.h ddefs.h dinttypes.h \
 dsize.h ddistancemap.h dthresholder.h

../obj/dmedianfilter.o: dmedianfilter.cpp dmedianfilter.h dimage.h ddefs.h \
 dinttypes.h dsize.h dprogress.h dinstancecounter.h dthreads.h

../obj/dminfilter.o: dminfilter.cpp dminfilter.h dimage.h ddefs.h dinttypes.h \
 dsize.h dprogress.h dinstancecounter.h dthreads.h

../obj/dmorphink.o: dmorphink.cpp dmorphink.h dimage.h ddefs.h dinttypes.h \
 dsize.h dmath.h dprofile.h dfeaturevector.h dinstancecounter.h \
 ddynamicprogramming.h ddistancemap.h dmedialaxis.h dtimer.h \
 dwordfeatures.h

../obj/dmorphology.o: dmorphology.cpp dmorphology.h dimage.h ddefs.h dinttypes.h \
 dsize.h dinstancecounter.h

../obj/dpoint.o: dpoint.cpp dpoint.h dinstancecounter.h

../obj/dprofile.o: dprofile.cpp dprofile.h dimage.h ddefs.h dinttypes.h dsize.h \
 dmath.h dinstancecounter.h

../obj/dprogress.o: dprogress.cpp dprogress.h dinstancecounter.h

../obj/drect.o: drect.cpp drect.h dpoint.h dinstancecounter.h

../obj/dsize.o: dsize.cpp dsize.h dinstancecounter.h

../obj/dslantangle.o: dslantangle.cpp dslantangle.h dimage.h ddefs.h dinttypes.h \
 dsize.h drect.h dpoint.h dmath.h dprofile.h dthreads.h

../obj/dthresholder.o: dthresholder.cpp dthresholder.h dimage.h ddefs.h \
 dinttypes.h dsize.h dprogress.h dkernel2d.h dconvolver.h \
 dconnectedcomplabeler.h dtimer.h dthreads.h

../obj/dtimer.o: dtimer.cpp dtimer.h dinstancecounter.h

../obj/dvariancefilter.o: dvariancefilter.cpp dvariancefilter.h dimage.h ddefs.h \
 dinttypes.h dsize.h

../obj/dwordfeatures.o: dwordfeatures.cpp dwordfeatures.h dimage.h ddefs.h \
 dinttypes.h dsize.h dfeaturevector.h dinstancecounter.h



