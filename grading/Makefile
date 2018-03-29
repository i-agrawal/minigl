CXX = g++
CXXFLAGS = -Wall -g -O3 -std=c++11 -w
CXXFLAGS2 = -Wall -g -O3 -std=c++11 -DNO_OPENGL -w
LIBS  = -lGL -lGLU -lglut -lX11 -lpng
LIBS2  = -lpng
SRCS = $(wildcard *.cpp)
OBJS := ${SRCS:.cpp=.o}
OBJS2 := main-nogl.o minigl.o parse-nogl.o dump_png.o

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBS = -framework OpenGL -framework GLUT -lGLEW -lpng
	CXXFLAGS += -Wno-deprecated-declarations
endif

TARGET = minigl
TARGET2 = minigl-nogl

.PHONY: all clean

all: $(TARGET) $(TARGET2)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(LIBS)
$(TARGET2): $(OBJS) 
	$(CXX) -c $(CXXFLAGS2) $(CPPFLAGS) -o main-nogl.o main.cpp
	$(CXX) -c $(CXXFLAGS2) $(CPPFLAGS) -o parse-nogl.o parse.cpp
	$(CXX) $(CXXFLAGS2) $(CPPFLAGS) -o $(TARGET2) $(OBJS2) $(LIBS2)

clean:
	rm -f *.o *.d *~ *.bak $(TARGET) $(TARGET2)

# Generate a list of dependencies for each cpp file
%.d: $(SRCS)
	@ $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $*.cpp | sed -e 's@^\(.*\)\.o:@\1.d \1.o:@' > $@
	@ $(CXX) $(CPPFLAGS) $(CXXFLAGS2) -MM $*.cpp | sed -e 's@^\(.*\)\.o:@\1-nogl.d \1.o:@' > $@

# At the end of the makefile
# Include the list of dependancies generated for each object file
# unless make was called with target clean
ifneq "$(MAKECMDGOALS)" "clean"
-include ${SRCS:.cpp=.d}
endif
