CXX=g++
CXXFLAGS=$(shell root-config --cflags)
LIBS=$(shell root-config --libs)
LD=g++
LDFLAGS=
DEBUGFLAGS= -O0 -D _DEBUG

XSECANALIBS=../xsecAna/libxsecAna.so ../xsecAna/libxsecAna_dict.so

TARGETNAMES=out_inspect selection potSum
FUNCTIONCLASSOBJECT=selection_functions.o

TARGETS=$(TARGETNAMES:%=%.exe)
OBJECTS=$(TARGETNAMES:%=%.o)
SOURCES=$(TARGETNAMES:%=%.cxx)

#$(info $(OBJECTS))

all: $(TARGETS)

clean:
	rm -f $(OBJECTS) $(TARGETS)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c -o "$@" $^

%.exe : %.o $(XSECANALIBS) $(FUNCTIONCLASSOBJECT)
	$(LD) $(LDFLAGS) $(LIBS) -o "$@" $^
