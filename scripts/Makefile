XSECANADIR=../xsecAna
CXX=g++
CXXFLAGS=$(shell root-config --cflags)
CPPFLAGS=-I$(XSECANADIR)
LIBS=$(shell root-config --libs) -L$(XSECANADIR) -luboonecode_uboone_NueXSecModules # -lxsecAna_dict
LD=g++
LDFLAGS=
DEBUGFLAGS= -O0 -D _DEBUG

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
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUGFLAGS) -c -o "$@" $^

%.exe : %.o $(FUNCTIONCLASSOBJECT)
	$(LD) $(LDFLAGS) $(LIBS) -o "$@" $^
