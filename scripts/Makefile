XSECANADIR=../xsecAna
CXX=g++
CXXFLAGS=$(shell root-config --cflags)
CPPFLAGS=-I$(XSECANADIR)
LIBS=$(shell root-config --libs) -L$(XSECANADIR) -luboonecode_uboone_NueXSecModules # -lxsecAna_dict
LD=g++
LDFLAGS=
DEBUGFLAGS= -O0 -g -D _DEBUG

#TARGETNAMES=out_inspect selection potSum selection_slim selection_dataplots selection_datamc
TARGETNAMES=custom_potSum main
FUNCTIONCLASSOBJECT=selection_functions.o selection_functions_data.o selection_cuts.o histogram_functions.o utility.o selection.o selection_slim.o variation_output.o variation_output_bkg.o mc_truth.o

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
