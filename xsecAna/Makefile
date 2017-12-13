CXX=g++
CXXFLAGS=$(shell root-config --cflags)
LIBS=$(shell root-config --libs)
LD=g++
LDFLAGS= -shared
DEBUGFLAGS= -O0 -D _DEBUG
COMPILE_FLAGS+=-fPIC

SOURCENAMES=TpcObjectContainer ParticleContainer
HEADERS=$(SOURCENAMES:%=%.h)
SOURCES=$(SOURCENAMES:%=%.cxx)
OBJECTS=$(SOURCENAMES:%=%.o)

DICTSOURCE=$(SOURCES)

# TARGETNAME=$(shell basename "`pwd`")
TARGETNAME=uboonecode_uboone_NueXSecModules
DICTNAME=$(TARGETNAME)_dict
DICTOBJECTS=$(DICTNAME).o $(DICTNAME).cxx $(DICTNAME)_rdict.pcm
TARGETS=lib$(TARGETNAME).so lib$(TARGETNAME)_dict.so

all: $(TARGETS)

clean:
	rm -f $(OBJECTS) $(TARGETS) $(DICTOBJECTS)

debug:
	@echo "Target name: '$(TARGETNAME)'"

$(TARGET): $(OBJECTS) $(HEADERS) LinkDef.h
	$(CXX) $(CXXFLAGS) $(LIBS) $(DEBUGFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

$(DICTNAME).cxx: $(HEADERS) LinkDef.h
	#rootcling -f $@ -rml lib$(TARGETNAME)_dict.so -rmf lib$(TARGETNAME)_dict.rootmap $^
	rootcling -f $@ $^

%.o: %.cxx
	$(CXX) $(COMPILE_FLAGS) $(CXXFLAGS) $(OPTFLAGS) -c -o "$@" $^

lib$(TARGETNAME).so: $(OBJECTS) 
	$(LD) $(LDFLAGS) $(LIBS) -o "$@" $^

lib$(TARGETNAME)_dict.so: lib$(TARGETNAME).so $(DICTNAME).o 
	$(LD) $(LDFLAGS) $(LIBS) -o "$@" $^

