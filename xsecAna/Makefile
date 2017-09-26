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

TARGETNAME=$(shell basename "`pwd`")
DICTNAME=$(TARGETNAME)_dict

TARGETS=lib$(TARGETNAME).so lib$(TARGETNAME)_dict.so

all: $(TARGETS)

clean:
	rm -f $(OBJECTS) $(TARGETS)

debug:
	@echo "Target name: '$(TARGETNAME)'"

$(TARGET): $(OBJECTS) $(HEADERS) LinkDef.h
	$(CXX) $(CXXFLAGS) $(LIBS) $(DEBUGFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

$(DICTNAME).cxx: $(HEADERS) LinkDef.h
	rootcling -f $@ $^

%.o: %.cxx
	$(CXX) $(COMPILE_FLAGS) $(CXXFLAGS) $(OPTFLAGS) -c -o "$@" $^

lib$(TARGETNAME).so: $(OBJECTS) 
	$(LD) $(LDFLAGS) $(LIBS) -o "$@" $^

lib$(TARGETNAME)_dict.so: lib$(TARGETNAME).so $(DICTNAME).o 
	$(LD) $(LDFLAGS) $(LIBS) -o "$@" $^

