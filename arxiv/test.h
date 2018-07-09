// Standard things to include
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

// These are the includes to use "Root" things 
#include "TInterpreter.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
// These are the larsoft includes that let you
// have access to data-products and the event 
// details
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"

// data-products themselves 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "uboone/EventWeight/MCEventWeight.h"

// Same as above but this lets you have access to 
/*
#include "/uboone/app/users/jaz8600/work/GALLERY/gallery-framework/UserDev/BasicTool/GeoAlgo/GeoAlgo.h"
#include "/uboone/app/users/jaz8600/work/GALLERY/gallery-framework/UserDev/BasicTool/GeoAlgo/GeoVector.h"
#include "/uboone/app/users/jaz8600/work/GALLERY/gallery-framework/UserDev/BasicTool/GeoAlgo/GeoHalfLine.h"
#include "/uboone/app/users/jaz8600/work/GALLERY/gallery-framework/UserDev/BasicTool/GeoAlgo/GeoTrajectory.h"
*/
//#include "GeoVector.h"
//#include "GeoAABox.h"
//#include "GeoHalfLine.h"
//#include "GeoAlgo.h"

// associations 
#include "canvas/Persistency/Common/FindMany.h"
