#ifndef TPCOBJECTMAKER_H
#define TPCOBJECTMAKER_H

#include "AnaHelper.h"
#include "RecoTrueHelper.h"
#include "TpcObjectHelper.h"
#include "TPCObject.h"
#include "MCGhost.h"

namespace xsec_ana {
class NueXSecProducer;
}

class xsec_ana::NueXSecProducer : public art::EDProducer {
public:
explicit NueXSecProducer(fhicl::ParameterSet const & p);

// Plugins should not be copied or assigned.
NueXSecProducer(NueXSecProducer const &) = delete;
NueXSecProducer(NueXSecProducer &&) = delete;
NueXSecProducer & operator = (NueXSecProducer const &) = delete;
NueXSecProducer & operator = (NueXSecProducer &&) = delete;

void reconfigure(fhicl::ParameterSet const &p) override;
void produce(art::Event & e) override;

private:

std::string _pfp_producer;
std::string _hitfinderLabel;
std::string _geantModuleLabel;
std::string _spacepointLabel;
std::string _neutrino_flash_match_producer;
std::string _cosmic_flash_match_producer;
std::string _opflash_producer_beam;
std::string _potsum_producer;
std::string _potsum_instance;
std::string _particle_id_producer;

std::string _vertexLabel;
std::string _trackLabel;
std::string _showerLabel;

bool _use_genie_info;
int _minimumHitRequirement;

double _beam_spill_start;
double _beam_spill_end;

bool _debug;
bool _verbose;
bool isMC;
bool isData;

int run;
int event;

};

#endif
