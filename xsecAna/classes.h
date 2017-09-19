#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

//#include "uboone/NueXSecModules/xsecAna/"
#include "uboone/NueXSecModules/xsecAna/TPCObject.h"
#include "uboone/NueXSecModules/xsecAna/MCGhost.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>
/*
template class art::Assns<anab::FlashMatch,recob::PFParticle>;
template class art::Assns<recob::PFParticle,anab::FlashMatch>;

template class art::Wrapper<art::Assns<anab::FlashMatch,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, anab::FlashMatch >>;


template class std::vector<ubana::FlashMatch>;

template class art::Assns<ubana::FlashMatch,recob::PFParticle>;
template class art::Assns<recob::PFParticle,ubana::FlashMatch>;
template class art::Assns<ubana::FlashMatch,recob::Track,void>;
template class art::Assns<recob::Track,ubana::FlashMatch,void>;

template class art::Wrapper<art::Assns<ubana::FlashMatch,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, ubana::FlashMatch >>;
template class art::Wrapper<art::Assns<ubana::FlashMatch,recob::Track,void> >;
template class art::Wrapper<art::Assns<recob::Track,ubana::FlashMatch,void> >;
*/



template class std::vector<xsecAna::TPCObject>;

template class art::Assns<xsecAna::TPCObject,recob::PFParticle>;
template class art::Assns<recob::PFParticle,xsecAna::TPCObject>;
template class art::Assns<xsecAna::TPCObject,recob::Track,void>;
template class art::Assns<recob::Track,xsecAna::TPCObject,void>;
template class art::Assns<xsecAna::TPCObject,recob::Shower,void>;
template class art::Assns<recob::Shower,xsecAna::TPCObject,void>;
//template class art::Assns<ubana::TPCObject,ubana::FlashMatch,void>;
//template class art::Assns<ubana::FlashMatch,ubana::TPCObject,void>;

template class art::Wrapper<art::Assns<xsecAna::TPCObject,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, xsecAna::TPCObject >>;
template class art::Wrapper<art::Assns<recob::Track,xsecAna::TPCObject,void> >;
template class art::Wrapper<art::Assns<xsecAna::TPCObject,recob::Track,void> >;
template class art::Wrapper<art::Assns<recob::Shower,xsecAna::TPCObject,void> >;
template class art::Wrapper<art::Assns<xsecAna::TPCObject,recob::Shower,void> >;
//template class art::Wrapper<art::Assns<ubana::TPCObject,ubana::FlashMatch,void> >;
//template class art::Wrapper<art::Assns<ubana::FlashMatch,ubana::TPCObject,void> >;



template class std::vector<xsecAna::MCGhost>;

template class art::Assns<simb::MCParticle,xsecAna::MCGhost>;
template class art::Assns<xsecAna::MCGhost,simb::MCParticle>;
template class art::Assns<recob::PFParticle,xsecAna::MCGhost>;
template class art::Assns<xsecAna::MCGhost,recob::PFParticle>;
template class art::Wrapper<art::Assns<simb::MCParticle,xsecAna::MCGhost>>;
template class art::Wrapper<art::Assns<xsecAna::MCGhost,simb::MCParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle,xsecAna::MCGhost>>;
template class art::Wrapper<art::Assns<xsecAna::MCGhost,recob::PFParticle>>;
