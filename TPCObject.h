/**
 * \class ubana::TPCObject
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a TPC Object
 *
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/03/02 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, May 15, 2017 at 16:34:34
 *
 */

#ifndef TPCObject_h
#define TPCObject_h

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include <vector>

namespace xsec_ana {

enum TPCObjectOrigin {
	kUnknown = -1,      // -1
	kBeamNeutrino = 0,  // 0
	kCosmicRay,         // 1
	kMixed,             // 2
};

enum TPCObjectOriginExtra {
	kNotSet = -1,       // -1
	kStoppingMuon = 0,  // 0
	kACPT,              // 1
	kNCPion,            // 2
	kNCProton,          // 3
};
}


namespace xsec_ana {

class TPCObject {

public:

TPCObject();
virtual ~TPCObject();

// Setter methods
void SetTracks(std::vector<recob::Track>);
void SetShowers(std::vector<recob::Shower>);
void SetPFPs(std::vector<recob::PFParticle>);
void SetVertex(recob::Vertex);
void SetOrigin(xsec_ana::TPCObjectOrigin);
void SetOriginExtra(xsec_ana::TPCObjectOriginExtra);
void SetParticleOrigins (std::vector< x_sec_ana::TPCObjectOrigin> origins)
void SetMultiplicity(int pfpMult, int trackMult, int showerMult);

// Getter methods
const std::vector<recob::Track>               & GetTracks()          const;
const std::vector<recob::Shower>              & GetShowers()         const;
const std::vector<recob::PFParticle>          & GetPFPs()            const;
const recob::Vertex                           & GetVertex()          const;
const xsec_ana::TPCObjectOrigin               & GetOrigin()          const;
const xsec_ana::TPCObjectOriginExtra          & GetOriginExtra()     const;
const std::vector<xsec_ana::TPCObjectOrigin>  & GetParticleOrigins() const;
const size_t                                  GetNTracks()           const;
const size_t                                  GetNShowers()          const;
const size_t                                  GetNPFP()              const;
const void                                    GetMultiplicity(int &, int &, int &) const;
const int                                     GetNTracksCloseToVertex(double)      const;

private:

std::vector<recob::Track>      fTracks;
std::vector<recob::Shower>     fShowers;
std::vector<recob::PFParticle> fPFParticles;
recob::Vertex fVertex;
std::vector< xsec_ana::TPCObjectOrigin > fOriginVector;
xsec_ana::TPCObjectOrigin fOrigin;
xsec_ana::TPCObjectOriginExtra fOriginExtra;
int fPfpMult;
int fTrackMult;
int fShowerMult;
};
} //end namespace

#endif /* TPCObject_h */
