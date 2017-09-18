#include "TPCObject.h"
#include <vector>
#include <iostream>
#include <stdexcept>

namespace xsecAna {

TPCObject::TPCObject() {
	fOrigin      = simb::kUnknown;
	fOriginExtra = simb::kUnknown;//was kUnSet before
}

// TPCObject::~TPCObject(){
// }

// Setter methods
void TPCObject::SetTracks      (std::vector<art::Ptr<recob::Track> > tracks)    {
	this->fTracks = tracks;
}
void TPCObject::SetShowers     (std::vector<art::Ptr<recob::Shower> > showers)  {
	this->fShowers = showers;
}
void TPCObject::SetPFPs        (std::vector<recob::PFParticle> pfps) {
	this->fPFParticles = pfps;
}
void TPCObject::SetVertex      (recob::Vertex vertex)                {
	this->fVertex = vertex;
}
void TPCObject::SetParticleOrigins (std::vector< simb::Origin_t > origins)
{
	this->fOriginVector = origins;
}
void TPCObject::SetOrigin      (simb::Origin_t origin)       {
	this->fOrigin = origin;
}
void TPCObject::SetOriginExtra (simb::Origin_t origin)  {
	this->fOriginExtra = origin;
}
void TPCObject::SetMultiplicity(int p, int t, int s)                 {
	this->fPfpMult = p; this->fTrackMult = t; this->fShowerMult = s;
}

// Getter methods
const std::vector<art::Ptr<recob::Track> >       & TPCObject::GetTracks()      const {
	return this->fTracks;
}
const std::vector<art::Ptr<recob::Shower> >      & TPCObject::GetShowers()     const {
	return this->fShowers;
}
const std::vector<recob::PFParticle>  & TPCObject::GetPFPs()        const {
	return this->fPFParticles;
}
const recob::Vertex                   & TPCObject::GetVertex()      const {
	return this->fVertex;
}
const std::vector< simb::Origin_t >   & TPCObject::GetParticleOrigins()  const {
	return this->fOriginVector;
}
const simb::Origin_t                  & TPCObject::GetOrigin()      const {
	return this->fOrigin;
}
const simb::Origin_t                  & TPCObject::GetOriginExtra() const {
	return this->fOriginExtra;
}


const size_t TPCObject::GetNTracks()  const {
	return (this->fTracks).size();
}
const size_t TPCObject::GetNPFP()     const {
	return (this->fPFParticles).size();
}
const size_t TPCObject::GetNShowers() const {
	return (this->fShowers).size();
}

const void TPCObject::GetMultiplicity(int &p, int &t, int &s) const {
	p = this->fPfpMult; t = this->fTrackMult; s = this->fShowerMult;
}

const int TPCObject::GetNTracksCloseToVertex(double tolerance) const {

	if (tolerance < 0.) {
		std::cerr << "[TPCObject] tolerance can't be less than zero" << std::endl;
		throw std::exception();
	}

	// Output
	int multiplicity = 0;

	// Vertex
	double v[3];
	this->fVertex.XYZ(v);
	TVector3 vtx(v[0], v[1], v[2]);

	// Loop over tracks and calulate multiplicity
	for (auto t : this->fTracks) {

		TVector3 start = t->Vertex();
		TVector3 end = t->End();

		if ( (vtx-start).Mag() < tolerance   ||   (vtx-end).Mag() < tolerance ) {
			multiplicity++;
		}
	}

	return multiplicity;
}

} //end namespace
