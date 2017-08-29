#include "AnaHelper.h"
#include "PostMCCorrecitons.h"

//_________________________________________________________________________________
double PostMCCorrecitons::GetCorrectedPhi(recob::Track t, recob::Vertex tpcobj_nu_vtx) {

	TVector3 pt1 = t.Vertex();
	TVector3 pt2 = t.End();

	double nu_vtx[3];
	tpcobj_nu_vtx.XYZ(nu_vtx);
	TVector3 nu(nu_vtx[0], nu_vtx[1], nu_vtx[2]);

	bool reverse = false;

	if ( (pt1-nu).Mag() > (pt2-nu).Mag() )
		reverse = true;

	/*
	   double original_phi = t.Phi();
	   if (reverse)
	   return ::TMath::Pi() - original_phi;
	   else
	   return original_phi;
	 */

	TVector3 dir = pt2 - pt1;
	if (reverse) dir = pt1 - pt2;

	// We are in the plane Z = 0
	dir.SetZ(0);
	TVector3 phi_versor (1, 0, 0);

	double phi = phi_versor.Angle(dir);

	// Just convention
	if (dir.Y() < 0)
		phi = -phi;

	std::cout << "My phi is " << phi << ", track phi is " << t.Phi() << ", reverse is " << (reverse ? "true" : "false") << std::endl;
	return phi;

}

//_________________________________________________________________________________
double PostMCCorrecitons::GetCorrectedCosTheta(recob::Track t, recob::Vertex tpcobj_nu_vtx) {

	TVector3 pt1 = t.Vertex();
	TVector3 pt2 = t.End();

	double nu_vtx[3];
	tpcobj_nu_vtx.XYZ(nu_vtx);
	TVector3 nu(nu_vtx[0], nu_vtx[1], nu_vtx[2]);

	bool reverse = false;

	if ( (pt1-nu).Mag() > (pt2-nu).Mag() )
		reverse = true;

	/*
	   double original_theta = t.Theta();
	   if (reverse)
	   return ::TMath::Pi() - original_theta;
	   else
	   return original_theta;
	 */

	TVector3 dir = pt2 - pt1;
	if (reverse) dir = pt1 - pt2;

	TVector3 theta_versor (0, 0, 1);

	double theta = theta_versor.Angle(dir);

	std::cout << "My theta is " << theta << ", track theta is " << t.Theta() << ", reverse is " << (reverse ? "true" : "false") << std::endl;

	return std::cos(theta);

}


//_________________________________________________________________________________
bool PostMCCorrecitons::PointIsCloseToDeadRegion(double *reco_nu_vtx, int plane_no){

	::art::ServiceHandle<geo::Geometry> geo;

	// Get nearest channel
	raw::ChannelID_t ch = geo->NearestChannel(reco_nu_vtx, plane_no);
	//std::cout << "nearest channel is " << ch << std::endl;

	// Get channel status
	const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
	//std::cout << "ch status " << chanFilt.Status(ch) << std::endl;
	if( chanFilt.Status(ch) < 4) return true;

	// Now check close wires
	for(unsigned int new_ch = ch - 5; new_ch < ch - 5 + 11; new_ch++) {
		//std::cout << "now trying with channel " << new_ch << std::endl;
		if( new_ch < 0 || new_ch >= 8256) continue;
		if( chanFilt.Status(new_ch) < 4 ) return true;
	}

/*
   double new_point_in_tpc[3];
   for (int z_off = -5; z_off < 5; z_off += 2.5) {
    new_point_in_tpc[0] = reco_nu_vtx[0];
    new_point_in_tpc[1] = reco_nu_vtx[1];
    new_point_in_tpc[2] = reco_nu_vtx[2] + z_off;
    //std::cout << "trying with point " << new_point_in_tpc[0] << " " << new_point_in_tpc[1] << " " << new_point_in_tpc[2] << std::endl;
    raw::ChannelID_t ch = geo->NearestChannel(new_point_in_tpc, 2);
    if( chanFilt.Status(ch) < 4) return true;
   }
 */

	return false;

}
