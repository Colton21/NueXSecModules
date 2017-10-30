#ifndef GEOMETRYHELPER_CXX
#define GEOMETRYHELPER_CXX

#include "GeometryHelper.h"

namespace xsecAna {

bool GeometryHelper::isFiducial(const std::vector<double> &x) const {
	if (x.size() != 3) {
		return false;
	}

	return this->isFiducial(&x[0]);

	// art::ServiceHandle<geo::Geometry> geo;
	// std::vector<double> bnd = {
	//     0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(),
	//     geo->DetHalfHeight(),
	//     0., geo->DetLength()};

	// bool is_x =
	//     x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
	// bool is_y =
	//     x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
	// bool is_z =
	//     x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
	// return is_x && is_y && is_z;
}

bool GeometryHelper::isFiducial(const TVector3 &x) const {
	std::vector<double> _x(3);

	_x[0] = x[0];
	_x[1] = x[1];
	_x[2] = x[2];

	return this->isFiducial(_x);

	// art::ServiceHandle<geo::Geometry> geo;
	// std::vector<double> bnd = {
	//     0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(),
	//     geo->DetHalfHeight(),
	//     0., geo->DetLength()};

	// bool is_x =
	//     x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
	// bool is_y =
	//     x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
	// bool is_z =
	//     x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
	// return is_x && is_y && is_z;
}

bool GeometryHelper::isFiducial(const double x[3]) const {

	art::ServiceHandle<geo::Geometry> geo;
	std::vector<double> bnd = {
		0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
		0., geo->DetLength()
	};


	bool is_x =
	        x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
	bool is_y =
	        x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
	bool is_z =
	        x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
	return is_x && is_y && is_z;
}

bool GeometryHelper::isActive(const std::vector<double> &x) const {
	if (x.size() != 3) {
		return false;
	}

	return this->isActive(&x[0]);
}

bool GeometryHelper::isActive(const double x[3]) const {

	art::ServiceHandle<geo::Geometry> geo;
	std::vector<double> bnd = {
		0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
		0., geo->DetLength()
	};

	bool is_x = x[0] > bnd[0] && x[0] < bnd[1];
	bool is_y = x[1] > bnd[2] && x[1] < bnd[3];
	bool is_z = x[2] > bnd[4] && x[2] < bnd[5];
	return is_x && is_y && is_z;
}

double GeometryHelper::distance(const std::vector<double> &a,
                                const std::vector<double> &b) const {
	if (a.size() != 3 || b.size() != 3) {
		return -1;
	}

	double d = 0;

	for (int i = 0; i < 3; i++) {
		d += pow((a[i] - b[i]), 2);
	}

	return sqrt(d);
}

double GeometryHelper::distance(const TVector3 &a, const TVector3 &b) const {
	return (a - b).Mag();
}

void GeometryHelper::setFiducialVolumeCuts(
        float m_fidvolXstart, float m_fidvolXend, float m_fidvolYstart,
        float m_fidvolYend, float m_fidvolZstart, float m_fidvolZend) {
	this->m_fidvolXstart = m_fidvolXstart;
	this->m_fidvolXend = m_fidvolXend;
	this->m_fidvolYstart = m_fidvolYstart;
	this->m_fidvolYend = m_fidvolYend;
	this->m_fidvolZstart = m_fidvolZstart;
	this->m_fidvolZend = m_fidvolZend;
	std::cout << "Fidvol x start-end " << m_fidvolXstart << " " << m_fidvolXend << std::endl;
	std::cout << "Fidvol y start-end " << m_fidvolYstart << " " << m_fidvolYend << std::endl;
	std::cout << "Fidvol z start-end " << m_fidvolZstart << " " << m_fidvolZend << std::endl;

}

TVector3 GeometryHelper::getAveragePosition(
        std::vector<art::Ptr<recob::SpacePoint> > &spcpnts) {
	double avg_x = 0;
	double avg_y = 0;
	double avg_z = 0;

	for (auto &spcpnt : spcpnts) {
		auto spcpnt_xyz = spcpnt->XYZ();
		avg_x += spcpnt_xyz[0];
		avg_y += spcpnt_xyz[1];
		avg_z += spcpnt_xyz[2];
	}

	return TVector3(avg_x / spcpnts.size(), avg_y / spcpnts.size(),
	                avg_z / spcpnts.size());
}


int GeometryHelper::isInside( std::vector<double> P,
                              std::vector< std::vector<double> > V) {

	int nvert = (int)V.size();

	int i, j, c = 0;
	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((V[i][1]>P[1]) != (V[j][1]>P[1])) &&
		     (P[0] < (V[j][0]-V[i][0]) * (P[1]-V[i][1]) / (V[j][1]-V[i][1]) + V[i][0]) )
			c = !c;
	}
	return c;
}


double GeometryHelper::getPitch(const TVector3 &direction,
                                const int &pl) const {
	// prepare a direction vector for the plane
	TVector3 wireDir = {0., 0., 0.};
	// the direction of the plane is the vector uniting two consecutive wires
	// such that this vector is perpendicular to both wires
	// basically this is the vector perpendicular to the wire length direction,
	// and still in the wire-plane direction
	if (pl == 0)
		wireDir = {0., -sqrt(3) / 2., 1 / 2.};
	else if (pl == 1)
		wireDir = {0., sqrt(3) / 2., 1 / 2.};
	else if (pl == 2)
		wireDir = {0., 0., 1.};

	// cosine between shower direction and plane direction gives the factor
	// by which to divide 0.3, the minimum wire-spacing
	double minWireSpacing = 0.3;
	double cos = wireDir.Dot(direction);
	if (cos < 0)
		cos *= -1;
	cos /= (wireDir.Mag() * direction.Mag());
	// if cosine is 0 the direction is perpendicular and the wire-spacing is
	// infinite
	if (cos == 0)
		return std::numeric_limits<double>::max();
	double pitch = minWireSpacing / cos;
	return pitch;
}


void GeometryHelper::buildRectangle(double length, double width,
                                    std::vector<double> &start,
                                    std::vector<double> &axis,
                                    std::vector<std::vector<double> > &points) {
	double perp_axis[2] = {-axis[1], axis[0]};

	std::vector<double> p1 = {start[0] + perp_axis[0] * width / 2,
		                  start[1] + perp_axis[1] * width / 2};
	std::vector<double> p2 = {p1[0] + axis[0] * length, p1[1] + axis[1] * length};

	std::vector<double> p3 = {start[0] - perp_axis[0] * width / 2,
		                  start[1] - perp_axis[1] * width / 2};
	std::vector<double> p4 = {p3[0] + axis[0] * length, p3[1] + axis[1] * length};




	points.insert(points.end(), {p1, p2, p4, p3});
}

// CORRECT SHOWER DIRECTION THAT SOMETIMES IS INVERTED
int GeometryHelper::correct_direction(size_t pfp_id, const art::Event &evt,
                                      std::string _pfp_producer) {

	auto const &pfparticle_handle =
	        evt.getValidHandle<std::vector<recob::PFParticle> >(_pfp_producer);

	art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt,
	                                                     _pfp_producer);
	std::vector<art::Ptr<recob::SpacePoint> > spcpnts =
	        spcpnts_per_pfpart.at(pfp_id);

	int direction = 1;

	if (spcpnts.size() > 0) {
		art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt,
		                                               _pfp_producer);
		auto const &vertex_obj = vertex_per_pfpart.at(pfp_id);
		double vertex_xyz[3];
		vertex_obj->XYZ(vertex_xyz);
		TVector3 start_vec(vertex_xyz[0], vertex_xyz[1], vertex_xyz[2]);

		art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt,
		                                               _pfp_producer);
		auto const &shower_obj = shower_per_pfpart.at(pfp_id);
		TVector3 shower_vec(shower_obj->Direction().X(),
		                    shower_obj->Direction().Y(),
		                    shower_obj->Direction().Z());

		TVector3 avg_spcpnt = getAveragePosition(spcpnts);

		TVector3 a;
		TVector3 b;
		a = shower_vec;
		b = avg_spcpnt - start_vec;
		double costheta = a.Dot(b);
		direction = costheta >= 0 ? 1 : -1;
	}

	return direction;
}

} // namespace lee

#endif
