////////////////////////////////////////////////////////////////////////
// Class:       GeometryHelper
// Module Type: filter
// File:        GeometryHelper.h
//
////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

#include "HelperBase.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "TVector3.h"

namespace xsecAna {

class GeometryHelper : public HelperBase {
public:
GeometryHelper() = default;
~GeometryHelper() = default;

int isInside( std::vector<double> P, std::vector< std::vector<double> > V);

/**
 * @brief Determine if the specified point is in the fiducial volume
 *
 * @param x vector of length 3
 * @return True or false
 */
bool isFiducial(const std::vector<double> &x) const;

/**
 * @brief Determine if the specified point is in the fiducial volume
 *
 * @param x TVector3 of 3D location
 * @return True or false
 */
bool isFiducial(const TVector3 &x) const;

/**
 * @brief Determine if the specified point is in the fiducial volume
 *        Not recommended, no array size checking is done.
 *
 * @param x array of 3D location
 * @return True or false
 */
bool isFiducial(const double x[3]) const;

/**
 * @brief Determine if the specified point is in the active volume
 *
 * @param x array of 3D location
 * @return True or false
 */
bool isActive(const std::vector<double> &x) const;

/**
 * @brief Determine if the specified point is in the active volume
 *        Not recommended, no array size checking is done.
 *
 * @param x array of 3D location
 * @return True or false
 */
bool isActive(const double x[3]) const;

/**
 * @brief Compute the 3D distance between two points
 *
 * @param a First Point
 * @param b Second Point
 * @return Returns SQRT( (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2 )
 */
double distance(const std::vector<double> &a,
                const std::vector<double> &b) const;

/**
 * @brief Compute the 3D distance between two points
 *
 * @param a First Point
 * @param b Second Point
 * @return Returns SQRT( (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2 )
 */
double distance(const TVector3 &a, const TVector3 &b) const;

/**
 * @brief      Returns the 3D effective pitch give a direction and a plane
 *
 * @param[in]  direction  The direction
 * @param[in]  pl         Plane of interest (0, 1, 2)
 *
 * @return     The pitch.
 */
double getPitch(const TVector3 &direction, const int &pl) const;

/**
 * @brief      Sets the fiducial volume cuts.
 *
 * @param[in]  m_fidvolXstart  The x start
 * @param[in]  m_fidvolXend    The x end
 * @param[in]  m_fidvolYstart  The y start
 * @param[in]  m_fidvolYend    The y end
 * @param[in]  m_fidvolZstart  The z start
 * @param[in]  m_fidvolZend    The z end
 */
void setFiducialVolumeCuts(float m_fidvolXstart, float m_fidvolXend,
                           float m_fidvolYstart, float m_fidvolYend,
                           float m_fidvolZstart, float m_fidvolZend);

void setRectangleSize(float rectangleLength, float rectangleWidth);

/**
 * @brief      Return the average position of a set of space points
 *
 * @param      spcpnts  The spcpnts
 *
 * @return     3D Vector
 */
TVector3
getAveragePosition(std::vector<art::Ptr<recob::SpacePoint> > &spcpnts);

/**
 * @brief      Builds a rectangle.
 *
 * @param[in]  length  The length
 * @param[in]  width   The width
 * @param      start   The start
 * @param      axis    The axis
 * @param      points  The points
 */
void buildRectangle(double length, double width, std::vector<double> &start,
                    std::vector<double> &axis,
                    std::vector<std::vector<double> > &points);

/**
 * @brief      { function_description }
 *
 * @param      points  The points
 * @param      point   The point
 *
 * @return     { description_of_the_return_value }
 */

int correct_direction(size_t pfp_id, const art::Event &evt, std::string _pfp_producer="pandoraNu");


private:
float m_fidvolXstart;
float m_fidvolXend;
float m_fidvolYstart;
float m_fidvolYend;
float m_fidvolZstart;
float m_fidvolZend;
};
}

#endif
