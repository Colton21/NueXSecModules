#include "AnaHelper.h"

#ifndef RECOTRUEMATCHER_H
#define RECOTRUEMATCHER_H

class PostMCCorrecitons {

public:

static double GetCorrectedPhi(recob::Track t, recob::Vertex tpcobj_nu_vtx);

/**
 *  @brief Returns the value of the cos(theta) angle for the track considering the nu vertex
 *
 *  @param track the recob::Track
 *  @param tpcobj_nu_vtx the recob::Vertex neutrino vertex */
static double GetCorrectedCosTheta(recob::Track t, recob::Vertex tpcobj_nu_vtx);

/**
 *  @brief Returns true if the point passed is close to a dead region
 *
 *  @param reco_nu_vtx the point to check
 *  @param plane_no the plane we want to check with  */
static bool PointIsCloseToDeadRegion(double *reco_nu_vtx, int plane_no);

};

#endif
