#include "RecoTrueHelper.h"

namespace nue_xsec
{
//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoNeutrinosToHits, const lar_pandora::HitsToMCTruth &trueHitsToNeutrinos,
                                          lar_pandora::MCTruthToPFParticles &matchedNeutrinos, lar_pandora::MCTruthToHits &matchedNeutrinoHits) const
{
	PFParticleSet recoVeto; MCTruthSet trueVeto;

	this->GetRecoToTrueMatches(recoNeutrinosToHits, trueHitsToNeutrinos, matchedNeutrinos, matchedNeutrinoHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoNeutrinosToHits, const lar_pandora::HitsToMCTruth &trueHitsToNeutrinos,
                                          lar_pandora::MCTruthToPFParticles &matchedNeutrinos, lar_pandora::MCTruthToHits &matchedNeutrinoHits,
                                          lar_pandora::PFParticleSet &vetoReco, lar_pandora::MCTruthSet &vetoTrue) const
{
	bool foundMatches(false);

	for (PFParticlesToHits::const_iterator iter1 = recoNeutrinosToHits.begin(), iterEnd1 = recoNeutrinosToHits.end();
	     iter1 != iterEnd1; ++iter1)
	{
		const art::Ptr<recob::PFParticle> recoNeutrino = iter1->first;
		if (vetoReco.count(recoNeutrino) > 0)
			continue;

		const HitVector &hitVector = iter1->second;

		MCTruthToHits truthContributionMap;

		for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
		{
			const art::Ptr<recob::Hit> hit = *iter2;

			HitsToMCTruth::const_iterator iter3 = trueHitsToNeutrinos.find(hit);
			if (trueHitsToNeutrinos.end() == iter3)
				continue;

			const art::Ptr<simb::MCTruth> trueNeutrino = iter3->second;
			if (vetoTrue.count(trueNeutrino) > 0)
				continue;

			truthContributionMap[trueNeutrino].push_back(hit);
		}

		MCTruthToHits::const_iterator mIter = truthContributionMap.end();

		for (MCTruthToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
		     iter4 != iterEnd4; ++iter4)
		{
			if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
			{
				mIter = iter4;
			}
		}

		if (truthContributionMap.end() != mIter)
		{
			const art::Ptr<simb::MCTruth> trueNeutrino = mIter->first;

			MCTruthToHits::const_iterator iter5 = matchedNeutrinoHits.find(trueNeutrino);

			if ((matchedNeutrinoHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
			{
				matchedNeutrinos[trueNeutrino] = recoNeutrino;
				matchedNeutrinoHits[trueNeutrino] = mIter->second;
				foundMatches = true;
			}
		}
	}

	if (!foundMatches)
		return;

	for (MCTruthToPFParticles::const_iterator pIter = matchedNeutrinos.begin(), pIterEnd = matchedNeutrinos.end();
	     pIter != pIterEnd; ++pIter)
	{
		vetoTrue.insert(pIter->first);
		vetoReco.insert(pIter->second);
	}

	if (_recursiveMatching)
		this->GetRecoToTrueMatches(recoNeutrinosToHits, trueHitsToNeutrinos, matchedNeutrinos, matchedNeutrinoHits, vetoReco, vetoTrue);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                          lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits) const
{
	PFParticleSet recoVeto; MCParticleSet trueVeto;

	this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                          lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits,
                                          lar_pandora::PFParticleSet &vetoReco, lar_pandora::MCParticleSet &vetoTrue) const
{
	bool foundMatches(false);

	for (PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
	     iter1 != iterEnd1; ++iter1)
	{
		const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
		if (vetoReco.count(recoParticle) > 0)
			continue;

		const HitVector &hitVector = iter1->second;

		MCParticlesToHits truthContributionMap;

		for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
		{
			const art::Ptr<recob::Hit> hit = *iter2;

			HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
			if (trueHitsToParticles.end() == iter3)
				continue;

			const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
			if (vetoTrue.count(trueParticle) > 0)
				continue;

			truthContributionMap[trueParticle].push_back(hit);
		}

		MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

		for (MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
		     iter4 != iterEnd4; ++iter4)
		{
			if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
			{
				mIter = iter4;
			}
		}

		if (truthContributionMap.end() != mIter)
		{
			const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

			MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);

			if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
			{
				matchedParticles[trueParticle] = recoParticle;
				matchedHits[trueParticle] = mIter->second;
				foundMatches = true;
			}
		}
	}

	if (!foundMatches)
		return;

	for (MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
	     pIter != pIterEnd; ++pIter)
	{
		vetoTrue.insert(pIter->first);
		vetoReco.insert(pIter->second);
	}

	if (_recursiveMatching)
		this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue);
}



//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::BuildRecoParticleMap(const lar_pandora::PFParticleVector &particleVector, lar_pandora::PFParticleMap &particleMap) const
{
	for (PFParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
	{
		const art::Ptr<recob::PFParticle> particle = *iter;
		particleMap[particle->Self()] = particle;
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::BuildTrueParticleMap(const lar_pandora::MCParticleVector &particleVector, lar_pandora::MCParticleMap &particleMap) const
{
	for (MCParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
	{
		const art::Ptr<simb::MCParticle> particle = *iter;
		particleMap[particle->TrackId()] = particle;
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------

int recotruehelper::CountHitsByType(const int view, const lar_pandora::HitVector &hitVector) const
{
	int nHits(0);

	for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
	{
		const art::Ptr<recob::Hit> hit = *iter;
		if (hit->View() == view)
			++nHits;
	}

	return nHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetStartAndEndPoints(const art::Ptr<simb::MCParticle> particle, int &startT, int &endT) const
{
	art::ServiceHandle<geo::Geometry> theGeometry;

	bool foundStartPosition(false);

	const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

	for (int nt = 0; nt < numTrajectoryPoints; ++nt)
	{
		try
		{
			double pos[3] = {particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
			unsigned int which_tpc(std::numeric_limits<unsigned int>::max());
			unsigned int which_cstat(std::numeric_limits<unsigned int>::max());
			theGeometry->PositionToTPC(pos, which_tpc, which_cstat);

			// TODO: Apply fiducial cut due to readout window

			endT = nt;
			if (!foundStartPosition)
			{
				startT = endT;
				foundStartPosition = true;
			}
		}
		catch (cet::exception &e) {
			continue;
		}
	}

	if (!foundStartPosition)
		throw cet::exception("LArPandora");
}

//------------------------------------------------------------------------------------------------------------------------------------------

double recotruehelper::GetLength(const art::Ptr<simb::MCParticle> particle, const int startT, const int endT) const
{
	if (endT <= startT)
		return 0.0;

	double length(0.0);

	for (int nt = startT; nt < endT; ++nt)
	{
		const double dx(particle->Vx(nt+1) - particle->Vx(nt));
		const double dy(particle->Vy(nt+1) - particle->Vy(nt));
		const double dz(particle->Vz(nt+1) - particle->Vz(nt));
		length += sqrt(dx * dx + dy * dy + dz * dz);
	}

	return length;
}

}//end namespace lar_pandora
