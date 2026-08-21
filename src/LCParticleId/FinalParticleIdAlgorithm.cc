/**
 *  @file   LCContent/src/LCParticleId/FinalParticleIdAlgorithm.cc
 * 
 *  @brief  Implementation of the final particle id algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCParticleId/FinalParticleIdAlgorithm.h"

using namespace pandora;

namespace lc_content
{

StatusCode FinalParticleIdAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    pdebug() << "Number of PFOs: " << pPfoList->size() << std::endl;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end();
        iter != iterEnd; ++iter)
    {
        pdebug() << "Processing PFO" << std::endl;

        const ParticleFlowObject *const pPfo = *iter;

        const TrackList &trackList(pPfo->GetTrackList());
        const ClusterList &clusterList(pPfo->GetClusterList());

        // Consider only pfos with a single cluster and no track sibling relationships
        if ((clusterList.size() != 1) || (trackList.empty()) || this->ContainsSiblingTrack(trackList)) {
            pdebug() << "n(cluster): " << clusterList.size() << std::endl;
            pdebug() << "track list is empty? " << trackList.empty() << std::endl;
            pdebug() << "PFO contains sibling track? " << this->ContainsSiblingTrack(trackList) << std::endl;
            pdebug() << "--> Skipping PFO object" << std::endl;
            continue;
        }
        const int charge(pPfo->GetCharge());

        if (0 == charge) {
            perror() << "Charge is zero" << std::endl;
            return STATUS_CODE_FAILURE;
        }

        // Ignore particle flow objects already tagged as electrons or muons
        if ((std::abs(pPfo->GetParticleId()) == E_MINUS) || (std::abs(pPfo->GetParticleId()) == MU_MINUS)) {
            pdebug() << "Charged PFO already tagged as electron or muon, skipping" << std::endl;
            continue;
        }

        const Cluster *const pCluster = *(clusterList.begin());
        const ParticleId *const pParticleId(PandoraContentApi::GetPlugins(*this)->GetParticleId());

        // Run electron id, followed by muon id
        PandoraContentApi::ParticleFlowObject::Metadata metadata;

        if (pParticleId->IsElectron(pCluster))
        {
            pdebug() << "Charged PFO tagged as electron" << std::endl;
            metadata.m_particleId = (charge < 0) ? E_MINUS : E_PLUS;
        }
        else if (pParticleId->IsMuon(pCluster))
        {
            pdebug() << "Charged PFO tagged as muon" << std::endl;
            metadata.m_particleId = (charge < 0) ? MU_MINUS : MU_PLUS;
        }
        else {
            pdebug() << "Charged PFO neither tagged as electron nor muon" << std::endl;
        }

        if (metadata.m_particleId.IsInitialized())
        {
            metadata.m_mass = PdgTable::GetParticleMass(metadata.m_particleId.Get());
            metadata.m_energy = std::sqrt(metadata.m_mass.Get() * metadata.m_mass.Get() + pPfo->GetMomentum().GetMagnitudeSquared());
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FinalParticleIdAlgorithm::ContainsSiblingTrack(const TrackList &trackList) const
{
    for (TrackList::const_iterator iter = trackList.begin(), iterEnd = trackList.end(); iter != iterEnd; ++iter)
    {
        if (!(*iter)->GetSiblingList().empty())
        {
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode FinalParticleIdAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
