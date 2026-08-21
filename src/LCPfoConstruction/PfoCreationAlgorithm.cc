/**
 *  @file   LCContent/src/LCPfoConstruction/PfoCreationAlgorithm.cc
 *
 *  @brief  Implementation of the pfo creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCPfoConstruction/PfoCreationAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lc_content
{

PfoCreationAlgorithm::PfoCreationAlgorithm() :
    m_shouldCreateTrackBasedPfos(true),
    m_shouldCreateNeutralPfos(true),
    m_minClusterHadronicEnergy(0.25f),
    m_minClusterElectromagneticEnergy(0.f),
    m_minHitsInCluster(5),
    m_allowSingleLayerClusters(false),
    m_photonPositionAlgorithm(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    if (m_shouldCreateTrackBasedPfos)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTrackBasedPfos());

    if (m_shouldCreateNeutralPfos)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateNeutralPfos());

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::CreateTrackBasedPfos() const
{
    pdebug() << "Creating track based PFOs" << std::endl;

    // Current track list should contain those tracks selected as "good" by the track preparation algorithm
    const TrackList *pTrackList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pTrackList));

    pdebug() << "Looping over track list" << std::endl;
    for (TrackList::const_iterator iter = pTrackList->begin(), iterEnd = pTrackList->end(); iter != iterEnd; ++iter)
    {
        const Track *const pTrack = *iter;
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;

        // Walk along list of associated daughter/sibling tracks and their cluster associations
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PopulateTrackBasedPfo(pTrack, pfoParameters));

        // Specify the pfo parameters
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SetTrackBasedPfoParameters(pTrack, pfoParameters));

        // Create the pfo
        const ParticleFlowObject *pPfo(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::PopulateTrackBasedPfo(const Track *const pTrack, PfoParameters &pfoParameters, const bool readSiblingInfo) const
{
    // Add track to the pfo
    if (pfoParameters.m_trackList.end() != std::find(pfoParameters.m_trackList.begin(), pfoParameters.m_trackList.end(), pTrack))
        return STATUS_CODE_SUCCESS;

    pfoParameters.m_trackList.push_back(pTrack);
    pdebug() << "Added track with momentum at DCA: " << pTrack->GetMomentumAtDca() << std::endl;

    // Add any cluster associated with this track to the pfo
    try
    {
        const Cluster *const pAssociatedCluster(pTrack->GetAssociatedCluster());
        pfoParameters.m_clusterList.push_back(pAssociatedCluster);
        pdebug() <<  "Added associated cluster " << pTrack->GetAssociatedCluster() << std::endl;
    }
    catch (StatusCodeException &)
    {
    }

    // Consider any sibling tracks
    if (readSiblingInfo)
    {
        pdebug() << "Adding sibling tracks" << std::endl;
        const TrackList &siblingTrackList(pTrack->GetSiblingList());

        for (TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PopulateTrackBasedPfo(*iter, pfoParameters, false));
        }
    }

    // Consider any daughter tracks
    pdebug() << "Adding daughter tracks" << std::endl;
    const TrackList &daughterTrackList(pTrack->GetDaughterList());
    for (TrackList::const_iterator iter = daughterTrackList.begin(), iterEnd = daughterTrackList.end(); iter != iterEnd; ++iter)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PopulateTrackBasedPfo(*iter, pfoParameters));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::SetTrackBasedPfoParameters(const Track *const pTrack, PfoParameters &pfoParameters) const
{
    const bool hasParent(!pTrack->GetParentList().empty());

    if (hasParent)
        return STATUS_CODE_NOT_ALLOWED;

    const bool hasSibling(!pTrack->GetSiblingList().empty());
    const bool hasDaughter(!pTrack->GetDaughterList().empty());

    if (hasSibling && hasDaughter)
        return STATUS_CODE_NOT_ALLOWED;

    if (hasSibling)
        return this->SetSiblingTrackBasedPfoParameters(pTrack, pfoParameters);

    if (hasDaughter)
        return this->SetDaughterTrackBasedPfoParameters(pTrack, pfoParameters);

    return this->SetSimpleTrackBasedPfoParameters(pTrack, pfoParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::SetSiblingTrackBasedPfoParameters(const Track *const pTrack, PfoParameters &pfoParameters) const
{
    pdebug() << "Setting PFO parameters based on sibling tracks" << std::endl;

    int charge(0);
    float energy(0.f);
    CartesianVector momentum(0.f, 0.f, 0.f);

    TrackList fullSiblingTrackList(pTrack->GetSiblingList());
    fullSiblingTrackList.push_back(pTrack);

    for (TrackList::const_iterator iter = fullSiblingTrackList.begin(), iterEnd = fullSiblingTrackList.end(); iter != iterEnd; ++iter)
    {
        const Track *const pSiblingTrack = *iter;
        charge += pSiblingTrack->GetCharge();

        if (!pSiblingTrack->CanFormPfo() && !pSiblingTrack->CanFormClusterlessPfo())
            continue;

        // ATTN Assume sibling-track-based pfos represent pair-production
        const float electronMass(PdgTable::GetParticleMass(E_MINUS));
        energy += std::sqrt(electronMass * electronMass + pSiblingTrack->GetMomentumAtDca().GetMagnitudeSquared());
        momentum += pSiblingTrack->GetMomentumAtDca();
    }

    if (energy < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_INITIALIZED;

    pfoParameters.m_energy = energy;
    pfoParameters.m_momentum = momentum;
    pfoParameters.m_mass = std::sqrt(std::max(energy * energy - momentum.GetDotProduct(momentum), 0.f));
    pfoParameters.m_charge = charge;
    pfoParameters.m_particleId = PHOTON;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::SetDaughterTrackBasedPfoParameters(const Track *const pTrack, PfoParameters &pfoParameters) const
{
    pdebug() << "Setting PFO parameters based on daughter tracks" << std::endl;

    int daughterCharge(0);
    float energy(0.f);
    CartesianVector momentum(0.f, 0.f, 0.f);

    const TrackList &daughterTrackList(pTrack->GetDaughterList());
    const unsigned int nDaughters(daughterTrackList.size());

    for (TrackList::const_iterator iter = daughterTrackList.begin(), iterEnd = daughterTrackList.end(); iter != iterEnd; ++iter)
    {
        const Track *const pDaughterTrack = *iter;

        if (!pDaughterTrack->CanFormPfo() && !pDaughterTrack->CanFormClusterlessPfo())
            continue;

        daughterCharge += pDaughterTrack->GetCharge();
        energy += pDaughterTrack->GetEnergyAtDca();
        momentum += pDaughterTrack->GetMomentumAtDca();
    }

    pfoParameters.m_energy = energy;
    pfoParameters.m_momentum = momentum;
    pfoParameters.m_mass = std::sqrt(std::max(energy * energy - momentum.GetDotProduct(momentum), 0.f));
    pfoParameters.m_charge = (nDaughters > 1) ? pTrack->GetCharge() : daughterCharge;
    pfoParameters.m_particleId = (pfoParameters.m_charge.Get() > 0) ? PI_PLUS : PI_MINUS;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::SetSimpleTrackBasedPfoParameters(const Track *const pTrack, PfoParameters &pfoParameters) const
{
    pdebug() << "Setting PFO parameters based on single track" << std::endl;
    pfoParameters.m_energy = pTrack->GetEnergyAtDca();
    pfoParameters.m_momentum = pTrack->GetMomentumAtDca();
    pfoParameters.m_mass = pTrack->GetMass();
    pfoParameters.m_charge = pTrack->GetCharge();
    pfoParameters.m_particleId = (pTrack->GetCharge() > 0) ? PI_PLUS : PI_MINUS;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::CreateNeutralPfos() const
{
    pdebug() << "Creating neutral PFOs from full list of clusters" << std::endl;

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Examine clusters with no associated tracks to form neutral pfos
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!pCluster->GetAssociatedTrackList().empty()) {
            pdebug() << "Skipping cluster with associated tracks" << std::endl;
            continue;
        }

        if (pCluster->GetNCaloHits() < m_minHitsInCluster) {
            pdebug() << "Skipping cluster with not enough hits: " << pCluster->GetNCaloHits() << std::endl;
            continue;
        }

        const bool isPhoton(pCluster->PassPhotonId(this->GetPandora()));
        float clusterEnergy(isPhoton ? pCluster->GetCorrectedElectromagneticEnergy(this->GetPandora()) : pCluster->GetCorrectedHadronicEnergy(this->GetPandora()));
        pdebug() << "Cluster is " << (isPhoton ? "photon" : "neutral hadron") << std::endl;
        pdebug() << "Corrected cluster energy: " << clusterEnergy << std::endl;

        // Veto non-photon clusters below hadronic energy threshold and those occupying a single layer
        if (!isPhoton)
        {
            if (clusterEnergy < m_minClusterHadronicEnergy) {
                pdebug() << "Skipping cluster flagged as hadron due to energy below threshold: " << m_minClusterHadronicEnergy << std::endl;
                continue;
            }

            if (!m_allowSingleLayerClusters && (pCluster->GetInnerPseudoLayer() == pCluster->GetOuterPseudoLayer())) {
                pdebug() << "Skipping cluster flagged as hadron due since inner and outer layer are the same" << std::endl;
                continue;
            }
        }
        else
        {
            if (clusterEnergy < m_minClusterElectromagneticEnergy) {
                pdebug() << "Skipping cluster flagged as photon due to energy below threshold: " << m_minClusterElectromagneticEnergy << std::endl;
                continue;
            }
        }

        pdebug() << "Good cluster found, setting parameters" << std::endl;

        // Specify the pfo parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = (isPhoton ? PHOTON : NEUTRON);
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = (isPhoton ? PdgTable::GetParticleMass(PHOTON) : PdgTable::GetParticleMass(NEUTRON));
        pfoParameters.m_energy = clusterEnergy;
        pfoParameters.m_clusterList.push_back(pCluster);

        // Photon position: 0) unweighted inner centroid, 1) energy-weighted inner centroid, 2+) energy-weighted centroid for all layers
        CartesianVector positionVector(0.f, 0.f, 0.f);
        const unsigned int clusterInnerLayer(pCluster->GetInnerPseudoLayer());

        if (!isPhoton || (0 == m_photonPositionAlgorithm))
        {
            positionVector = pCluster->GetCentroid(clusterInnerLayer);
        }
        else if (1 == m_photonPositionAlgorithm)
        {
            positionVector = this->GetEnergyWeightedCentroid(pCluster, clusterInnerLayer, clusterInnerLayer);
        }
        else
        {
            positionVector = this->GetEnergyWeightedCentroid(pCluster, clusterInnerLayer, pCluster->GetOuterPseudoLayer());
        }

        const CartesianVector momentum(positionVector.GetUnitVector() * clusterEnergy);
        pfoParameters.m_momentum = momentum;

        const ParticleFlowObject *pPfo(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector PfoCreationAlgorithm::GetEnergyWeightedCentroid(const Cluster *const pCluster, const unsigned int innerPseudoLayer,
    const unsigned int outerPseudoLayer) const
{
    float energySum(0.f);
    CartesianVector energyPositionSum(0.f, 0.f, 0.f);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        if (iter->first > outerPseudoLayer)
            break;

        if (iter->first < innerPseudoLayer)
            continue;

        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const float electromagneticEnergy((*hitIter)->GetElectromagneticEnergy());
            energySum += electromagneticEnergy;
            energyPositionSum += ((*hitIter)->GetPositionVector() * electromagneticEnergy);
        }
    }

    if (energySum < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return (energyPositionSum * (1.f / energySum));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldCreateTrackBasedPfos", m_shouldCreateTrackBasedPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldCreateNeutralPfos", m_shouldCreateNeutralPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterHadronicEnergy", m_minClusterHadronicEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterElectromagneticEnergy", m_minClusterElectromagneticEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AllowSingleLayerClusters", m_allowSingleLayerClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PhotonPositionAlgorithm", m_photonPositionAlgorithm));

    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
