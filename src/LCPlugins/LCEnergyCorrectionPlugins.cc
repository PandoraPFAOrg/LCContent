/**
 *  @file   LCContent/src/LCPlugins/LCEnergyCorrectionPlugins.cc
 * 
 *  @brief  Implementation of the lc energy correction plugins class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCHelpers/ReclusterHelper.h"

#include "LCPlugins/LCEnergyCorrectionPlugins.h"

#include <map>

using namespace pandora;

namespace lc_content
{

namespace
{

typedef std::pair<std::string, EnergyCorrectionType> ThetaEnergyCorrectionKey;

class ThetaEnergyCorrectionTable
{
public:
    ThetaEnergyCorrectionTable() = default;

    ThetaEnergyCorrectionTable(const FloatVector &thetaBinEdges, const FloatVector &energyBinEdges, const FloatVector &scaleFactors) :
        m_thetaBinEdges(thetaBinEdges),
        m_energyBinEdges(energyBinEdges),
        m_scaleFactors(scaleFactors)
    {
    }

    FloatVector m_thetaBinEdges;
    FloatVector m_energyBinEdges;
    FloatVector m_scaleFactors;
};

typedef std::map<ThetaEnergyCorrectionKey, ThetaEnergyCorrectionTable> ThetaEnergyCorrectionTableMap;

ThetaEnergyCorrectionTableMap &GetThetaEnergyCorrectionTableMap()
{
    static ThetaEnergyCorrectionTableMap thetaEnergyCorrectionTableMap;
    return thetaEnergyCorrectionTableMap;
}

bool IsStrictlyIncreasing(const FloatVector &values)
{
    if (values.size() < 2)
        return false;

    for (unsigned int i = 1; i < values.size(); ++i)
    {
        if (values.at(i) <= values.at(i - 1))
            return false;
    }

    return true;
}

int FindBin(const FloatVector &edges, const float value)
{
    if (edges.size() < 2)
        return -1;

    if ((value < edges.front()) || (value >= edges.back()))
        return -1;

    for (unsigned int i = 0; i + 1 < edges.size(); ++i)
    {
        if ((edges.at(i) <= value) && (value < edges.at(i + 1)))
            return static_cast<int>(i);
    }

    return -1;
}

float GetCorrection(const FloatVector &thetaBinEdges, const FloatVector &energyBinEdges, const FloatVector &scaleFactors,
    const float theta, const float energy)
{
    const int thetaBin(FindBin(thetaBinEdges, theta));
    const int energyBin(FindBin(energyBinEdges, energy));

    if ((thetaBin < 0) || (energyBin < 0))
        return 1.f;

    const unsigned int nEnergyBins(energyBinEdges.size() - 1);
    return scaleFactors.at(static_cast<unsigned int>(thetaBin) * nEnergyBins + static_cast<unsigned int>(energyBin));
}

} // namespace

void LCEnergyCorrectionPlugins::RegisterThetaEnergyCorrection(const std::string &name, const EnergyCorrectionType energyCorrectionType,
    const FloatVector &thetaBinEdges, const FloatVector &energyBinEdges, const FloatVector &scaleFactors)
{
    if (!IsStrictlyIncreasing(thetaBinEdges) || !IsStrictlyIncreasing(energyBinEdges))
        return;

    const unsigned int nThetaBins(thetaBinEdges.size() - 1);
    const unsigned int nEnergyBins(energyBinEdges.size() - 1);

    if ((0 == nThetaBins) || (0 == nEnergyBins) || (nThetaBins * nEnergyBins != scaleFactors.size()))
        return;

    GetThetaEnergyCorrectionTableMap()[ThetaEnergyCorrectionKey(name, energyCorrectionType)] =
        ThetaEnergyCorrectionTable(thetaBinEdges, energyBinEdges, scaleFactors);
}

//--------------------------------------------------------------------------

float LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(const EnergyCorrectionType energyCorrectionType,
    const CartesianVector &direction, const float energy)
{
    if (direction.GetMagnitude() < std::numeric_limits<float>::epsilon())
        return energy;

    const float cosTheta(std::max(-1.f, std::min(1.f, direction.GetCosOpeningAngle(CartesianVector(0.f, 0.f, 1.f)))));
    const float theta(std::acos(cosTheta));

    for (const ThetaEnergyCorrectionTableMap::value_type &mapEntry : GetThetaEnergyCorrectionTableMap())
    {
        if (mapEntry.first.second != energyCorrectionType)
            continue;

        const ThetaEnergyCorrectionTable &table(mapEntry.second);
        return energy * GetCorrection(table.m_thetaBinEdges, table.m_energyBinEdges, table.m_scaleFactors, theta, energy);
    }

    return energy;
}

LCEnergyCorrectionPlugins::NonLinearityCorrection::NonLinearityCorrection(const FloatVector &inputEnergyCorrectionPoints,
        const FloatVector &outputEnergyCorrectionPoints) :
    m_inputEnergyCorrectionPoints(inputEnergyCorrectionPoints)
{
    const unsigned int nEnergyBins(m_inputEnergyCorrectionPoints.size());

    if (nEnergyBins != outputEnergyCorrectionPoints.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int i = 0; i < nEnergyBins; ++i)
    {
        const float inputEnergy(m_inputEnergyCorrectionPoints.at(i));
        const float outputEnergy(outputEnergyCorrectionPoints.at(i));

        if (std::fabs(inputEnergy) < std::numeric_limits<float>::epsilon())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        m_energyCorrections.push_back(outputEnergy / inputEnergy);
    }

    if (nEnergyBins != m_energyCorrections.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::NonLinearityCorrection::NonLinearityCorrection(const FloatVector &thetaBinEdges, const FloatVector &energyBinEdges,
        const FloatVector &scaleFactors) :
    m_useTwoDimensionalCorrection(true),
    m_thetaBinEdges(thetaBinEdges),
    m_energyBinEdges(energyBinEdges),
    m_energyCorrections(scaleFactors)
{
    if (!IsStrictlyIncreasing(m_thetaBinEdges) || !IsStrictlyIncreasing(m_energyBinEdges))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    const unsigned int nThetaBins(m_thetaBinEdges.size() - 1);
    const unsigned int nEnergyBins(m_energyBinEdges.size() - 1);

    if ((0 == nThetaBins) || (0 == nEnergyBins) || (nThetaBins * nEnergyBins != m_energyCorrections.size()))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCEnergyCorrectionPlugins::NonLinearityCorrection::MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const
{
    if (m_useTwoDimensionalCorrection)
    {
        if (NULL == pCluster)
            return pandora::STATUS_CODE_SUCCESS;

        const CartesianVector &clusterDirection(pCluster->GetFitToAllHitsResult().IsFitSuccessful() ?
            pCluster->GetFitToAllHitsResult().GetDirection() : pCluster->GetInitialDirection());

        if (clusterDirection.GetMagnitude() < std::numeric_limits<float>::epsilon())
            return pandora::STATUS_CODE_SUCCESS;

        const float cosTheta(std::max(-1.f, std::min(1.f, clusterDirection.GetCosOpeningAngle(CartesianVector(0.f, 0.f, 1.f)))));
        correctedEnergy *= this->GetCorrection(std::acos(cosTheta), correctedEnergy);
        return pandora::STATUS_CODE_SUCCESS;
    }

    const unsigned int nEnergyBins(m_energyCorrections.size());

    if (0 == nEnergyBins)
        return pandora::STATUS_CODE_SUCCESS;

    unsigned int index(nEnergyBins);

    for (unsigned int i = 0; i < nEnergyBins; ++i)
    {
        if (correctedEnergy < m_inputEnergyCorrectionPoints.at(i))
        {
            index = i;
            break;
        }
    }

    float correction(1.f);

    if ((0 == index) || (nEnergyBins == index))
    {
        correction = m_energyCorrections.at(std::min(index, nEnergyBins - 1));
    }
    else
    {
        const float lowCorrection(m_energyCorrections.at(index - 1)), highCorrection(m_energyCorrections.at(index));
        const float lowEnergy(m_inputEnergyCorrectionPoints.at(index - 1)), highEnergy(m_inputEnergyCorrectionPoints.at(index));
        correction = lowCorrection + (correctedEnergy - lowEnergy) * (highCorrection - lowCorrection) / (highEnergy - lowEnergy);
    }

    correctedEnergy *= correction;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCEnergyCorrectionPlugins::NonLinearityCorrection::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LCEnergyCorrectionPlugins::NonLinearityCorrection::GetCorrection(const float theta, const float energy) const
{
    const int thetaBin(FindBin(m_thetaBinEdges, theta));
    const int energyBin(FindBin(m_energyBinEdges, energy));

    if ((thetaBin < 0) || (energyBin < 0))
        return 1.f;

    const unsigned int nEnergyBins(m_energyBinEdges.size() - 1);
    return m_energyCorrections.at(static_cast<unsigned int>(thetaBin) * nEnergyBins + static_cast<unsigned int>(energyBin));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::CleanCluster::CleanCluster() :
    m_minCleanHitEnergy(0.5f),
    m_minCleanHitEnergyFraction(0.01f),
    m_minCleanCorrectedHitEnergy(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::CleanCluster::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    const unsigned int firstPseudoLayer(this->GetPandora().GetPlugins()->GetPseudoLayerPlugin()->GetPseudoLayerAtIp());

    const float clusterHadronicEnergy(pCluster->GetHadronicEnergy());

    if (std::fabs(clusterHadronicEnergy) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    bool isFineGranularity(true);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    // Loop over all constituent inner layer fine granularity hits, looking for anomalies
    for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end();
        (layerIter != layerIterEnd) && isFineGranularity; ++layerIter)
    {
        const unsigned int pseudoLayer(layerIter->first);

        for (CaloHitList::const_iterator hitIter = layerIter->second->begin(), hitIterEnd = layerIter->second->end();
            hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            if (this->GetPandora().GetGeometry()->GetHitTypeGranularity((*hitIter)->GetHitType()) > FINE)
            {
                isFineGranularity = false;
                break;
            }

            const float hitHadronicEnergy(pCaloHit->GetHadronicEnergy());

            if ((hitHadronicEnergy > m_minCleanHitEnergy) && (hitHadronicEnergy / clusterHadronicEnergy > m_minCleanHitEnergyFraction))
            {
                // Calculate new energy from surrounding layers
                float energyInPreviousLayer(0.);

                if (pseudoLayer > firstPseudoLayer)
                    energyInPreviousLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer - 1);

                float energyInNextLayer(0.);

                if (pseudoLayer < std::numeric_limits<unsigned int>::max())
                    energyInNextLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer + 1);

                const float energyInCurrentLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer);

                // Calculate new energy estimate for hit and update cluster best energy estimate
                float energyInAdjacentLayers(energyInPreviousLayer + energyInNextLayer);

                if (pseudoLayer > firstPseudoLayer)
                    energyInAdjacentLayers /= 2.f;

                float newHitHadronicEnergy(energyInAdjacentLayers - energyInCurrentLayer + hitHadronicEnergy);
                newHitHadronicEnergy = std::max(newHitHadronicEnergy, m_minCleanCorrectedHitEnergy);

                if (newHitHadronicEnergy < hitHadronicEnergy)
                    correctedHadronicEnergy += newHitHadronicEnergy - hitHadronicEnergy;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LCEnergyCorrectionPlugins::CleanCluster::GetHadronicEnergyInLayer(const OrderedCaloHitList &orderedCaloHitList, const unsigned int pseudoLayer) const
{
    OrderedCaloHitList::const_iterator iter = orderedCaloHitList.find(pseudoLayer);

    float hadronicEnergy(0.f);

    if (iter != orderedCaloHitList.end())
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            hadronicEnergy += (*hitIter)->GetHadronicEnergy();
        }
    }

    return hadronicEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::CleanCluster::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCleanHitEnergy", m_minCleanHitEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCleanHitEnergyFraction", m_minCleanHitEnergyFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCleanCorrectedHitEnergy", m_minCleanCorrectedHitEnergy));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::ScaleHotHadrons::ScaleHotHadrons() :
    m_minHitsForHotHadron(5),
    m_maxHitsForHotHadron(100),
    m_hotHadronInnerLayerCut(10),
    m_hotHadronMipFractionCut(0.4),
    m_hotHadronNHitsCut(50),
    m_hotHadronMipsPerHit(15.f),
    m_scaledHotHadronMipsPerHit(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::ScaleHotHadrons::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    const unsigned int nHitsInCluster(pCluster->GetNCaloHits());

    // Initial hot hadron cuts
    if ((nHitsInCluster < m_minHitsForHotHadron) || (nHitsInCluster > m_maxHitsForHotHadron))
        return STATUS_CODE_SUCCESS;

    if ((pCluster->GetInnerPseudoLayer() < m_hotHadronInnerLayerCut) && (pCluster->GetMipFraction() < m_hotHadronMipFractionCut) &&
        (nHitsInCluster > m_hotHadronNHitsCut))
    {
        return STATUS_CODE_SUCCESS;
    }

    // Finally, check the number of mips per hit
    float clusterMipEnergy(0.);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end();
        layerIter != layerIterEnd; ++layerIter)
    {
        for (CaloHitList::const_iterator hitIter = layerIter->second->begin(), hitIterEnd = layerIter->second->end();
            hitIter != hitIterEnd; ++hitIter)
        {
            clusterMipEnergy += (*hitIter)->GetMipEquivalentEnergy();
        }
    }

    const float meanMipsPerHit(clusterMipEnergy / static_cast<float>(nHitsInCluster));

    if ((meanMipsPerHit > 0.f) && (meanMipsPerHit > m_hotHadronMipsPerHit))
        correctedHadronicEnergy *= m_scaledHotHadronMipsPerHit / meanMipsPerHit;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::ScaleHotHadrons::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForHotHadron", m_minHitsForHotHadron));

    if (0 == m_minHitsForHotHadron)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitsForHotHadron", m_maxHitsForHotHadron));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronInnerLayerCut", m_hotHadronInnerLayerCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronMipFractionCut", m_hotHadronMipFractionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronNHitsCut", m_hotHadronNHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HotHadronMipsPerHit", m_hotHadronMipsPerHit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ScaledHotHadronMipsPerHit", m_scaledHotHadronMipsPerHit));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCEnergyCorrectionPlugins::MuonCoilCorrection::MuonCoilCorrection() :
    m_muonHitEnergy(0.5f),
    m_coilEnergyLossCorrection(10.f),
    m_minMuonHitsInInnerLayer(3),
    m_coilEnergyCorrectionChi(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::MuonCoilCorrection::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    bool containsMuonHit(false);
    unsigned int nMuonHitsInInnerLayer(0);
    unsigned int muonInnerLayer(std::numeric_limits<unsigned int>::max());

    // Extract muon-based properties from the cluster
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            if ((*hitIter)->GetHitType() == MUON)
            {
                containsMuonHit = true;
                ++nMuonHitsInInnerLayer;
            }
        }

        if (containsMuonHit)
        {
            muonInnerLayer = iter->first;
            break;
        }
    }

    if (!containsMuonHit)
        return STATUS_CODE_SUCCESS;;

    // Check whether energy deposits are likely to have been lost in coil region
    const CartesianVector muonInnerLayerCentroid(pCluster->GetCentroid(muonInnerLayer));
    const float centroidX(muonInnerLayerCentroid.GetX()), centroidY(muonInnerLayerCentroid.GetY());

    const float muonInnerLayerRadius(std::sqrt(centroidX * centroidX + centroidY * centroidY));
    const float coilInnerRadius(this->GetPandora().GetGeometry()->GetSubDetector(COIL).GetInnerRCoordinate());

    if (muonInnerLayerRadius < coilInnerRadius)
        return STATUS_CODE_SUCCESS;;

    const TrackList &trackList(pCluster->GetAssociatedTrackList());

    if (pCluster->GetInnerPseudoLayer() == muonInnerLayer)
    {
        // Energy correction for standalone muon cluster
        correctedHadronicEnergy += m_muonHitEnergy * static_cast<float>(nMuonHitsInInnerLayer);
    }
    else if (trackList.empty())
    {
        // Energy correction for neutral hadron cluster spilling into coil and muon detectors
        correctedHadronicEnergy += m_coilEnergyLossCorrection;
    }
    else
    {
        // Energy correction for charged hadron cluster spilling into coil and muon detectors
        if (nMuonHitsInInnerLayer < m_minMuonHitsInInnerLayer)
            return STATUS_CODE_SUCCESS;;

        float trackEnergySum(0.);

        for (TrackList::const_iterator iter = trackList.begin(), iterEnd = trackList.end(); iter != iterEnd; ++iter)
        {
            trackEnergySum += (*iter)->GetEnergyAtDca();
        }

        const float oldChi(ReclusterHelper::GetTrackClusterCompatibility(this->GetPandora(), correctedHadronicEnergy, trackEnergySum));
        const float newChi(ReclusterHelper::GetTrackClusterCompatibility(this->GetPandora(), correctedHadronicEnergy + m_coilEnergyLossCorrection, trackEnergySum));

        if ((oldChi < m_coilEnergyCorrectionChi) && (std::fabs(newChi) < std::fabs(oldChi)))
        {
            correctedHadronicEnergy += m_coilEnergyLossCorrection;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEnergyCorrectionPlugins::MuonCoilCorrection::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonHitEnergy", m_muonHitEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CoilEnergyLossCorrection", m_coilEnergyLossCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMuonHitsInInnerLayer", m_minMuonHitsInInnerLayer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CoilEnergyCorrectionChi", m_coilEnergyCorrectionChi));

    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
