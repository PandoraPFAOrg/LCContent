/**
 *  @file   LCContent/src/LCUtility/EventPreparationAlgorithm.cc
 * 
 *  @brief  Implementation of the event preparation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCUtility/EventPreparationAlgorithm.h"

using namespace pandora;

namespace lc_content
{

StatusCode EventPreparationAlgorithm::Run()
{
    // Filter current track list to select tracks to be used during clustering
    const TrackList *pCurrentTrackList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCurrentTrackList));

    TrackList clusteringTrackList;

    for (TrackList::const_iterator iter = pCurrentTrackList->begin(), iterEnd = pCurrentTrackList->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->GetDaughterList().empty())
            clusteringTrackList.push_back(*iter);
    }

    // Save the filtered list and set it to be the current list for subsequent algorithms
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, clusteringTrackList, m_outputTrackListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Track>(*this, m_replacementTrackListName));

    // Split input calo hit list into ecal/hcal and muon calo hits
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    CaloHitList caloHitList, muonCaloHitList;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        if (MUON == (*hitIter)->GetHitType())
        {
            muonCaloHitList.push_back(*hitIter);
        }
        else
        {
            caloHitList.push_back(*hitIter);
        }
    }

    // Save the lists, setting the ecal/hcal list to be the current list for subsequent algorithms
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, muonCaloHitList, m_outputMuonCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList, m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_replacementCaloHitListName));

    pdebug() << "Initial number of tracks: " << pCurrentTrackList->size() << std::endl;
    pdebug() << "Tracks for clustering: " << clusteringTrackList.size() << std::endl;
    pdebug() << "Initial number of calo hits: " << pCaloHitList->size() << std::endl;
    pdebug() << "- ecal/hcal hits: " << caloHitList.size() << std::endl;
    pdebug() << "- muon hits: " << muonCaloHitList.size() << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputTrackListName", m_outputTrackListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListName", m_outputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMuonCaloHitListName", m_outputMuonCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplacementTrackListName", m_replacementTrackListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplacementCaloHitListName", m_replacementCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
