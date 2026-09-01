/**
 *  @file   LCContent/src/LCUtility/TrainingSoftwareCompensation.cc
 *
 *  @brief  Implementation of the training software compensation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCUtility/TrainingSoftwareCompensation.h"

using namespace pandora;

namespace lc_content {

TrainingSoftwareCompensation::TrainingSoftwareCompensation()
    : m_myRootFileName(""), m_trainingTreeName("SoftwareCompensationTrainingTree") {}

//------------------------------------------------------------------------------------------------------------------------------------------

TrainingSoftwareCompensation::~TrainingSoftwareCompensation() {
  try {
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_myRootFileName, "UPDATE"));
  } catch (const StatusCodeException&) {
    std::cout << "TrainingSoftwareCompensation: Unable to write tree " << m_trainingTreeName << " to file "
              << m_myRootFileName << std::endl;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrainingSoftwareCompensation::Run() {
#ifdef MONITORING
  const pandora::PfoList* pPfoList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

  if (pPfoList->size() != 1)
    return STATUS_CODE_SUCCESS;

  const pandora::Pfo* const pPfo = *pPfoList->begin();
  const pandora::ClusterList* pClusterList = &pPfo->GetClusterList();

  if (pClusterList->size() != 1)
    return STATUS_CODE_SUCCESS;

  const Cluster* const pCluster = *pClusterList->begin();
  pandora::CaloHitList clusterCaloHitList;
  pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);
  clusterCaloHitList.insert(clusterCaloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(),
                            pCluster->GetIsolatedCaloHitList().end());

  const float rawEnergyOfCluster(pCluster->GetHadronicEnergy());
  const float pfoEnergy(pPfo->GetEnergy());

  FloatVector cellSize0, cellSize1, cellThickness, hitEnergies;
  IntVector hitType, pseudoLayers, isolatedHit;

  for (pandora::CaloHitList::const_iterator hitIter = clusterCaloHitList.begin(), endhitIter = clusterCaloHitList.end();
       endhitIter != hitIter; ++hitIter) {
    const pandora::CaloHit* pCaloHit = *hitIter;
    const float cellSize0ToAdd(pCaloHit->GetCellSize0());
    const float cellSize1ToAdd(pCaloHit->GetCellSize1());
    const float cellThicknessToAdd(pCaloHit->GetCellThickness());
    const float cellHadronicEnergy(pCaloHit->GetHadronicEnergy());
    const int pseudoLayer(pCaloHit->GetPseudoLayer());
    const int isIsolated(pCaloHit->IsIsolated());

    cellSize0.push_back(cellSize0ToAdd);
    cellSize1.push_back(cellSize1ToAdd);
    cellThickness.push_back(cellThicknessToAdd);
    hitEnergies.push_back(cellHadronicEnergy);
    pseudoLayers.push_back(pseudoLayer);
    isolatedHit.push_back(isIsolated);

    if (HCAL == pCaloHit->GetHitType()) {
      hitType.push_back(2);
    } else if (ECAL == pCaloHit->GetHitType()) {
      hitType.push_back(1);
    } else {
      hitType.push_back(3);
    }
  }

  const int firstPseudoLayer(this->GetPandora().GetPlugins()->GetPseudoLayerPlugin()->GetPseudoLayerAtIp());

  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "FirstPseudoLayer", firstPseudoLayer));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "EnergyOfPfo", pfoEnergy));
  PANDORA_MONITORING_API(
      SetTreeVariable(this->GetPandora(), m_trainingTreeName, "RawEnergyOfCluster", rawEnergyOfCluster));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "HitEnergies", &hitEnergies));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "CellSize0", &cellSize0));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "CellSize1", &cellSize1));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "CellThickness", &cellThickness));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "HitType", &hitType));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "PseudoLayer", &pseudoLayers));
  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "IsIsolated", &isolatedHit));
  PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));

#endif
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrainingSoftwareCompensation::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                           XmlHelper::ReadValue(xmlHandle, "MyRootFileName", m_myRootFileName));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "SoftCompTrainingTreeName", m_trainingTreeName));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
