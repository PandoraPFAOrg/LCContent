/**
 *  @file   larpandoracontent/LArContent.cc
 * 
 *  @brief  Factory implementations for content intended for use with particle flow reconstruction at an e+e- linear collider
 * 
 *  $Log: $
 */

#include "Pandora/Algorithm.h"
#include "Pandora/Pandora.h"

#include "LCCheating/CheatingClusterCleaningAlgorithm.h"
#include "LCCheating/CheatingParticleIDAlgorithm.h"
#include "LCCheating/CheatingTrackToClusterMatching.h"
#include "LCCheating/PerfectClusteringAlgorithm.h"
#include "LCCheating/PerfectFragmentRemovalAlgorithm.h"
#include "LCCheating/PerfectParticleFlowAlgorithm.h"

#include "LCClustering/ClusteringParentAlgorithm.h"
#include "LCClustering/ConeClusteringAlgorithm.h"
#include "LCClustering/ForcedClusteringAlgorithm.h"

#include "LCFragmentRemoval/BeamHaloMuonRemovalAlgorithm.h"
#include "LCFragmentRemoval/MainFragmentRemovalAlgorithm.h"
#include "LCFragmentRemoval/MergeSplitPhotonsAlgorithm.h"
#include "LCFragmentRemoval/NeutralFragmentRemovalAlgorithm.h"
#include "LCFragmentRemoval/PhotonFragmentMergingAlgorithm.h"
#include "LCFragmentRemoval/PhotonFragmentRemovalAlgorithm.h"
#include "LCFragmentRemoval/RecoPhotonFragmentMergingAlgorithm.h"

#include "LCMonitoring/ClusterComparisonAlgorithm.h"
#include "LCMonitoring/DumpPfosMonitoringAlgorithm.h"
#include "LCMonitoring/EfficiencyMonitoringAlgorithm.h"
#include "LCMonitoring/VisualMonitoringAlgorithm.h"

#include "LCParticleId/FinalParticleIdAlgorithm.h"
#include "LCParticleId/MuonReconstructionAlgorithm.h"
#include "LCParticleId/PhotonReconstructionAlgorithm.h"
#include "LCParticleId/PhotonRecoveryAlgorithm.h"
#include "LCParticleId/PhotonSplittingAlgorithm.h"

#include "LCPersistency/LCEventReadingAlgorithm.h"
#include "LCPersistency/LCEventWritingAlgorithm.h"

#include "LCPfoConstruction/CLICPfoSelectionAlgorithm.h"
#include "LCPfoConstruction/PfoCreationAlgorithm.h"
#include "LCPfoConstruction/PfoCreationParentAlgorithm.h"
#include "LCPfoConstruction/V0PfoCreationAlgorithm.h"

#include "LCPlugins/LCBFieldPlugin.h"
#include "LCPlugins/LCEnergyCorrectionPlugins.h"
#include "LCPlugins/LCParticleIdPlugins.h"
#include "LCPlugins/LCPseudoLayerPlugin.h"
#include "LCPlugins/ALLEGROPseudoLayerPlugin.h"
#include "LCPlugins/LCShowerProfilePlugin.h"
#include "LCPlugins/LCSoftwareCompensation.h"

#include "LCReclustering/ExitingTrackAlg.h"
#include "LCReclustering/ForceSplitTrackAssociationsAlg.h"
#include "LCReclustering/ResolveTrackAssociationsAlg.h"
#include "LCReclustering/SplitMergedClustersAlg.h"
#include "LCReclustering/SplitTrackAssociationsAlg.h"
#include "LCReclustering/TrackDrivenAssociationAlg.h"
#include "LCReclustering/TrackDrivenMergingAlg.h"

#include "LCTopologicalAssociation/BackscatteredTracksAlgorithm.h"
#include "LCTopologicalAssociation/BackscatteredTracks2Algorithm.h"
#include "LCTopologicalAssociation/BrokenTracksAlgorithm.h"
#include "LCTopologicalAssociation/ConeBasedMergingAlgorithm.h"
#include "LCTopologicalAssociation/HighEnergyPhotonRecoveryAlgorithm.h"
#include "LCTopologicalAssociation/IsolatedHitMergingAlgorithm.h"
#include "LCTopologicalAssociation/LoopingTracksAlgorithm.h"
#include "LCTopologicalAssociation/MipPhotonSeparationAlgorithm.h"
#include "LCTopologicalAssociation/MuonPhotonSeparationAlgorithm.h"
#include "LCTopologicalAssociation/MuonClusterAssociationAlgorithm.h"
#include "LCTopologicalAssociation/ProximityBasedMergingAlgorithm.h"
#include "LCTopologicalAssociation/ShowerMipMergingAlgorithm.h"
#include "LCTopologicalAssociation/ShowerMipMerging2Algorithm.h"
#include "LCTopologicalAssociation/ShowerMipMerging3Algorithm.h"
#include "LCTopologicalAssociation/ShowerMipMerging4Algorithm.h"
#include "LCTopologicalAssociation/SoftClusterMergingAlgorithm.h"
#include "LCTopologicalAssociation/TopologicalAssociationParentAlgorithm.h"

#include "LCTrackClusterAssociation/LoopingTrackAssociationAlgorithm.h"
#include "LCTrackClusterAssociation/TrackClusterAssociationAlgorithm.h"
#include "LCTrackClusterAssociation/TrackRecoveryAlgorithm.h"
#include "LCTrackClusterAssociation/TrackRecoveryHelixAlgorithm.h"
#include "LCTrackClusterAssociation/TrackRecoveryInteractionsAlgorithm.h"

#include "LCUtility/CaloHitPreparationAlgorithm.h"
#include "LCUtility/ClusterPreparationAlgorithm.h"
#include "LCUtility/EventPreparationAlgorithm.h"
#include "LCUtility/PfoPreparationAlgorithm.h"
#include "LCUtility/TrackPreparationAlgorithm.h"
#include "LCUtility/TrainingSoftwareCompensation.h"

#include "LCContent.h"

#define LC_ALGORITHM_LIST(d)                                                                                                    \
    d("CheatingClusterCleaning",                CheatingClusterCleaningAlgorithm)                                               \
    d("CheatingParticleID",                     CheatingParticleIDAlgorithm)                                                    \
    d("CheatingTrackToClusterMatching",         CheatingTrackToClusterMatching)                                                 \
    d("PerfectClustering",                      PerfectClusteringAlgorithm)                                                     \
    d("PerfectFragmentRemoval",                 PerfectFragmentRemovalAlgorithm)                                                \
    d("PerfectParticleFlow",                    PerfectParticleFlowAlgorithm)                                                   \
    d("ClusteringParent",                       ClusteringParentAlgorithm)                                                      \
    d("ConeClustering",                         ConeClusteringAlgorithm)                                                        \
    d("ForcedClustering",                       ForcedClusteringAlgorithm)                                                      \
    d("BeamHaloMuonRemoval",                    BeamHaloMuonRemovalAlgorithm)                                                   \
    d("MainFragmentRemoval",                    MainFragmentRemovalAlgorithm)                                                   \
    d("MergeSplitPhotons",                      MergeSplitPhotonsAlgorithm)                                                     \
    d("NeutralFragmentRemoval",                 NeutralFragmentRemovalAlgorithm)                                                \
    d("PhotonFragmentMerging",                  PhotonFragmentMergingAlgorithm)                                                 \
    d("PhotonFragmentRemoval",                  PhotonFragmentRemovalAlgorithm)                                                 \
    d("RecoPhotonFragmentMerging",              RecoPhotonFragmentMergingAlgorithm)                                             \
    d("ClusterComparison",                      ClusterComparisonAlgorithm)                                                     \
    d("DumpPfosMonitoring",                     DumpPfosMonitoringAlgorithm)                                                    \
    d("EfficiencyMonitoring",                   EfficiencyMonitoringAlgorithm)                                                  \
    d("VisualMonitoring",                       VisualMonitoringAlgorithm)                                                      \
    d("FinalParticleId",                        FinalParticleIdAlgorithm)                                                       \
    d("MuonReconstruction",                     MuonReconstructionAlgorithm)                                                    \
    d("PhotonReconstruction",                   PhotonReconstructionAlgorithm)                                                  \
    d("PhotonRecovery",                         PhotonRecoveryAlgorithm)                                                        \
    d("PhotonSplitting",                        PhotonSplittingAlgorithm)                                                       \
    d("LCEventReading",                         LCEventReadingAlgorithm)                                                        \
    d("LCEventWriting",                         LCEventWritingAlgorithm)                                                        \
    d("CLICPfoSelection",                       CLICPfoSelectionAlgorithm)                                                      \
    d("PfoCreation",                            PfoCreationAlgorithm)                                                           \
    d("PfoCreationParent",                      PfoCreationParentAlgorithm)                                                     \
    d("V0PfoCreation",                          V0PfoCreationAlgorithm)                                                         \
    d("ExitingTrack",                           ExitingTrackAlg)                                                                \
    d("ForceSplitTrackAssociations",            ForceSplitTrackAssociationsAlg)                                                 \
    d("ResolveTrackAssociations",               ResolveTrackAssociationsAlg)                                                    \
    d("SplitMergedClusters",                    SplitMergedClustersAlg)                                                         \
    d("SplitTrackAssociations",                 SplitTrackAssociationsAlg)                                                      \
    d("TrackDrivenAssociation",                 TrackDrivenAssociationAlg)                                                      \
    d("TrackDrivenMerging",                     TrackDrivenMergingAlg)                                                          \
    d("BackscatteredTracks",                    BackscatteredTracksAlgorithm)                                                   \
    d("BackscatteredTracks2",                   BackscatteredTracks2Algorithm)                                                  \
    d("BrokenTracks",                           BrokenTracksAlgorithm)                                                          \
    d("ConeBasedMerging",                       ConeBasedMergingAlgorithm)                                                      \
    d("HighEnergyPhotonRecovery",               HighEnergyPhotonRecoveryAlgorithm)                                              \
    d("IsolatedHitMerging",                     IsolatedHitMergingAlgorithm)                                                    \
    d("LoopingTracks",                          LoopingTracksAlgorithm)                                                         \
    d("MipPhotonSeparation",                    MipPhotonSeparationAlgorithm)                                                   \
    d("MuonPhotonSeparation",                   MuonPhotonSeparationAlgorithm)                                                  \
    d("MuonClusterAssociation",                 MuonClusterAssociationAlgorithm)                                                \
    d("ProximityBasedMerging",                  ProximityBasedMergingAlgorithm)                                                 \
    d("ShowerMipMerging",                       ShowerMipMergingAlgorithm)                                                      \
    d("ShowerMipMerging2",                      ShowerMipMerging2Algorithm)                                                     \
    d("ShowerMipMerging3",                      ShowerMipMerging3Algorithm)                                                     \
    d("ShowerMipMerging4",                      ShowerMipMerging4Algorithm)                                                     \
    d("SoftClusterMerging",                     SoftClusterMergingAlgorithm)                                                    \
    d("TopologicalAssociationParent",           TopologicalAssociationParentAlgorithm)                                          \
    d("LoopingTrackAssociation",                LoopingTrackAssociationAlgorithm)                                               \
    d("TrackRecovery",                          TrackRecoveryAlgorithm)                                                         \
    d("TrackRecoveryHelix",                     TrackRecoveryHelixAlgorithm)                                                    \
    d("TrackRecoveryInteractions",              TrackRecoveryInteractionsAlgorithm)                                             \
    d("TrackClusterAssociation",                TrackClusterAssociationAlgorithm)                                               \
    d("CaloHitPreparation",                     CaloHitPreparationAlgorithm)                                                    \
    d("ClusterPreparation",                     ClusterPreparationAlgorithm)                                                    \
    d("EventPreparation",                       EventPreparationAlgorithm)                                                      \
    d("PfoPreparation",                         PfoPreparationAlgorithm)                                                        \
    d("TrackPreparation",                       TrackPreparationAlgorithm)                                                      \
    d("TrainingSoftwareCompensation",           TrainingSoftwareCompensation)

#define LC_ENERGY_CORRECTION_LIST(d)                                                                                            \
    d("CleanClusters",          pandora::HADRONIC,      lc_content::LCEnergyCorrectionPlugins::CleanCluster)                    \
    d("ScaleHotHadrons",        pandora::HADRONIC,      lc_content::LCEnergyCorrectionPlugins::ScaleHotHadrons)                 \
    d("MuonCoilCorrection",     pandora::HADRONIC,      lc_content::LCEnergyCorrectionPlugins::MuonCoilCorrection)

#define LC_PARTICLE_ID_LIST(d)                                                                                                  \
    d("LCEmShowerId",                           lc_content::LCParticleIdPlugins::LCEmShowerId)                                  \
    d("LCPhotonId",                             lc_content::LCParticleIdPlugins::LCPhotonId)                                    \
    d("LCElectronId",                           lc_content::LCParticleIdPlugins::LCElectronId)                                  \
    d("LCMuonId",                               lc_content::LCParticleIdPlugins::LCMuonId)

#define FACTORY Factory

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lc_content
{

#define LC_CONTENT_CREATE_ALGORITHM_FACTORY(a, b)                                                                               \
class b##FACTORY : public pandora::AlgorithmFactory                                                                             \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::Algorithm *CreateAlgorithm() const {return new b;};                                                                \
};

LC_ALGORITHM_LIST(LC_CONTENT_CREATE_ALGORITHM_FACTORY)

} // namespace lc_content

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

#define LC_CONTENT_REGISTER_ALGORITHM(a, b)                                                                                     \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmFactory(pandora, a, new lc_content::b##FACTORY));         \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LCContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LC_ALGORITHM_LIST(LC_CONTENT_REGISTER_ALGORITHM);
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCContent::RegisterBasicPlugins(const pandora::Pandora &pandora, const std::string detectorName)
{
    LC_ENERGY_CORRECTION_LIST(PANDORA_REGISTER_ENERGY_CORRECTION);
    LC_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);

    if (detectorName == "ALLEGRO")
    {
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, new allegro_content::ALLEGROPseudoLayerPlugin));
    }
    else {
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, new lc_content::LCPseudoLayerPlugin));
    }
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetShowerProfilePlugin(pandora, new lc_content::LCShowerProfilePlugin));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCContent::RegisterBFieldPlugin(const pandora::Pandora &pandora, const float innerBField,
    const float muonBarrelBField, const float muonEndCapBField)
{
    return PandoraApi::SetBFieldPlugin(pandora, new lc_content::LCBFieldPlugin(innerBField, muonBarrelBField, muonEndCapBField));
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCContent::RegisterNonLinearityEnergyCorrection(const pandora::Pandora &pandora, const std::string &name,
    const pandora::EnergyCorrectionType energyCorrectionType, const pandora::FloatVector &inputEnergyCorrectionPoints,
    const pandora::FloatVector &outputEnergyCorrectionPoints)
{
    return PandoraApi::RegisterEnergyCorrectionPlugin(pandora, name, energyCorrectionType,
        new lc_content::LCEnergyCorrectionPlugins::NonLinearityCorrection(inputEnergyCorrectionPoints, outputEnergyCorrectionPoints));
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LCContent::RegisterSoftwareCompensationEnergyCorrection(const pandora::Pandora &pandora, const std::string &name,
    const lc_content::LCSoftwareCompensationParameters &parameters)
{
    return PandoraApi::RegisterEnergyCorrectionPlugin(pandora, name, pandora::HADRONIC,
        new lc_content::LCSoftwareCompensation(parameters));
}
