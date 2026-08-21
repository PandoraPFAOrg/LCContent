/**
 *  @file   LCContent/include/LCPlugins/ALLEGROPseudoLayerPlugin.h
 *
 *  @brief  Header file for the ALLEGRO pseudo layer plugin class.
 *
 *  @author Giovanni Marchiori
 */
#ifndef ALLEGRO_PSEUDO_LAYER_PLUGIN_H
#define ALLEGRO_PSEUDO_LAYER_PLUGIN_H 1

#include "Pandora/PandoraInputTypes.h"

#include "Plugins/PseudoLayerPlugin.h"

namespace allegro_content
{

/**
 *  @brief  ALLEGROPseudoLayerPlugin class
 */
class ALLEGROPseudoLayerPlugin : public pandora::PseudoLayerPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    ALLEGROPseudoLayerPlugin();

private:
    pandora::StatusCode Initialize();
    unsigned int GetPseudoLayer(const pandora::CartesianVector &positionVector) const;
    unsigned int GetPseudoLayerAtIp() const;

    /**
     *  @brief  Get the appropriate pseudolayer for a specified parameters
     *
     *  @param  rCoordinate the radial coordinate
     *  @param  zCoordinate the z coordinate
     *  @param  rCorrection the barrel/endcap overlap r correction
     *  @param  zCorrection the barrel/endcap overlap z correction
     *  @param  barrelInnerR the barrel inner r coordinate
     *  @param  endCapInnerZ the endcap inner z coordinate
     *  @param  pseudoLayer to receive the appropriate pseudolayer
     */
    pandora::StatusCode GetPseudoLayer(const float rCoordinate, const float zCoordinate, const float rCorrection, const float zCorrection,
        const float barrelInnerR, const float endCapInnerZ, unsigned int &pseudoLayer) const;

    typedef std::vector<float> LayerPositionList;

    /**
     *  @brief  Find the layer number corresponding to a specified position, via reference to a specified layer position list
     *
     *  @param  position the specified position
     *  @param  layerPositionList the specified layer position list
     *  @param  layer to receive the layer number
     */
     pandora::StatusCode FindMatchingLayer(const float position, const LayerPositionList &layerPositionList, unsigned int &layer) const;

    /**
     *  @brief  Store all revelevant barrel and endcap layer positions upon initialization
     */
    void StoreLayerPositions();

    /**
     *  @brief  Store subdetector layer positions upon initialization
     *
     *  @param  subDetector the sub detector
     *  @param  layerParametersList the layer parameters list
     */
    void StoreLayerPositions(const pandora::SubDetector &subDetector, LayerPositionList &LayerPositionList);

    /**
     *  @brief  Store positions of barrel and endcap outer edges upon initialization
     */
    //   void StoreDetectorOuterEdge();


    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Print retrieved geometry information
     */
    void Print();

    float m_ecalBarrelOuterZ;                       ///< Maximum z of ECAL barrel
    LayerPositionList m_ecalBarrelLayerRadii;       ///< List of radial positions of the ECAL barrel layers
    int m_ecalBarrelNLayers;                        ///< Number of ECAL barrel layers

    float m_ecalEndCapInnerZ;                       ///< Minimum z of ECAL endcap (added temporarily until we get info about endcap layers)
    float m_ecalEndCapInnerR;                       ///< Minimum r of ECAL endcap
    float m_ecalEndCapOuterR;                       ///< Maximum r of ECAL endcap
    LayerPositionList m_ecalEndCapLayerZ;           ///< List of z positions of the ECAL endcap layers
    int m_ecalEndCapNLayers;                        ///< Number of ECAL endcap layers

    float m_hcalBarrelOuterZ;                       ///< Maximum z of HCAL barrel
    LayerPositionList m_hcalBarrelLayerRadii;       ///< List of radial positions of the HCAL barrel layers
    int m_hcalBarrelNLayers;                        ///< Number of ECAL barrel layers

    // GM: this will have to be revisited later to handle properly the 3-part HCAL endcap
    // we will need vectors of floats for innerR/outerR and vectors of layerpositionlist for the layer depths
    // we might even want to have something a bit more complicated if the first two parts of the endcap
    // have R-layers and the 3rd one has z-layers
    // will do that later
    float m_hcalEndCapInnerR;                       ///< Minimum r of HCAL endcap
    float  m_hcalEndCapOuterR;                      ///< Maximum r of HCAL endcap
    LayerPositionList m_hcalEndCapLayerZ;          ///< List of z positions of the HCAL endcap layers
    int m_hcalEndCapNLayers;                        ///< Number of HCAL endcap layers

    float m_muonBarrelOuterZ;                       ///< Maximum z of MUON barrel
    LayerPositionList m_muonBarrelLayerRadii;       ///< List of radial positions of the MUON barrel layers
    int m_muonBarrelNLayers;                        ///< Number of MUON barrel layers

    float m_muonEndCapInnerR;                       ///< Minimum r of MUON endcap
    float  m_muonEndCapOuterR;                      ///< Maximum r of MUON endcap
    LayerPositionList m_muonEndCapLayerZ;           ///< List of z positions of the ECAL endcap layers
    int m_muonEndCapNLayers;                        ///< Number of MUON endcap layers
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ALLEGROPseudoLayerPlugin::GetPseudoLayerAtIp() const
{
    const unsigned int pseudoLayerAtIp(this->GetPseudoLayer(pandora::CartesianVector(0.f, 0.f, 0.f)));
    return pseudoLayerAtIp;
}

} // namespace allegro_content

#endif // #ifndef ALLEGRO_PSEUDO_LAYER_PLUGIN_H
