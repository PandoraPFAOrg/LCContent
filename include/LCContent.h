/**
 *  @file   LCContent/include/LCContent.h
 * 
 *  @brief  Header file detailing content for use with particle flow reconstruction at an e+e- linear collider
 * 
 *  $Log: $
 */
#ifndef LINEAR_COLLIDER_CONTENT_H
#define LINEAR_COLLIDER_CONTENT_H 1

#include "Api/PandoraApi.h"

namespace lc_content {class LCSoftwareCompensationParameters;}

/**
 *  @brief  LCContent class
 */
class LCContent
{
public:
    /**
     *  @brief  Register all the linear collider algorithms with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);

    /**
     *  @brief  Register the basic (no configuration required on user side) linear collider plugins with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  detectorName the name of the detector (only relevant for ALLEGRO)
     */
    static pandora::StatusCode RegisterBasicPlugins(const pandora::Pandora &pandora, const std::string detectorName="");

    /**
     *  @brief  Register the b field plugin (note user side configuration) with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  innerBField the bfield in the main tracker, ecal and hcal, units Tesla
     *  @param  muonBarrelBField the bfield in the muon barrel, units Tesla
     *  @param  muonEndCapBField the bfield in the muon endcap, units Tesla
     */
    static pandora::StatusCode RegisterBFieldPlugin(const pandora::Pandora &pandora, const float innerBField, const float muonBarrelBField,
        const float muonEndCapBField);

    /**
     *  @brief  Register the non linearity energy correction plugin (note user side configuration) with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  name the name/label associated with the energy correction plugin
     *  @param  energyCorrectionType the energy correction type
     *  @param  inputEnergyCorrectionPoints the input energy points for energy correction
     *  @param  outputEnergyCorrectionPoints the output energy points for energy correction
     */
    static pandora::StatusCode RegisterNonLinearityEnergyCorrection(const pandora::Pandora &pandora, const std::string &name,
        const pandora::EnergyCorrectionType energyCorrectionType, const pandora::FloatVector &inputEnergyCorrectionPoints,
        const pandora::FloatVector &outputEnergyCorrectionPoints);
    
    /**
     *  @brief  Register the software compensation energy correction plugin (note user side configuration) with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  name the name/label associated with the energy correction plugin
     *  @param  energyCorrectionType the energy correction type
     *  @param  parameters the software compensation input parameters
     */
    static pandora::StatusCode RegisterSoftwareCompensationEnergyCorrection(const pandora::Pandora &pandora, const std::string &name,
        const lc_content::LCSoftwareCompensationParameters &parameters);
};

#endif // #ifndef LINEAR_COLLIDER_CONTENT_H
