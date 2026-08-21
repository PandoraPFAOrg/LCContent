/**
 *  @file   LCContent/src/LCPlugins/ALLEGROPseudoLayerPlugin.cc
 *
 *  @brief  Implementation of the ALLEGRO pseudo layer plugin class.
 *
 *  @author Giovanni Marchiori
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCPlugins/ALLEGROPseudoLayerPlugin.h"

using namespace pandora;

namespace allegro_content
{

ALLEGROPseudoLayerPlugin::ALLEGROPseudoLayerPlugin() :
    m_ecalBarrelOuterZ(0.),
    m_ecalBarrelNLayers(0),
    m_ecalEndCapInnerZ(0.),
    m_ecalEndCapInnerR(0.),
    m_ecalEndCapOuterR(0.),
    m_ecalEndCapNLayers(0),
    m_hcalBarrelOuterZ(0.),
    m_hcalBarrelLayerRadii(0.),
    m_hcalBarrelNLayers(0),
    m_hcalEndCapInnerR(0.),
    m_hcalEndCapOuterR(0.),
    m_hcalEndCapNLayers(0),
    m_muonBarrelOuterZ(0.),
    m_muonBarrelNLayers(0),
    m_muonEndCapInnerR(0.),
    m_muonEndCapOuterR(0.),
    m_muonEndCapNLayers(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ALLEGROPseudoLayerPlugin::Initialize()
{
    try
    {
        this->StoreLayerPositions();
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "ALLEGROPseudoLayerPlugin: Incomplete geometry - consider using a different PseudoLayerCalculator." << std::endl;
        return statusCodeException.GetStatusCode();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int ALLEGROPseudoLayerPlugin::GetPseudoLayer(const CartesianVector &positionVector) const
{
    const float tolerance = 5e-2; // 50 um

    // retrieve |z| and check that it is lower than the outer z of the detector (from the last muon endcap layer)
    const float z(std::fabs(positionVector.GetZ()));

    if (z > (m_muonEndCapLayerZ[m_muonEndCapNLayers-1] + tolerance))
    {
        std::cout << "Hit is beyond muon endcap z" << std::endl;
        std::cout << "Hit coordinates are: " << positionVector << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    // retrieve |r| and check that it is not outside detector envelope
    const float x(std::fabs(positionVector.GetX()));
    const float y(std::fabs(positionVector.GetY()));
    const float r(sqrt(x*x + y*y));

    if (r > (m_muonBarrelLayerRadii[m_muonBarrelNLayers-1] + tolerance) && r > m_muonEndCapOuterR)
    {
        std::cout << "Hit is beyond muon outer r" << std::endl;
        std::cout << "Hit coordinates are: " << positionVector << " r: " << r << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    if (r < (m_ecalBarrelLayerRadii[0] - tolerance) && z < m_ecalBarrelOuterZ)
    {
        // Hit is upstream of calo volume, set layer to 0 (could be track extrapolation)
        // std::cout << "Hit is upstream of calorimeter" << std::endl;
        // std::cout << "Hit coordinates are: " << positionVector << " r: " << r << std::endl;
        return 0;
    }

    // Reserve a pseudo layer for track projections, etc.
    unsigned int pseudoLayer(1);

    // GM debug code to be removed
    /*
    std::cout << " z : " << z << std::endl;
    std::cout << " m_ecalBarrelOuterZ : " << m_ecalBarrelOuterZ << std::endl;
    std::cout << " r : " << r << std::endl;
    std::cout << " m_ecalBarrelLayerRadii[0] " << m_ecalBarrelLayerRadii[0] << std::endl;
    std::cout << " m_ecalBarrelLayerRadii[m_ecalBarrelNLayers-1] " << m_ecalBarrelLayerRadii[m_ecalBarrelNLayers-1] << std::endl;
    */
    // ecal barrel
    // if (z < m_ecalBarrelOuterZ && r > (m_ecalBarrelLayerRadii[0]-tolerance) && r < (m_ecalBarrelLayerRadii[m_ecalBarrelNLayers-1]+tolerance))
    // using this instead because the cells at the edge of the barrel are positioned outside of the barrel volume
    // and thus fail the z < m_ecalBarrelOuterZ check..
    if (z < m_ecalEndCapInnerZ && r > (m_ecalBarrelLayerRadii[0]-tolerance) && r < (m_ecalBarrelLayerRadii[m_ecalBarrelNLayers-1]+tolerance))
    {
        unsigned int iLayer(0);
        StatusCode statusCode = FindMatchingLayer(r, m_ecalBarrelLayerRadii, iLayer);
        if (STATUS_CODE_SUCCESS != statusCode)
        {
            std::cout << "Layer for ecal hit not found" << std::endl;
            std::cout << "Hit coordinates are: " << positionVector << " r: " << r << std::endl;
            throw StatusCodeException(statusCode);
        }
        pseudoLayer += iLayer;
        return pseudoLayer;
    }

    // hcal barrel
    if (z < m_hcalBarrelOuterZ && r > (m_hcalBarrelLayerRadii[0]-tolerance) && r < (m_hcalBarrelLayerRadii[m_hcalBarrelNLayers-1]+tolerance))
    {
        pseudoLayer += m_ecalBarrelNLayers;
        unsigned int iLayer(0);
        StatusCode statusCode = FindMatchingLayer(r, m_hcalBarrelLayerRadii, iLayer);
        if (STATUS_CODE_SUCCESS != statusCode)
        {
            std::cout << "Layer for hcal hit not found" << std::endl;
            std::cout << "Hit coordinates are: " << positionVector << " r: " << r << std::endl;
            throw StatusCodeException(statusCode);
        }
        pseudoLayer += iLayer;
        return pseudoLayer;
    }

    // muon barrel
    if (z < m_muonBarrelOuterZ && r > (m_muonBarrelLayerRadii[0]-tolerance) && r < (m_muonBarrelLayerRadii[m_muonBarrelNLayers-1]+tolerance))
    {
        pseudoLayer += m_ecalBarrelNLayers;
        pseudoLayer += m_hcalBarrelNLayers;
        unsigned int iLayer(0);
        StatusCode statusCode = FindMatchingLayer(r, m_muonBarrelLayerRadii, iLayer);
        if (STATUS_CODE_SUCCESS != statusCode)
        {
            std::cout << "Layer for muon hit not found" << std::endl;
            std::cout << "Hit coordinates are: " << positionVector << " r: " << r << std::endl;
            throw StatusCodeException(statusCode);
        }
        pseudoLayer += iLayer;
        return pseudoLayer;
    }

    // endcaps to be implemented

    // if layer not found
    std::cout << "Layer not found" << std::endl;
    std::cout << "Hit coordinates are: " << positionVector << " r: " << r << std::endl;
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ALLEGROPseudoLayerPlugin::FindMatchingLayer(const float position, const LayerPositionList &layerPositionList,
                                                       unsigned int &layer) const
{
    // GM debug code to be removed
    // std::cout << "FindMatchingLayer: position = " << position << std::endl;
    const float tolerance = 5e-2; // 50 um
    for (int i=0; i<layerPositionList.size(); i++) {
        // GM debug code to be removed
        // std::cout << "FindMatchingLayer: layer = " << i << " position = " << layerPositionList[i] << std::endl;
        if (std::fabs(position - layerPositionList[i]) < tolerance) {
            layer = i;
            return STATUS_CODE_SUCCESS;
        }
    }

    layer = -1;
    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ALLEGROPseudoLayerPlugin::StoreLayerPositions()
{
    // GM: we might want to revisit this a bit due to the HCAL endcap 3-part geometry
    // in case we split the HCAL endcap in 3 different subdetectors (at least logically, for pandora)
    // we would then retrieve the subdetectors by name rather than type, unless we use something different
    // than HCAL_ENDCAP for the other two wheels of the EC (like plug and ring?)

    // also for the ECAL endcap we will have to deal properly with the layer information..

    // retrieve the layer positions of the various subdetectors
    // GM: check if we get both min and max...
    const GeometryManager *const pGeometryManager(this->GetPandora().GetGeometry());
    this->StoreLayerPositions(pGeometryManager->GetSubDetector(ECAL_BARREL), m_ecalBarrelLayerRadii);
    this->StoreLayerPositions(pGeometryManager->GetSubDetector(HCAL_BARREL), m_hcalBarrelLayerRadii);
    this->StoreLayerPositions(pGeometryManager->GetSubDetector(MUON_BARREL), m_muonBarrelLayerRadii);
    this->StoreLayerPositions(pGeometryManager->GetSubDetector(ECAL_ENDCAP), m_ecalEndCapLayerZ);
    this->StoreLayerPositions(pGeometryManager->GetSubDetector(HCAL_ENDCAP), m_hcalEndCapLayerZ);
    this->StoreLayerPositions(pGeometryManager->GetSubDetector(MUON_ENDCAP), m_muonEndCapLayerZ);
    m_ecalBarrelNLayers = m_ecalBarrelLayerRadii.size();
    m_hcalBarrelNLayers = m_hcalBarrelLayerRadii.size();
    m_muonBarrelNLayers = m_muonBarrelLayerRadii.size();
    m_ecalEndCapNLayers = m_ecalEndCapLayerZ.size();
    m_hcalEndCapNLayers = m_hcalEndCapLayerZ.size();
    m_muonEndCapNLayers = m_muonEndCapLayerZ.size();

    // retrieve the needed envelope info of the various subdetectors
    m_ecalBarrelOuterZ = pGeometryManager->GetSubDetector(ECAL_BARREL).GetOuterZCoordinate();
    m_ecalEndCapInnerZ = pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetInnerZCoordinate();
    m_ecalEndCapInnerR = pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetInnerRCoordinate();
    m_ecalEndCapOuterR = pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetOuterRCoordinate();
    m_hcalBarrelOuterZ = pGeometryManager->GetSubDetector(HCAL_BARREL).GetOuterZCoordinate();
    m_hcalEndCapInnerR = pGeometryManager->GetSubDetector(HCAL_ENDCAP).GetInnerRCoordinate();
    m_hcalEndCapOuterR = pGeometryManager->GetSubDetector(HCAL_ENDCAP).GetOuterRCoordinate();
    m_muonBarrelOuterZ = pGeometryManager->GetSubDetector(MUON_BARREL).GetOuterZCoordinate();
    m_muonEndCapInnerR = pGeometryManager->GetSubDetector(MUON_ENDCAP).GetInnerRCoordinate();
    m_muonEndCapOuterR = pGeometryManager->GetSubDetector(MUON_ENDCAP).GetOuterRCoordinate();

    // print retrieved information
    Print();

    // check that layer information is filled
    if (m_ecalBarrelLayerRadii.empty() ||
        m_hcalBarrelLayerRadii.empty() ||
        m_muonBarrelLayerRadii.empty() ||
        // GM: ecal endcap layer information needs to be filled
        // m_ecalEndCapLayerZ.empty() ||
        m_hcalEndCapLayerZ.empty() ||
        m_muonEndCapLayerZ.empty() ||
        m_ecalBarrelOuterZ == 0. ||
        m_hcalBarrelOuterZ == 0. ||
        m_muonBarrelOuterZ == 0. ||
        m_ecalEndCapInnerZ == 0. ||
        m_ecalEndCapInnerR == 0. ||
        m_hcalEndCapInnerR == 0. ||
        m_muonEndCapInnerR == 0. ||
        m_ecalEndCapOuterR == 0. ||
        m_hcalEndCapOuterR == 0. ||
        m_muonEndCapOuterR == 0.
        )
    {
        std::cout << "ALLEGROPseudoLayerPlugin: some position information are not specified." << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ALLEGROPseudoLayerPlugin::StoreLayerPositions(const SubDetector &subDetector, LayerPositionList &layerPositionList)
{
    if (!subDetector.IsMirroredInZ())
    {
        std::cout << "ALLEGROPseudoLayerPlugin: Error, detector must be symmetrical about z=0 plane." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const SubDetector::SubDetectorLayerVector &subDetectorLayerVector(subDetector.GetSubDetectorLayerVector());

    for (SubDetector::SubDetectorLayerVector::const_iterator iter = subDetectorLayerVector.begin(), iterEnd = subDetectorLayerVector.end(); iter != iterEnd; ++iter)
    {
        layerPositionList.push_back(iter->GetClosestDistanceToIp());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
void ALLEGROPseudoLayerPlugin::StoreDetectorOuterEdge()
{
    const GeometryManager *const pGeometryManager(this->GetPandora().GetGeometry());

    m_barrelEdgeR = (std::max(pGeometryManager->GetSubDetector(ECAL_BARREL).GetOuterRCoordinate(), std::max(
        pGeometryManager->GetSubDetector(HCAL_BARREL).GetOuterRCoordinate(),
        pGeometryManager->GetSubDetector(MUON_BARREL).GetOuterRCoordinate()) ));

    m_endCapEdgeZ = (std::max(std::fabs(pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetOuterZCoordinate()), std::max(
        std::fabs(pGeometryManager->GetSubDetector(HCAL_ENDCAP).GetOuterZCoordinate()),
        std::fabs(pGeometryManager->GetSubDetector(MUON_ENDCAP).GetOuterZCoordinate())) ));

    if ((m_barrelLayerPositions.end() != std::upper_bound(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end(), m_barrelEdgeR)) )
    // FIXME: AD: why this line is checking if all layers inner radii are lower than endcap outer z coordinate??? ALLEGRO detector can not satisfy this -> disabled
    // || (m_endCapLayerPositions.end() != std::upper_bound(m_endCapLayerPositions.begin(), m_endCapLayerPositions.end(), m_endCapEdgeZ)))
    {
        std::cout << "ALLEGROPseudoLayerPlugin: Layers specified outside detector edge." << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    m_barrelLayerPositions.push_back(m_barrelEdgeR);
    m_endCapLayerPositions.push_back(m_endCapEdgeZ);
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

void ALLEGROPseudoLayerPlugin::Print()
{
    std::cout << "ALLEGROPseudoLayerPlugin: information read from geometry" << std::endl;
    std::cout << "ECAL barrel" << std::endl;
    std::cout << "- outer z: " << m_ecalBarrelOuterZ << std::endl;
    std::cout << "- layers: " << m_ecalBarrelNLayers << std::endl;
    std::cout << "- layer radii: " << std::endl;
    for (int i=0; i<m_ecalBarrelNLayers; i++)
    {
        std::cout << "    " << i << " : " << m_ecalBarrelLayerRadii[i] << std::endl;
    }
    std::cout << "ECAL endcap" << std::endl;
    std::cout << "- inner r: " << m_ecalEndCapInnerR << std::endl;
    std::cout << "- outer r: " << m_ecalEndCapOuterR << std::endl;
    std::cout << "- layers " << m_ecalEndCapNLayers << std::endl;
    std::cout << "- layer |z|: " << std::endl;
    for (int i=0; i<m_ecalEndCapNLayers; i++)
    {
        std::cout << "    " << i << " : " << m_ecalEndCapLayerZ[i] << std::endl;
    }
    std::cout << "HCAL barrel" << std::endl;
    std::cout << "- outer z: " << m_hcalBarrelOuterZ << std::endl;
    std::cout << "- layers: " << m_hcalBarrelNLayers << std::endl;
    std::cout << "- layer radii: " << std::endl;
    for (int i=0; i<m_hcalBarrelNLayers; i++)
    {
        std::cout << "    " << i << " : " << m_hcalBarrelLayerRadii[i] << std::endl;
    }
    std::cout << "HCAL endcap" << std::endl;
    std::cout << "- inner r: " << m_hcalEndCapInnerR << std::endl;
    std::cout << "- outer r: " << m_hcalEndCapOuterR << std::endl;
    std::cout << "- layers " << m_hcalEndCapNLayers << std::endl;
    std::cout << "- layer |z|: " << std::endl;
    for (int i=0; i<m_hcalEndCapNLayers; i++)
    {
        std::cout << "    " << i << " : " << m_hcalEndCapLayerZ[i] << std::endl;
    }
    std::cout << "MUON barrel" << std::endl;
    std::cout << "- outer z: " << m_muonBarrelOuterZ << std::endl;
    std::cout << "- layers: " << m_muonBarrelNLayers << std::endl;
    std::cout << "- layer radii: " << std::endl;
    for (int i=0; i<m_muonBarrelNLayers; i++)
    {
        std::cout << "    " << i << " : " << m_muonBarrelLayerRadii[i] << std::endl;
    }
    std::cout << "MUON endcap" << std::endl;
    std::cout << "- inner r: " << m_muonEndCapInnerR << std::endl;
    std::cout << "- outer r: " << m_muonEndCapOuterR << std::endl;
    std::cout << "- layers " << m_muonEndCapNLayers << std::endl;
    std::cout << "- layer |z|: " << std::endl;
    for (int i=0; i<m_muonEndCapNLayers; i++)
    {
        std::cout << "    " << i << " : " << m_muonEndCapLayerZ[i] << std::endl;
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ALLEGROPseudoLayerPlugin::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace allegro_content
