/**
 *  @file   LCContent/test/LCContentSmokeTest.cc
 * 
 *  @brief  Implementation of a smoke test, registering the full LC content with a pandora instance
 * 
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"

#include "LCContent.h"

#include <iostream>
#include <string>

#define RETURN_ON_FAILURE(operation)                                                            \
{                                                                                               \
    const pandora::StatusCode statusCode(operation);                                            \
                                                                                                \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                             \
    {                                                                                           \
        std::cerr << "LCContentSmokeTest: " << #operation << " returned "                       \
                  << pandora::StatusCodeToString(statusCode) << std::endl;                      \
        return 1;                                                                               \
    }                                                                                           \
}

int main()
{
    try
    {
        const pandora::Pandora pandoraInstance;

        RETURN_ON_FAILURE(LCContent::RegisterAlgorithms(pandoraInstance));
        RETURN_ON_FAILURE(LCContent::RegisterBasicPlugins(pandoraInstance));
        RETURN_ON_FAILURE(LCContent::RegisterBFieldPlugin(pandoraInstance, 3.5f, 0.f, 0.f));
        RETURN_ON_FAILURE(LCContent::RegisterNonLinearityEnergyCorrection(pandoraInstance, "NonLinearity",
            pandora::HADRONIC, pandora::FloatVector(1, 1.f), pandora::FloatVector(1, 1.f)));
    }
    catch (const pandora::StatusCodeException &statusCodeException)
    {
        std::cerr << "LCContentSmokeTest: pandora exception, " << statusCodeException.ToString() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "LCContentSmokeTest: unknown exception" << std::endl;
        return 1;
    }

    std::cout << "LCContentSmokeTest: all LC content registered successfully" << std::endl;

    return 0;
}
