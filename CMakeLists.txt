cmake_minimum_required(VERSION 3.23 FATAL_ERROR)


if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "LCContent requires an out-of-source build.")
endif()

# project name
if(NOT LC_CONTENT_LIBRARY_NAME STREQUAL "LCPandoraContent")
    set(LC_CONTENT_LIBRARY_NAME "LCContent")
endif()
project(${LC_CONTENT_LIBRARY_NAME} LANGUAGES CXX)

set(${PROJECT_NAME}_VERSION_MAJOR 03)
set(${PROJECT_NAME}_VERSION_MINOR 02)
set(${PROJECT_NAME}_VERSION_PATCH 00)
set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

include(GNUInstallDirs)

find_package(PandoraSDK 03.00.00 REQUIRED)

if(PANDORA_MONITORING)
    find_package(PandoraMonitoring 03.00.00 REQUIRED)
endif()

# Optional ONNX ML-inference sub-library.  The ONNX-dependent algorithms live in src/MLInference and
# are built into a SEPARATE library (libLCContentML), so the core LCContent has NO onnxruntime
# dependency and can be linked without ONNX in the stack.  Enable with -DLCContent_BUILD_ML=ON.
option(LCContent_BUILD_ML "Build the ONNX ML-inference algorithms into libLCContentML (needs onnxruntime)" OFF)

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
# C++17 required by the ONNX Runtime C++ API (also matches the build flags used downstream).
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

# Prevent CMake falls back to the latest standard the compiler does support
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wno-long-long -Wno-sign-compare -Wshadow -fno-strict-aliasing ${CMAKE_CXX_FLAGS}")

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
endif()


file(GLOB_RECURSE LC_CONTENT_SRCS RELATIVE ${PROJECT_SOURCE_DIR} "src/*.cc")
file(GLOB_RECURSE headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "include/*.h")

# Keep the ONNX-dependent MLInference sources/headers out of the core library.
list(FILTER LC_CONTENT_SRCS EXCLUDE REGEX "^src/MLInference/")
list(FILTER headers        EXCLUDE REGEX "^include/MLInference/")

add_library(${PROJECT_NAME} SHARED ${LC_CONTENT_SRCS})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
set(LCContent_LIBRARIES ${PROJECT_NAME}::${PROJECT_NAME})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET headers TYPE HEADERS BASE_DIRS "include" FILES ${headers})
if(TARGET PandoraPFA::PandoraSDK)
  target_link_libraries(${PROJECT_NAME} PUBLIC PandoraPFA::PandoraSDK)
else()
  # SYSTEM: the headers of the external dependency are not maintained here, so
  # diagnostics raised inside them must not fail a -Werror build of LCContent
  target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${PandoraSDK_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${PandoraSDK_LIBRARIES})
endif()
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(PANDORA_MONITORING)
  if(TARGET PandoraPFA::PandoraMonitoring)
      target_link_libraries(${PROJECT_NAME} PUBLIC PandoraPFA::PandoraMonitoring)
  else()
      target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${PandoraMonitoring_INCLUDE_DIRS})
      target_link_libraries(${PROJECT_NAME} PUBLIC ${PandoraMonitoring_LIBRARIES})
  endif()
  target_compile_definitions(${PROJECT_NAME} PRIVATE -DMONITORING ${PandoraMonitoring_DEFINITIONS})
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(test)
endif()

# ---- Optional ONNX ML-inference library (libLCContentML) ----
if(LCContent_BUILD_ML)
    # New onnxruntime packages ship onnxruntime-Config.cmake (target onnxruntime::onnxruntime);
    # older ones ship ONNXRuntime.
    find_package(onnxruntime)
    if(NOT onnxruntime_FOUND)
        find_package(ONNXRuntime)
    endif()
    if(NOT (onnxruntime_FOUND OR ONNXRuntime_FOUND))
        message(FATAL_ERROR "LCContent_BUILD_ML=ON but ONNX Runtime was not found.")
    endif()

    file(GLOB_RECURSE LC_CONTENT_ML_SRCS    RELATIVE ${PROJECT_SOURCE_DIR}        "src/MLInference/*.cc")
    file(GLOB_RECURSE LC_CONTENT_ML_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "include/MLInference/*.h")

    add_library(LCContentML SHARED ${LC_CONTENT_ML_SRCS})
    add_library(${PROJECT_NAME}::LCContentML ALIAS LCContentML)
    target_sources(LCContentML PUBLIC FILE_SET headers TYPE HEADERS BASE_DIRS "include" FILES ${LC_CONTENT_ML_HEADERS})
    target_link_libraries(LCContentML PUBLIC ${PROJECT_NAME})   # brings PandoraSDK + LCContent includes
    if(TARGET onnxruntime::onnxruntime)
        target_link_libraries(LCContentML PUBLIC onnxruntime::onnxruntime)
    else()
        target_include_directories(LCContentML PUBLIC ${ONNXRuntime_INCLUDE_DIRS})
        target_link_libraries(LCContentML PUBLIC ${ONNXRuntime_LIBRARIES})
    endif()

    install(TARGETS LCContentML
      EXPORT ${PROJECT_NAME}Targets
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
      FILE_SET headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT dev
    )
endif()

option(LCContent_BUILD_DOCS "Build documentation for ${PROJECT_NAME}" OFF)
if(LCContent_BUILD_DOCS)
    add_subdirectory(doc)
endif()

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/LCContent
  FILE_SET headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT dev
)

include(cmake/LCContentCreateConfig.cmake)

