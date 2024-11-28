cmake_minimum_required(VERSION 3.2)

project(continuum-multirobot-estimation)

# Force Release mode unless specified otherwise
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

#Important CMAKE stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#VTK
set(VTK_DIR "usr/local/lib")
find_package(VTK REQUIRED)

#Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR})

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(continuum-multirobot-estimation src/main.cpp ${H_FILES} ${CPP_FILES})
target_link_libraries(continuum-multirobot-estimation ${VTK_LIBRARIES} Eigen3::Eigen)
