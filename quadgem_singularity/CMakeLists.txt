# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(quadgem)

find_package(ROOT REQUIRED)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(quadgem quadgem.C)
target_link_libraries(quadgem Garfield::Garfield ROOT::Tree ROOT::Physics)


add_executable(rdo_iter rdo_iter.C)
target_link_libraries(rdo_iter Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(quadgem_trackheed quadgem_trackheed.C)
target_link_libraries(quadgem_trackheed Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(extrackheed pion_ms.C)
target_link_libraries(extrackheed Garfield::Garfield ROOT::Tree ROOT::Physics)