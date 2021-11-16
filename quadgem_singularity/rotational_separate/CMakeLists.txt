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
add_executable(qg1 quadgem_rota1.C)
target_link_libraries(qg1 Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(qg2 quadgem_rota2.C)
target_link_libraries(qg2 Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(qg3 quadgem_rota3.C)
target_link_libraries(qg3 Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(qg4 quadgem_rota4.C)
target_link_libraries(qg4 Garfield::Garfield ROOT::Tree ROOT::Physics)