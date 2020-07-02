### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###

#USE_CYCLUS("multiisotopeenrichment" "multi_isotope_enrich")
USE_CYCLUS("multiisotopeenrichment" "enrichment_calculator")
USE_CYCLUS("multiisotopeenrichment" "multi_isotope_helper")

INSTALL_CYCLUS_MODULE("multiisotopeenrichment" "./")

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
