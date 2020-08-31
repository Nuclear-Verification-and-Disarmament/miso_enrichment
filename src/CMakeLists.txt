### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###
USE_CYCLUS("misoenrichment" "miso_enrich")
USE_CYCLUS("misoenrichment" "enrichment_calculator")
USE_CYCLUS("misoenrichment" "miso_helper")
USE_CYCLUS("misoenrichment" "flexible_input")

INSTALL_CYCLUS_MODULE("misoenrichment" "./")

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
