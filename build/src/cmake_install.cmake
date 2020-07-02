# Install script for directory: /Users/test/Uni/Masterarbeit/multi_isotope_enrichment/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/test/.local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmultiisotopeenrichmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cyclus/." TYPE SHARED_LIBRARY FILES "/Users/test/Uni/Masterarbeit/multi_isotope_enrichment/build/lib/cyclus/multiisotopeenrichment/libmultiisotopeenrichment.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cyclus/./libmultiisotopeenrichment.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cyclus/./libmultiisotopeenrichment.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "/Users/test/.local/lib/cyclus/.//libmultiisotopeenrichment.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cyclus/./libmultiisotopeenrichment.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/test/.local/lib:/Users/test/.local/lib/cyclus:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/cyclus"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cyclus/./libmultiisotopeenrichment.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/test/.local/lib:/Users/test/.local/lib/cyclus:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/cyclus"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cyclus/./libmultiisotopeenrichment.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cyclus/./libmultiisotopeenrichment.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmultiisotopeenrichmentx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmultiisotopeenrichment_testingx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/test/Uni/Masterarbeit/multi_isotope_enrichment/build/bin/multiisotopeenrichment_unit_tests")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/multiisotopeenrichment_unit_tests" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/multiisotopeenrichment_unit_tests")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -change "/Users/test/Uni/Masterarbeit/multi_isotope_enrichment/build/lib/cyclus/multiisotopeenrichment/libmultiisotopeenrichment.dylib" "/Users/test/.local/lib/cyclus/.//libmultiisotopeenrichment.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/multiisotopeenrichment_unit_tests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/test/.local/lib:/Users/test/.local/lib/cyclus:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/cyclus"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/multiisotopeenrichment_unit_tests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/test/.local/lib:/Users/test/.local/lib/cyclus:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/cyclus"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/multiisotopeenrichment_unit_tests")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/multiisotopeenrichment_unit_tests")
    endif()
  endif()
endif()

