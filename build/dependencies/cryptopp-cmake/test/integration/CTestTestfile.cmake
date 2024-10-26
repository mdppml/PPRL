# CMake generated Testfile for 
# Source directory: /Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration
# Build directory: /Users/sesame/Desktop/CECILIA/build/dependencies/cryptopp-cmake/test/integration
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[cryptopp-int-install-default-configure]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "-GUnix Makefiles" "-D" "TEST_CMAKE_FILES_DIR=/Users/sesame/Desktop/CECILIA/cmake" "-D" "INT_TEST_CMAKE_FILES_DIR=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/cmake" "-D" "TEST_EXAMPLE_SOURCES_DIR=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/example-src" "-D" "USE_CCACHE=TRUE" "-D" "CMAKE_USER_MAKE_RULES_OVERRIDE=/Users/sesame/Desktop/CECILIA/cmake/ResetInitialCompilerOptions.cmake" "-D" "CMAKE_INSTALL_PREFIX=/Users/sesame/Desktop/CECILIA/build/test-dirs/install/int-install-default" "-S" "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/int-install-default" "-B" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-default" "-D" "CPM_cryptopp-cmake_SOURCE=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/../.." "-D" "CMAKE_VERBOSE_MAKEFILE=ON" "-D" "CRYPTOPP_MINIMUM_CMAKE_VERSION=3.20" "-D" "CMAKE_BUILD_TYPE=Release")
set_tests_properties([=[cryptopp-int-install-default-configure]=] PROPERTIES  FIXTURES_SETUP "int-install-default-config" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" install_dir "/Users/sesame/Desktop/CECILIA/build/test-dirs/install/int-install-default" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;21;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-default-build]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-default" "--config" "Release")
set_tests_properties([=[cryptopp-int-install-default-build]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-config" FIXTURES_SETUP "int-install-default-build" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;55;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-default-install]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-default" "--target" "install" "--config" "Release")
set_tests_properties([=[cryptopp-int-install-default-install]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-build" FIXTURES_SETUP "int-install-default-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;65;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-default-checks]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-default" "--target" "do-checks" "--config" "Release")
set_tests_properties([=[cryptopp-int-install-default-checks]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;76;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-configure]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "-GUnix Makefiles" "-D" "TEST_CMAKE_FILES_DIR=/Users/sesame/Desktop/CECILIA/cmake" "-D" "INT_TEST_CMAKE_FILES_DIR=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/cmake" "-D" "TEST_EXAMPLE_SOURCES_DIR=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/example-src" "-D" "USE_CCACHE=TRUE" "-D" "CMAKE_USER_MAKE_RULES_OVERRIDE=/Users/sesame/Desktop/CECILIA/cmake/ResetInitialCompilerOptions.cmake" "-D" "CMAKE_INSTALL_PREFIX=/Users/sesame/Desktop/CECILIA/build/test-dirs/install/int-install-prefix" "-S" "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/int-install-prefix" "-B" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-prefix" "-D" "CPM_cryptopp-cmake_SOURCE=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/../.." "-D" "CMAKE_VERBOSE_MAKEFILE=ON" "-D" "CRYPTOPP_MINIMUM_CMAKE_VERSION=3.20" "-D" "CMAKE_BUILD_TYPE=Release")
set_tests_properties([=[cryptopp-int-install-prefix-configure]=] PROPERTIES  FIXTURES_SETUP "int-install-prefix-config" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" install_dir "/Users/sesame/Desktop/CECILIA/build/test-dirs/install/int-install-prefix" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;21;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-build]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-prefix" "--config" "Release")
set_tests_properties([=[cryptopp-int-install-prefix-build]=] PROPERTIES  FIXTURES_REQUIRED "int-install-prefix-config" FIXTURES_SETUP "int-install-prefix-build" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;55;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-install]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-prefix" "--target" "install" "--config" "Release")
set_tests_properties([=[cryptopp-int-install-prefix-install]=] PROPERTIES  FIXTURES_REQUIRED "int-install-prefix-build" FIXTURES_SETUP "int-install-prefix-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;65;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-checks]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-install-prefix" "--target" "do-checks" "--config" "Release")
set_tests_properties([=[cryptopp-int-install-prefix-checks]=] PROPERTIES  FIXTURES_REQUIRED "int-install-prefix-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;76;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-find-package-configure]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "-GUnix Makefiles" "-D" "TEST_EXAMPLE_SOURCES_DIR=/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/example-src" "-D" "CRYPTOPP_CMAKE_INSTALL_ROOT=/Users/sesame/Desktop/CECILIA/build/test-dirs/install" "-S" "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/int-find-package" "-B" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-find-package" "-D" "CMAKE_VERBOSE_MAKEFILE=ON" "-D" "CRYPTOPP_MINIMUM_CMAKE_VERSION=3.20" "-D" "CMAKE_BUILD_TYPE=Release" "-D" "cryptopp_DIR=/Users/sesame/Desktop/CECILIA/build/test-dirs/install/int-install-default/share/cmake/cryptopp")
set_tests_properties([=[cryptopp-int-find-package-configure]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-install" FIXTURES_SETUP "int-find-package-config" LABELS "cryptopp;cryptopp-find_package" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;95;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-find-package-build]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/sesame/Desktop/CECILIA/build/test-dirs/int-find-package" "--config" "Release")
set_tests_properties([=[cryptopp-int-find-package-build]=] PROPERTIES  FIXTURES_REQUIRED "int-find-package-config" LABELS "cryptopp;cryptopp-find_package" _BACKTRACE_TRIPLES "/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;122;add_test;/Users/sesame/Desktop/CECILIA/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")