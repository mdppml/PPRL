cmake_minimum_required(VERSION 3.21)
project(CECILIA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE  "-Wall -Wextra -O3 -flto")

# Libraries:
find_package(Threads REQUIRED)

# record linkage
add_executable(record_linkage apps/record_linkage/main_record_linkage.cpp)
target_link_libraries(record_linkage PRIVATE Threads::Threads cryptopp::cryptopp)


# include dependencies using git submodule:
add_subdirectory(dependencies/eigen)
add_subdirectory(dependencies/cryptopp-cmake)
