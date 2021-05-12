cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(poly VERSION 1.0 DESCRIPTION "Single tethered polymer dynamics"
        LANGUAGES CXX)




find_package(Eigen3 3.3 REQUIRED NO_MODULE)


add_subdirectory(src)
add_subdirectory(main)
