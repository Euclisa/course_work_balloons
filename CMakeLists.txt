cmake_minimum_required(VERSION 3.5)

project(cw LANGUAGES C)

set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11 -Wall") # -fsanitize=address,undefined
add_definitions(-DCW_PROJECT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK3 REQUIRED gtk+-3.0)

include_directories(include)

file(GLOB SOURCES CMAKE_CONFIGURE_DEPENDS
        "${CMAKE_SOURCE_DIR}/src/equations/*.c"
        "${CMAKE_SOURCE_DIR}/src/*.c")

add_executable(cw ${SOURCES})

find_package(GSL REQUIRED)
target_link_libraries(cw PUBLIC gsl gslcblas)
target_link_libraries(cw PRIVATE m)
target_link_libraries(cw PRIVATE ${GTK3_LIBRARIES})
target_include_directories(cw PRIVATE ${GTK3_INCLUDE_DIRS})