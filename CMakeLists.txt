cmake_minimum_required(VERSION 3.5)

project(cw LANGUAGES C)

set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11 -Wall -fsanitize=address,undefined")

include_directories(include)

file(GLOB SOURCES CMAKE_CONFIGURE_DEPENDS
        "${CMAKE_SOURCE_DIR}/src/equations/*.c"
        "${CMAKE_SOURCE_DIR}/src/*.c")

add_executable(cw ${SOURCES})

find_package(GSL REQUIRED)
target_link_libraries(cw PUBLIC gsl gslcblas)
target_link_libraries(cw PRIVATE m)