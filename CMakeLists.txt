cmake_minimum_required(VERSION 3.10)

# Set the project name
project(NA)

# Set the source directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Find all source files in the directory
file(GLOB SOURCES "${SOURCE_DIR}/*.c")

# Create an executable for each source file
foreach(SOURCE_FILE ${SOURCES})
    get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)
    add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE})
    target_link_libraries(${EXECUTABLE_NAME} m) # Link the math library
endforeach()