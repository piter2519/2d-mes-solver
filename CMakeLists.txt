cmake_minimum_required(VERSION 3.10) # zalecane 4.2.0

project(MESProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(MESProject
    src/main.cpp
    src/calculate.cpp
    src/elemUniv.cpp
    src/gauss.cpp
    src/input.cpp
)

target_include_directories(MESProject PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

add_subdirectory(lib/NumLib)
target_link_libraries(MESProject PRIVATE NumLib) 

set_target_properties(MESProject PROPERTIES     
RUNTIME_OUTPUT_DIRECTORY               "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG         "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE       "${CMAKE_SOURCE_DIR}")

if (MSVC)
  target_compile_options(MESProject PRIVATE /utf-8)
endif()

if(MSVC)
    set_property(TARGET MESProject PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()
