cmake_minimum_required(VERSION 3.22)

project(DNAStorage)
message("project name is automaically set as: ${PROJECT_NAME}")

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

enable_language(CXX)
enable_language(CUDA)
set(PROJECT_VERSION 3.10.4)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /Qvec-report:1")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# add_subdirectory(src/common)
# add_subdirectory(src/ldpc)
add_subdirectory(src/idec)
# add_subdirectory(src/dnas)
add_subdirectory(src/ldpcidec)

target_link_libraries(ldpcidec PRIVATE ldpc idec)

add_executable(aux ${PROJECT_SOURCE_DIR}/src/auxiliary.cpp)

add_executable(cawgn ${PROJECT_SOURCE_DIR}/src/mainawgn.cu)
add_executable(awgn ${PROJECT_SOURCE_DIR}/src/mainawgn.cpp)
add_executable(mids ${PROJECT_SOURCE_DIR}/src/mainmids.cpp)
add_executable(lids ${PROJECT_SOURCE_DIR}/src/mainlids.cpp)

add_executable(dvnb ${PROJECT_SOURCE_DIR}/src/maindvnb.cpp)
add_executable(dvdb ${PROJECT_SOURCE_DIR}/src/maindvdb.cpp)
add_executable(ddnb ${PROJECT_SOURCE_DIR}/src/mainddnb.cpp)
add_executable(dddb ${PROJECT_SOURCE_DIR}/src/maindddb.cpp)

target_link_libraries(mids PRIVATE idec)
target_link_libraries(lids PRIVATE ldpcidec)
