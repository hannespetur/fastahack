cmake_minimum_required(VERSION 2.8)

SET_SOURCE_FILES_PROPERTIES(disorder.c PROPERTIES LANGUAGE CXX)

# Create a fastahack library
set(fastahack_SOURCE_FILES
    Fasta.cpp
    FastaHackAPI.cpp
    disorder.c
    split.cpp
)

add_library(fastahack_core ${fastahack_SOURCE_FILES} )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(fastahack FastaHack.cpp)
# find_package( Threads REQUIRED )
target_link_libraries(fastahack fastahack_core pthread ${ZLIB_LIBRARIES})


# Add the executable
target_link_libraries(fastahack ${LIBS})

# Installs fastahack to /usr/local/bin using [sudo] make install
install (TARGETS fastahack DESTINATION bin)
