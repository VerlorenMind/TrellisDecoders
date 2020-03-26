cmake_minimum_required(VERSION 3.9)
project(beast)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

set(SOURCE_FILES src/BeastDecoder.cpp include/BeastDecoder.h include/myutil.h include/TrellisDecoder.h src/TrellisDecoder.cpp src/myutil.cpp src/BSDDecoder.cpp include/BSDDecoder.h)

if(TESTING)
    set(TEST_FILES tests/maintest.cpp include/catch.hpp tests/test_beast.cpp ${SOURCE_FILES} tests/test_bsd.cpp tests/test_trellis_constructor.cpp include/testutil.h src/testutil.cpp)
    add_executable(test_beast ${TEST_FILES})
else()
    add_executable(beast src/main.cpp ${SOURCE_FILES} )
endif()