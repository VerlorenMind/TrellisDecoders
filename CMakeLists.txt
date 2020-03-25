cmake_minimum_required(VERSION 3.9)
project(beast)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

set(SOURCE_FILES src/BeastDecoder.cpp include/BeastDecoder.h include/myutil.h src/TrellisDecoder.cpp src/myutil.cpp)

if(TESTING)
    set(TEST_FILES tests/maintest.cpp include/catch.hpp tests/test_constructor.cpp ${SOURCE_FILES})
    add_executable(test_beast ${TEST_FILES})
else()
    add_executable(beast src/main.cpp ${SOURCE_FILES})
endif()