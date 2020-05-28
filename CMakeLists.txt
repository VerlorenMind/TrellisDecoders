cmake_minimum_required(VERSION 3.9)
project(beast)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

set(SOURCE_FILES src/BeastDecoder.cpp include/BeastDecoder.h include/Utility.h include/TrellisDecoder.h src/TrellisDecoder.cpp src/Utility.cpp src/BSDDecoder.cpp include/BSDDecoder.h include/Simulation.h src/Simulation.cpp include/SoftDecoder.h include/DecoderID.h include/OrderedStatisticsDecoder.h src/OrderedStatisticsDecoder.cpp src/SoftDecoder.cpp include/Trellis.h src/Trellis.cpp include/ViterbiDecoder.h src/ViterbiDecoder.cpp)

if(TESTING)
    add_compile_definitions(CATCH_TESTING)
    set(TEST_FILES tests/maintest.cpp include/catch.hpp tests/test_beast.cpp ${SOURCE_FILES} tests/test_bsd.cpp tests/test_trellis_constructor.cpp include/CatchWrap.h tests/test_simulation.cpp tests/test_utility.cpp tests/test_osd.cpp)
    add_executable(test_beast ${TEST_FILES})
else()
    add_executable(beast src/main.cpp ${SOURCE_FILES})
endif()