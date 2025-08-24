cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -march=native")
link_libraries(/usr/local/lib/libntl.a)
link_directories(/usr/local/lib)
link_directories(/usr/local/bin)

include_directories(/usr/local/include)


project(DioAppr)

# aux_source_directory(. SRC_LIST)

# add_executable(DioAppr ${SRC_LIST})


# add_executable(test CDT.cpp)

# add_executable(test testsampler.cpp)

# add_executable(test test.cpp)

add_executable(test testFrodo.cpp)
# add_executable(test testFalcon.cpp)

# add_executable(test testSVP.cpp)
# add_executable(test testLinfSVP.cpp)
# add_executable(test test0515.cpp)

# add_executable(test cfrac.cpp)


# add_executable(test SDSVP.cpp)


target_link_libraries(test -lm -lgmp -lntl -lpthread -lssl -lcrypto )