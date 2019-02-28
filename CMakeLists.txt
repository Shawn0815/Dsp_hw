# Cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

# Include
INCLUDE_DIRECTORIES(dsp_hw1/c_cpp)

# Project Name
PROJECT(dsp_hw1)

#
#ADD_CUSTOM_TARGET(hw1 COMMAND make -C ${hw1_SOURCE_DIR}
#        CLION_EXE_DIR=${PROJECT_BINARY_DIR} && cp ${hw1_SOURCE_DIR}/test_hmm ./})

# ADD
ADD_EXECUTABLE(test_hmm dsp_hw1/c_cpp/test_hmm.c dsp_hw1/c_cpp/train.c dsp_hw1/c_cpp/test.c)