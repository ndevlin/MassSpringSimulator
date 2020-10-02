cmake_minimum_required(VERSION 2.8.2)

project(partio-download NONE)

include(ExternalProject)
ExternalProject_Add(partio-download
                    GIT_REPOSITORY;https://github.com/wdas/partio.git;GIT_TAG;v1.1.0;UPDATE_DISCONNECTED;1;PATCH_COMMAND;/usr/local/bin/git;apply;/Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Deps/partio.patch
                    SOURCE_DIR          "/Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Release/partio-src"
                    BINARY_DIR          "/Users/Nathan/Desktop/PennHwork/563PhysBased/Project1/Project1Repo/563Project1/cispba-master/Build/Release/partio-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
