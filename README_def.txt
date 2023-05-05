# def_app --help
#Usage: ./cmake-build-release/app/def_app [OPTIONS]
#
#Options:
#  -h,--help                   Print this help message and exit
#  -i,--input FILE REQUIRED    Input mesh.
#  -o,--output TEXT            Output mesh.
#  -j,--thread INT             thread.

cmake-build-release/app/def_app --input ./data_ahmed/Nefertiti.obj -j8
cmake-build-release/app/def_app --input ./data/armadillo.obj -j8
