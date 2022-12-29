#include <def/DelaunayEdgeFlip.h>

#include <CLI/CLI.hpp>

#include <wmtk/utils/Reader.hpp>

#include <igl/Timer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeDMAT.h>


#include <stdlib.h>
#include <chrono>
#include <cstdlib>
#include <iostream>

using namespace wmtk;
using namespace app::def;
using namespace std::chrono;

extern "C" {
#include <wmtk/utils/getRSS.c>
}

int main(int argc, char** argv)
{
    std::cout << "Hello from def\n";
}