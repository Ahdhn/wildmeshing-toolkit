#include <def/DelaunayEdgeFlip.h>

#include <CLI/CLI.hpp>

#include <wmtk/utils/Reader.hpp>

#include <igl/Timer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeDMAT.h>


#include <chrono>
#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace wmtk;
using namespace app::def;
using namespace std::chrono;

extern "C" {
#include <wmtk/utils/getRSS.c>
}

int main(int argc, char** argv)
{
    std::cout << "Hello from def\n";

    std::string input_mesh_name;
    std::string output_mesh_name = "out.obj";
    int thread = 1;

    CLI::App app{argv[0]};

    auto input_option = app.add_option("-i, --input", input_mesh_name, "Input mesh.")->check(CLI::ExistingFile);
    input_option->required();

    auto output_option = app.add_option("-o, --output", output_mesh_name, "Output mesh.");
    app.add_option("-j, --thread", thread, "thread.");

    CLI11_PARSE(app, argc, argv);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    bool ok = igl::read_triangle_mesh(input_mesh_name, V, F);
    assert(ok != false);

    ok = igl::write_triangle_mesh(output_mesh_name, V, F);
    assert(ok != false);

}