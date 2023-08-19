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
    int num_threads = 1;

    CLI::App app{argv[0]};

    auto input_option = app.add_option("-i, --input", input_mesh_name, "Input mesh.")->check(CLI::ExistingFile);
    input_option->required();

    auto output_option = app.add_option("-o, --output", output_mesh_name, "Output mesh.");
    app.add_option("-j, --thread", num_threads, "thread.");

    CLI11_PARSE(app, argc, argv);

    wmtk::logger().set_level(spdlog::level::off);
    wmtk::logger().info("def on {}", input_mesh_name);
    wmtk::logger().info("def output to {}", output_mesh_name);

    // Read in the mesh data
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    bool ok = igl::read_triangle_mesh(input_mesh_name, V, F);
    assert(ok != false);

    // Convert the mesh data into a cleaner, more usable format
    //
    wmtk::logger().info("Before_vertices#: {} \n Before_tris#: {}", V.rows(), F.rows());

    Eigen::VectorXi SVI, SVJ;
    Eigen::MatrixXd temp_V = V; // for STL file
    igl::remove_duplicate_vertices(temp_V, 0, V, SVI, SVJ);
    for (int i = 0; i < F.rows(); i++)
        for (int j : {0, 1, 2}) F(i, j) = SVJ[F(i, j)];

    wmtk::logger().info("After_vertices#: {} \n After_tris#: {}", V.rows(), F.rows());

    std::vector<Eigen::Vector3d> v(V.rows());
    std::vector<std::array<size_t, 3>> tri(F.rows());
    for (int i = 0; i < V.rows(); i++) {
        v[i] = V.row(i);
    }
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) tri[i][j] = (size_t)F(i, j);
    }

    // Delaunay Edge Flip
    //
    DelaunayEdgeFlip def_mesh(v, tri, num_threads);

    igl::Timer timer;
//    for(int i = 0; i < 8; i++) {
//        wmtk::logger().info("*** Delaunay Flip Iteration: {}", i);
//
//        timer.start();
//        def_mesh.swap_all_edges();
//        timer.stop();
//
//        wmtk::logger().info(
//            "***** Delaunay Flip Time *****: {} ms\n any_edges_flipped = {}\n",
//            timer.getElapsedTimeInMilliSec(), def_mesh.any_edges_flipped);
//    }

    int iteration = 0;
    timer.start();
    do {
        wmtk::logger().info("*** Delaunay Flip Iteration: {}", iteration++);


        def_mesh.swap_all_edges();


    } while(def_mesh.any_edges_flipped == true);
    timer.stop();
    wmtk::logger().info(
        "***** Delaunay Flip Time *****: {} ms\n any_edges_flipped = {}\n",
        timer.getElapsedTimeInMilliSec(), def_mesh.any_edges_flipped);  
    std::cout << "***** Delaunay Flip Time *****: " << timer.getElapsedTimeInMilliSec()
              << "ms\n any_edges_flipped = " << def_mesh.any_edges_flipped;
              

    def_mesh.consolidate_mesh();
    def_mesh.write_triangle_mesh(output_mesh_name);
}