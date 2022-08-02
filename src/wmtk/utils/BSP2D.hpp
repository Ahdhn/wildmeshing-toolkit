#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <wmtk/utils/Rational.hpp>

namespace wmtk {
namespace bsp2d {

typedef std::array<Rational, 2> Point;

// class Point
// {
// public:
//     Rational x;
//     Rational y;
// };

/**
 * Compute the exact intersection of a segment ab and a segment cd
 * @param[a] first point of segment ab
 * @param[b] second point of segment ab
 * @param[c] first point of segment cd
 * @param[d] second point of segment cd
 * @param[is_cross_c] true if ab intersects c
 * @param[is_cross_d] true if cd intersects d
 * @param[intersection] the intersection
 *
 * @returns true if the two segments intersect
 */
bool segment_intersection(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d,
    bool& is_cross_c,
    bool& is_cross_d,
    Point& intersection);


// void triwild::triangulation::BSP_subdivision(const Eigen::MatrixXd& V, const
// std::vector<std::array<int, 2>>& edges,
//         MeshData& mesh, std::vector<std::vector<int>>& tag_boundary_es, GEO::MeshFacetsAABB&
//         b_tree) {
//     ////delaunay
//     Eigen::MatrixXd new_V = V;
//     const int n = new_V.rows();
// //    new_V.conservativeResize(n + 4, 2);
//     Eigen::RowVector2d min(args.box_min(0) - args.target_edge_len, args.box_min(1) -
//     args.target_edge_len); Eigen::RowVector2d max(args.box_max(0) + args.target_edge_len,
//     args.box_max(1) + args.target_edge_len);
// //    new_V.row(n) = min;
// //    new_V.row(n + 1) << min(0), max(1);
// //    new_V.row(n + 2) = max;
// //    new_V.row(n + 3) << max(0), min(1);

//     //add voxel points
//     int nx = (max[0] - min[0])/(args.diagonal_len/10) + 1.5;
//     int ny = (max[1] - min[1])/(args.diagonal_len/10) + 1.5;
//     double threshold_2 = (args.diagonal_len/1000)*(args.diagonal_len/1000);
//     for(double i=0;i<=nx;i++) {
//         for (double j = 0; j <= ny; j++) {
//             Eigen::RowVector2d p(i / nx * min[0] + (1 - i / nx) * max[0], j / ny * min[1] + (1 - j / ny) * max[1]);
//             GEO::vec3 p0(p[0], p[1], 0);
//             GEO::vec3 nearest_point;
//             double sq_dist;
//             GEO::index_t prev_facet = b_tree.nearest_facet(p0, nearest_point, sq_dist);
//             if (sq_dist < threshold_2)
//                 continue;

//             new_V.conservativeResize(new_V.rows() + 1, 2);
//             new_V.row(new_V.rows() - 1) = p;
//         }
//     }

//     GEO::Delaunay::initialize();
//     GEO::Delaunay_var delaunay2d = GEO::Delaunay::create(2, "BDEL2d");
// //    delaunay2d->set_vertices(new_V.rows(), new_V.data());
//     auto *flatten = new double[new_V.rows() * 2];
//     for (int i = 0; i < new_V.rows(); i++) {
//         flatten[i * 2] = new_V(i, 0);
//         flatten[i * 2 + 1] = new_V(i, 1);
//     }
//     delaunay2d->set_vertices(new_V.rows(), flatten);
//     delete[] flatten;
//     auto cell_adj = delaunay2d->cell_to_v();

//     // Copy triangles to initialize the bsp tree
//     std::vector<std::vector<int>> bsp_faces;
//     bsp_faces.resize(delaunay2d->nb_cells());
//     for (int i = 0; i < delaunay2d->nb_cells(); i++)
//         bsp_faces[i] = {{cell_adj[i * 3], cell_adj[i * 3 + 1], cell_adj[i * 3 + 2]}};

//     // Counts how many input edges are incident on a given vertex
//     std::vector<int> cnt_conn_es(V.rows(), 0);
//     for (int i = 0; i < edges.size(); i++) {
//         cnt_conn_es[edges[i][0]]++;
//         cnt_conn_es[edges[i][1]]++;
//     }
//     // The vertices of the BSP are the same as the input mesh
//     // 4 corners are fixed
//     // initially all vertices have the on point property
//     auto &bsp_vertices = mesh.tri_vertices;
//     bsp_vertices.resize(new_V.rows());
//     for (int i = 0; i < new_V.rows(); i++) {
//         bsp_vertices[i].pos[0] = new_V(i, 0);
//         bsp_vertices[i].pos[1] = new_V(i, 1);
//         if ((bsp_vertices[i].pos[0] == min[0] || bsp_vertices[i].pos[0] == max[0])
//             && (bsp_vertices[i].pos[1] == min[1] || bsp_vertices[i].pos[1] == max[1]))
//             bsp_vertices[i].is_freezed = true;//freeze bbox corners
//         if (i < n && cnt_conn_es[i] == 1) {
//             bsp_vertices[i].is_on_point = true;
//             bsp_vertices[i].input_posf.set(new_V(i, 0), new_V(i, 1));
//         }
//     }

//     // conn_fs stores all faces incident to each vertex
//     std::vector<std::unordered_set<int>> conn_fs(bsp_vertices.size());
//     for (int i = 0; i < bsp_faces.size(); i++) {
//         for (int j = 0; j < 3; j++)
//             conn_fs[bsp_faces[i][j]].insert(i);
//     }

//     // for every vertex, store the id of the edges to insert that are touching it
//     tag_boundary_es.resize(bsp_vertices.size());

//     int old_cnt_v = bsp_vertices.size();
//     int old_cnt_f = bsp_faces.size();

//     //// Cuts a face f_id into two
//     // r_v_id is the vertex before the cut i think
//     // cnt is how many vertices to keep in the first face
//     auto split_a_face = [&](int r_v_id, int cnt, int n_size, int f_id) -> int {
//         std::vector<int> new_face;

//         // Creates the first half face
//         for (int k = 0; k <= cnt; k++)//p1->
//             new_face.push_back(bsp_faces[f_id][(r_v_id + k) % n_size]);
//         bsp_faces.push_back(new_face);
//         int new_f_id = bsp_faces.size() - 1;

//         // Creates the second half face
//         new_face.clear();
//         for (int k = cnt; k <= n_size; k++)//->p2
//             new_face.push_back(bsp_faces[f_id][(r_v_id + k) % n_size]);
//         bsp_faces[f_id] = new_face;

//         // Updates the vertex connectivity
//         for (int v_id:bsp_faces[new_f_id]) {
//             conn_fs[v_id].erase(f_id);
//             conn_fs[v_id].insert(new_f_id);
//         }
//         for (int v_id:bsp_faces[f_id]) {
//             conn_fs[v_id].insert(f_id);
//         }

//         return new_f_id;
//     };

//     TriVertex intersection_v;
//     for (int e_id = 0; e_id < edges.size(); e_id++) {

//         int v1_id = edges[e_id][0];
//         int v2_id = edges[e_id][1];
//         tag_boundary_es[v1_id].push_back(e_id);

//         // First make sure that the edge does not exist, if it does skip
//         std::vector<int> tmp = optimization::set_intersection(conn_fs[v1_id],
//         conn_fs[v2_id]);
//         if (tmp.size() > 0) { //edge exist in the init mesh
//             tag_boundary_es[v2_id].push_back(e_id);
//             continue;
//         }

//         // Check to make sure the edge is not degenerate (and if it is? skip?)
//         assert(v1_id != v2_id);
//         assert(bsp_vertices[v1_id].pos != bsp_vertices[v2_id].pos);

//         int v_id = v1_id;
//         std::unordered_set<int> tmp_conn_fs = conn_fs[v1_id];
//         while (true) {
//             bool is_finished = false;
//             bool is_intersected = false;
//             // For all faces touching the first vertex of the edge (or later on the last inserted)
//             for (int f_id: tmp_conn_fs) {
//                 is_intersected = false;
//                 int n_size = bsp_faces[f_id].size();
//                 // find the starting vertex inside the face
//                 int r_v_id = std::find(bsp_faces[f_id].begin(), bsp_faces[f_id].end(), v_id)
//                 - bsp_faces[f_id].begin();
//                 // for every edge i,j of the face (skipping first and last)
//                 for (int i = 1; i < n_size - 1; i++) {
//                     const int j = (r_v_id + i) % n_size;
//                     // get the edge
//                     int p1_id = bsp_faces[f_id][j];
//                     int p2_id = bsp_faces[f_id][(j + 1) % n_size];
//                     // This should never happen! Assert this instead
//                     if (p1_id == v_id || p2_id == v_id)
//                         continue;

//                     // check for intersection between the edge e_id and the side of the polygon
//                     bool is_cross_p1 = false;
//                     bool is_cross_p2 = false;
//                     if (!segment_intersection(mesh, v_id, v2_id, p1_id, p2_id, is_cross_p1,
//                     is_cross_p2,
//                                               intersection_v))
//                         continue;

//                     // if the edge passes exactly through p1
//                     if (is_cross_p1) {
//                         // we reached the end of the tracing, stop after this
//                         if (p1_id == v2_id) {
//                             is_finished = true;
//                         }
//                         tag_boundary_es[p1_id].push_back(e_id);

//                         tmp_conn_fs = conn_fs[p1_id];
//                         std::vector<int> n_f_ids =
//                         optimization::set_intersection(conn_fs[v_id], conn_fs[p1_id]); for
//                         (int n_f_id:n_f_ids)
//                             tmp_conn_fs.erase(n_f_id);

//                         // if it is not the first segment? not sure how or why this can happen, assert?
//                         if (i != 1) {
//                             int new_f_id = split_a_face(r_v_id, i, n_size, f_id);
//                         }
//                         v_id = p1_id;
//                     // the edge is touching exactly p2
//                     } else if (is_cross_p2) {
//                         // if we reached the end stop
//                         if (p2_id == v2_id) {
//                             is_finished = true;
//                         }
//                         tag_boundary_es[p2_id].push_back(e_id);

//                         tmp_conn_fs = conn_fs[p2_id];
//                         std::vector<int> n_f_ids =
//                         optimization::set_intersection(conn_fs[v_id], conn_fs[p2_id]); for
//                         (int n_f_id:n_f_ids)
//                             tmp_conn_fs.erase(n_f_id);
//                         // not sure how this could be, assert?
//                         if (i != n_size - 2) {
//                             int new_f_id = split_a_face(r_v_id, i + 1, n_size, f_id);
//                         }
//                         v_id = p2_id;
//                     // generic case, the two edges do not intersect at a vertex
//                     } else {
//                         // Add a new vertex to the BSP (and update data-structures)
//                         bsp_vertices.push_back(intersection_v);
//                         int new_v_id = bsp_vertices.size() - 1;
//                         conn_fs.resize(conn_fs.size() + 1);
//                         tag_boundary_es.resize(tag_boundary_es.size() + 1);

//                         // Update the tagging of boundary with all edges touched by both vertices of the edge
//                         std::vector<int> p1p2_e_ids =
//                         optimization::set_intersection(tag_boundary_es[p1_id],
//                         tag_boundary_es[p2_id]); for(int p1p2_e_id:p1p2_e_ids)
//                             tag_boundary_es[new_v_id].push_back(p1p2_e_id);

//                         // Update the existing face with the additional point to get ready for a split
//                         bsp_faces[f_id].insert(bsp_faces[f_id].begin() + j + 1, new_v_id);
//                         std::vector<int> n_f_ids =
//                         optimization::set_intersection(conn_fs[p1_id], conn_fs[p2_id]); int
//                         n_f_id = n_f_ids[0] == f_id ? n_f_ids[1] : n_f_ids[0]; for (int k =
//                         0; k < bsp_faces[n_f_id].size(); k++) {
//                             if (bsp_faces[n_f_id][k] == p2_id
//                                 && bsp_faces[n_f_id][(k + 1) % bsp_faces[n_f_id].size()] ==
//                                 p1_id) { bsp_faces[n_f_id].insert(bsp_faces[n_f_id].begin() +
//                                 k + 1, new_v_id); break;
//                             }
//                         }

//                         // I don't understand this
//                         if (j + 1 <= r_v_id)
//                             r_v_id++;

//                         // Split the face
//                         int new_f_id = split_a_face(r_v_id, i + 1, n_size + 1, f_id);

//                         // Update connectivity and tags
//                         conn_fs[new_v_id].insert(n_f_id);
//                         tag_boundary_es[new_v_id].push_back(e_id);

//                         // We added a new vertex, recursively continue splitting
//                         v_id = new_v_id;
//                         tmp_conn_fs.clear();
//                         tmp_conn_fs.insert(n_f_id);
//                     }

//                     is_intersected = true;
//                     break;
//                 }

//                 if (is_intersected)
//                     break;
//             }
//             if (is_finished)
//                 break;
//         }
//     }

//     // Why converting rational to double?
//     for (auto &v: mesh.tri_vertices) {
//         v.posf[0] = v.pos[0].to_double();
//         v.posf[1] = v.pos[1].to_double();
//     }


//     mesh.tris.reserve(bsp_faces.size());
//     for (int f_id = 0; f_id < bsp_faces.size(); f_id++) {
//         assert(bsp_faces[f_id].size() >= 3);
//         // Single triangle can be copied
//         if (bsp_faces[f_id].size() == 3) {
//             mesh.tris.push_back({{bsp_faces[f_id][0], bsp_faces[f_id][1],
//             bsp_faces[f_id][2]}}); continue;
//         }
//         // Polygon is triangolated (probably adding a vertex is a better way)
//         for (int j = 1; j < int(bsp_faces[f_id].size()) - 1; j++) {
//             mesh.tris.push_back({{bsp_faces[f_id][0], bsp_faces[f_id][j], bsp_faces[f_id][j +
//             1]}});
//         }
//     }

//     //output and check

//     cout << "#v " << old_cnt_v << "->" << mesh.tri_vertices.size() << endl;
//     cout << "#f " << old_cnt_f << "->" << mesh.tris.size() << endl;
// }

} // namespace bsp2d
} // namespace wmtk
