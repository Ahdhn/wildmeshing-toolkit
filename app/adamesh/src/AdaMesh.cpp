#include "AdaMesh.hpp"

#include <array>
#include <cassert>
#include <vector>
#include <wmtk/TetMeshOperations.hpp>
#include "prism/geogram/AABB.hpp"
#include "prism/geogram/AABB_tet.hpp"
#include "prism/predicates/inside_prism_tetra.hpp"
#include "spdlog/spdlog.h"
#include "tetwild/TetWild.h"
#include "wmtk/ConcurrentTetMesh.h"
#include "wmtk/TetMesh.h"
#include "wmtk/utils/InsertTriangleUtils.hpp"
#include "wmtk/utils/Logger.hpp"
#include "TetOperations.hpp"

#include <geogram/numerics/predicates.h>


bool adamesh::SplitFace::before(const wmtk::TetMesh::Tuple& loc)
{
    cache.face_tag = m_app.m_face_attribute[loc.fid(m)];

    for (auto& t : {loc, loc.switch_tetrahedron(m).value()}) {
        auto vs = m.oriented_tet_vids(t);
        for (int j = 0; j < 4; j++) {
            std::array<size_t, 3> f_vids = {{vs[(j + 1) % 4], vs[(j + 2) % 4], vs[(j + 3) % 4]}};
            std::sort(f_vids.begin(), f_vids.end());
            auto [_, global_fid] = m.tuple_from_face(f_vids);
            cache.f_attrs.emplace(f_vids, m_app.m_face_attribute[global_fid]);
        }
    }
    return true;
}

bool adamesh::SplitFace::after(const std::vector<wmtk::TetMesh::Tuple>& locs)
{
    // vertex
    m_app.m_vertex_attribute[ux] = tetwild::VertexAttributes(cache.pos);
    if (cache.face_tag.m_is_surface_fs) m_app.m_vertex_attribute[ux].m_is_on_surface = true;
    if (cache.face_tag.m_is_bbox_fs >= 0)
        m_app.m_vertex_attribute[ux].on_bbox_faces.push_back(cache.face_tag.m_is_bbox_fs);

    // face: 2 tet -> 6 tet, bounding faces directly inherit,
    for (auto& [old_vids, attr] : cache.f_attrs) {
        std::vector<int> j_vn; // vertex apart from  tri
        for (int j = 0; j < 3; j++) {
            if (old_vids[j] != tri_verts[0] && old_vids[j] != tri_verts[1] &&
                old_vids[j] != tri_verts[2]) {
                j_vn.push_back(j);
            }
        }

        if (j_vn.size() == 0) { // tri 0-1-2 of interest, split to 3.
            for (auto j = 0; j < 3; j++) {
                auto vids = old_vids;
                vids[j] = ux;
                auto [_, global_fid] = m.tuple_from_face(vids);
                m_app.m_face_attribute[global_fid] = attr;
            }
        } else if (j_vn.size() == 1) { // tri 0-1-T/B
            auto [_, global_fid] = m.tuple_from_face(old_vids); // inherit boundary
            m_app.m_face_attribute[global_fid] = attr;

            for (auto j = 0; j < 3; j++) { // internal (0, x, T), duplicate appear here.
                auto [_, global_fid] =
                    m.tuple_from_face({{ux, tri_verts[j], old_vids[j_vn.front()]}});
                m_app.m_face_attribute[global_fid].reset();
            }
        } else { // j_vn.size() == 2
            assert(false);
        }
    }


    // tet
    for (auto& loc : locs) {
        m_app.m_tet_attribute[loc.tid(m)].m_quality = m_app.get_quality(loc);
    }

    return true;
}


bool adamesh::DivideTet::before(const wmtk::TetMesh::Tuple& t)
{
    auto vs = m.oriented_tet_vids(t);
    for (int j = 0; j < 4; j++) {
        std::array<size_t, 3> f_vids = {{vs[(j + 1) % 4], vs[(j + 2) % 4], vs[(j + 3) % 4]}};
        std::sort(f_vids.begin(), f_vids.end());
        auto [_, global_fid] = m.tuple_from_face(f_vids);
        cache.old_fids.push_back(global_fid);
        cache.f_attrs.emplace_back(f_vids, m_app.m_face_attribute[global_fid]);
        cache.old_fids.push_back(global_fid);
    }
    return true;
}

bool adamesh::DivideTet::after(const std::vector<wmtk::TetMesh::Tuple>& locs)
{
    // vertex
    m_app.m_vertex_attribute[ux] = tetwild::VertexAttributes(cache.pos);

    // face: 2 tet -> 6 tet, bounding faces directly inherit,
    for (auto f : cache.old_fids) {
        m_app.m_face_attribute[f].reset();
    }
    for (auto& [old_vids, attr] : cache.f_attrs) {
        auto [_, global_fid] = m.tuple_from_face(old_vids);
        m_app.m_face_attribute[global_fid] = attr;
    }

    // tet
    for (auto& loc : locs) {
        m_app.m_tet_attribute[loc.tid(m)].m_quality = m_app.get_quality(loc);
    }

    return true;
}
