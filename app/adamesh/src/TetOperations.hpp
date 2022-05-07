#pragma once

#include <wmtk/TetMesh.h>
#include <wmtk/TetMeshOperations.hpp>
#include "wmtk/utils/VectorUtils.h"

constexpr auto replace = [](auto& arr, size_t a, size_t b) {
    for (auto i = 0; i < arr.size(); i++) {
        if (arr[i] == a) {
            arr[i] = b;
            return i;
        }
    }
    return -1;
};

namespace wmtk {
struct SplitFace : public TetMesh::OperationBuilder
{
    const TetMesh& m;
    std::vector<size_t> affected;
    std::array<size_t, 3> tri_verts;
    size_t ux;

    SplitFace(const TetMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        auto oppo = m.switch_tetrahedron(t).value();
        affected = {t.tid(m), oppo.tid(m)};

        return affected;
    }
    int request_vert_slots() { return 1; }
    std::vector<std::array<size_t, 4>> replacing_tets(const std::vector<size_t>& slots)
    {
        assert(slots.size() == 1);
        ux = slots.front();

        auto new_tets = std::vector<std::array<size_t, 4>>();

        new_tets.reserve(2 * 3);
        for (auto i = 0; i < 2; i++) {
            auto t_conn = m.oriented_tet_vids(m.tuple_from_tet(affected[i]));
            for (auto j = 0; j < 3; j++) {
                new_tets.push_back(t_conn);
                replace(new_tets.back(), tri_verts[j], ux);
            }
        }
        return new_tets;
    }
};


struct DivideTet : public TetMesh::OperationBuilder
{
    const TetMesh& m;
    TetMesh::Tuple tet;
    size_t ux;

    DivideTet(const TetMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        tet = t;
        return {t.tid(m)};
    }
    int request_vert_slots() { return 1; }
    std::vector<std::array<size_t, 4>> replacing_tets(const std::vector<size_t>& slots)
    {
        assert(slots.size() == 1);
        ux = slots.front();

        std::array<size_t, 4> t_conn;
        auto vs = m.oriented_tet_vertices(tet);
        for (auto i = 0; i < 4; i++) t_conn[i] = vs[i].vid(m);
        auto new_tets = std::vector<std::array<size_t, 4>>(4, t_conn);
        for (auto i = 0; i < 4; i++) new_tets[i][i] = ux;

        return new_tets;
    }
};

struct CollapseEdge : public TetMesh::OperationBuilder
{
    const TetMesh& m;
	
    std::vector<size_t> affected; // nb of v1
	std::array<size_t, 2> edge_verts;

public:
    CollapseEdge(const TetMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        auto t1 = t.switch_vertex(m);
		affected = m.get_one_ring_tids_for_vertex(t);
		edge_verts = {t.vid(m), t1.vid(m)};
        return affected;
    }
    std::vector<std::array<size_t, 4>> replacing_tets(const std::vector<size_t>&)
    {
        auto new_tet_conn = std::vector<std::array<size_t, 4>>();
        std::vector<size_t> preserved_tids;
		auto [v1_id, v2_id] = edge_verts;
        for (auto t_id : affected) {
            auto t_conn = m.oriented_tet_vids(m.tuple_from_tet(t_id));
			if (replace(t_conn, v2_id, v2_id) != -1) { // already has v2
				continue;
			}
			replace(t_conn, v1_id, v2_id);
            
            new_tet_conn.push_back(t_conn);
            preserved_tids.push_back(t_id);
        }
		return new_tet_conn;
    }
};
} // namespace wmtk

