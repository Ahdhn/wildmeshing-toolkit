#pragma once

#include "TetMesh.h"

namespace wmtk {

    class VertexAttributes {
    public:
        Vector3 m_pos;
        Vector3f m_posf;

        bool m_is_on_surface;
        bool m_is_on_boundary;
        bool m_is_on_bbox;
        bool m_is_outside;

        Scalar m_sizing_scalars;
        Scalar m_scalars;
        bool m_is_freezed;
    };

    class EdgeAttributes {
    public:
        // Scalar length;
    };

    class FaceAttributes {
    public:
        Scalar tag;

        int m_is_surface_fs;
        int m_is_bbox_fs;
        int m_opp_t_ids;
        int m_surface_tags;
    };

    class TetrahedronAttributes {
    public:
        Scalar m_qualities;
        Scalar m_scalars;
        bool m_is_outside;
    };

    class TetWild : public TetMesh {
    public:
        Parameters& m_params;
        Envelope& m_envelope;

        TetWild(Parameters& _m_params, Envelope& _m_envelope): m_params(_m_params), m_envelope(_m_envelope){}

        // Stores the attributes attached to simplices
        std::vector<VertexAttributes> m_vertex_attribute;
        std::vector<EdgeAttributes> m_edge_attribute;
        std::vector<FaceAttributes> m_face_attribute;
        std::vector<TetrahedronAttributes> m_tet_attribute;

        void resize_attributes(size_t v, size_t e, size_t t, size_t tt) {
            m_vertex_attribute.resize(v);
            m_edge_attribute.resize(e);
            m_face_attribute.resize(t);
            m_tet_attribute.resize(tt);
        }

        void smoothing(const Tuple &t);

        // all the other functions
    private:
        bool split_before(const Tuple &t) override;
        bool split_after(const Tuple &t) override;

        bool is_inverted(size_t t_id);
        double get_quality(size_t t_id);

    };

}