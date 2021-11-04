//
// Encodes the topology and topology editing operations of a triangle mesh
//

#pragma once
#include <vector>

namespace wmtk {

    class TriMesh {
    public:
        // Cell Tuple Navigator
        class Tuple {
        public:
            size_t vid;
            size_t eid;
            size_t fid;

            void flip_vertex(const TetMesh &m);//along edge

            void flip_edge(const TetMesh &m);//along face

            void flip_face(const TetMesh &m);//along tet

            size_t get_vertex_attribute_id(const TetMesh &m);
            size_t get_edge_attribute_id(const TetMesh &m);
            size_t get_face_attribute_id(const TetMesh &m);
        };

        class VertexConnectivity {
        public:
            std::vector<size_t> m_conn_tris;
            bool m_is_removed = false;

            inline size_t &operator[](const size_t index) {
                assert(index >= 0 && index < m_conn_tris.size());
                return m_conn_tris[index];
            }

            inline size_t operator[](const size_t index) const {
                assert(index >= 0 && index < m_conn_tris.size());
                return m_conn_tris[index];
            }
        };

        class TriangleConnectivity {
        public:
            std::array<size_t, 3> m_indices;
            bool m_is_removed = false;

            inline size_t &operator[](size_t index) {
                assert(index >= 0 && index < 3);
                return m_indices[index];
            }

            inline size_t operator[](size_t index) const {
                assert(index >= 0 && index < 3);
                return m_indices[index];
            }

            inline int find(int v_id) const {
                for(int j=0;j<3;j++) {
                    if (v_id == m_indices[j])
                        return j;
                }
                return -1;
            }
        };

        virtual ~TriMesh(){};

        void split_edge(const Tuple &t);
        void collapse_edge(const Tuple &t);
        void swapping_edge(const Tuple &t);

    protected:
        // Stores the connectivity of the mesh
        std::vector<VertexConnectivity> m_vertex_connectivity;
        std::vector<TriangleConnectivity> m_tet_connectivity;

        int t_empty_slot = 0;
        int v_empty_slot = 0;
        int find_next_empty_slot_t();
        int find_next_empty_slot_v();

        //// Split the edge in the tuple
        // Checks if the split should be performed or not (user controlled)
        virtual bool split_before(const Tuple &t);//check edge condition
        // This function computes the attributes for the added simplices
        // if it returns false then the operation is undone
        virtual bool split_after(const Tuple &t);//check tet condition

        //// Collapse the edge in the tuple
        // Checks if the collapse should be performed or not (user controlled)
        virtual bool collapse_before(const Tuple &t);
        // If it returns false then the operation is undone (the tuple indexes a vertex and tet that survived)
        virtual bool collapse_after(const Tuple &t);
        //todo: quality, inversion, envelope: change v1 pos before this, only need to change partial attributes

        //// Swap the edge in the tuple
        // Checks if the swapping should be performed or not (user controlled)
        virtual bool swapping_before(const Tuple &t);
        // If it returns false then the operation is undone (the tuple indexes TODO)
        virtual bool swapping_after(const Tuple &t);
        //quality, inversion

        // Invariants that are called on all the new or modified elements after an operation is performed
        virtual bool VertexInvariant(const Tuple &t);
        virtual bool EdgeInvariant(const Tuple &t);
        virtual bool FaceInvariant(const Tuple &t);

        virtual void resize_attributes(size_t v, size_t e, size_t t);

        void compact(); // cleans up the deleted vertices or tetrahedra, and fixes the corresponding indices
    };

}
