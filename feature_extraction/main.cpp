#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <iostream>
#include <string>
#include <CGAL/draw_surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef Mesh::Vertex_index vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

void neighbours_add(halfedge_descriptor &h, Mesh &mesh, Mesh &newMesh, int k,
                    std::vector<face_descriptor> &facesAdded) {
    if (k <= 0) return;
    CGAL::Halfedge_around_face_iterator<Mesh> fbegin, fend;
    for (boost::tie(fbegin, fend) = halfedges_around_face(h, mesh); fbegin != fend; ++fbegin) {
        std::vector<vertex_descriptor> point_list;
        std::vector<K::Point_3> point_list_raw;
        CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
        if (!is_border(opposite(*fbegin, mesh), mesh)) {
            for (boost::tie(vbegin, vend) = vertices_around_face(opposite(*fbegin, mesh), mesh); vbegin != vend; ++
                 vbegin) {
                point_list_raw.push_back(mesh.point(*vbegin));
            }
            halfedge_descriptor nextHe = opposite(*fbegin, mesh);
            neighbours_add(nextHe, mesh, newMesh, --k, facesAdded);
            face_descriptor f = face(opposite(*fbegin, mesh), mesh);
            if (find(facesAdded.begin(), facesAdded.end(), f) == facesAdded.end()) {
                for (auto pt: point_list_raw) {
                    point_list.push_back(newMesh.add_vertex(pt));
                }
                facesAdded.push_back(f);
                newMesh.add_face(point_list);
            }
        }
    }
}


void add_mesh(halfedge_descriptor &h, Mesh &mesh, Mesh &newMesh, int k, std::vector<face_descriptor> &facesAdded) {
    if (!is_border(h, mesh)) {
        neighbours_add(h, mesh, newMesh, k, facesAdded);
        face_descriptor f = face(h, mesh);
        if (find(facesAdded.begin(), facesAdded.end(), f) == facesAdded.end()) {
            facesAdded.push_back(f);

            std::vector<vertex_descriptor> point_list;
            CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
            for (boost::tie(vbegin, vend) = vertices_around_face(h, mesh); vbegin != vend; ++vbegin) {
                point_list.push_back(newMesh.add_vertex(mesh.point(*vbegin)));
            }
            newMesh.add_face(point_list);
        }
    }
}


int main(int argc, char *argv[]) {
    const std::string filename = (argc > 1)
                                     ? argv[1]
                                     : CGAL::data_file_path("../test/neo.obj");
    Mesh mesh;
    if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
        std::cerr << "Invalid input." << std::endl;
        return 1;
    }
    typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;
    EIFMap eif = get(CGAL::edge_is_feature, mesh);
    PMP::detect_sharp_edges(mesh, 60, eif);
    std::size_t nb_sharp_edges = 0;
    std::vector<face_descriptor> facesAdded;
    Mesh newMesh;
    for (boost::graph_traits<Mesh>::edge_descriptor e: edges(mesh)) {
        if (get(eif, e)) {
            ++nb_sharp_edges;
            halfedge_descriptor h = halfedge(e, mesh);
            add_mesh(h, mesh, newMesh, 0, facesAdded);
            halfedge_descriptor oppositeH = opposite(h, mesh);
            add_mesh(oppositeH, mesh, newMesh, 0, facesAdded);
        }
    }

    std::ofstream outfile;
    outfile.open("../results/neo_wire_90.ply");
    CGAL::IO::write_PLY(outfile, newMesh);

    return 0;
}
