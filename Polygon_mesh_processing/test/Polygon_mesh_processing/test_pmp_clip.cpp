#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <iostream>
#include <fstream>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template <class TriangleMesh>
void test()
{
  // test with a clipper mesh
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/elephant.off") >> tm1;
    std::ifstream("data-coref/sphere.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false).face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(!CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/elephant.off") >> tm1;
    std::ifstream("data-coref/sphere.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(true).face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(CGAL::is_closed(tm1));
  }

  // test with a iso-cuboid
  {
    TriangleMesh tm1;
    std::ifstream("data-coref/elephant.off") >> tm1;
    K::Iso_cuboid_3 iso_cuboid(K::Point_3(0,0,0), K::Point_3(0.4, 0.6, 0.4));

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);

    PMP::clip(tm1, iso_cuboid, params::clip_volume(true).face_index_map(custom_face_index_map_1));
    assert(CGAL::is_closed(tm1));
  }

  // test with a plane
  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;
    K::Plane_3 plane(0, 0, 1, -1);

    PMP::clip(tm1, plane, params::clip_volume(true));
    assert(CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;
    K::Plane_3 plane(0, 0, 1, -1);

    PMP::clip(tm1, plane, params::clip_volume(false).use_compact_clipper(false));
    assert(!CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;
    K::Plane_3 plane(0, 0, 1, -1);

    PMP::clip(tm1, plane, params::clip_volume(false).use_compact_clipper(true));
    assert(CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    PMP::clip(tm1, K::Plane_3(-0.236474, 0.437732, 0.867451, -0.838791), params::clip_volume(true));
    assert(CGAL::is_closed(tm1));
    assert(!CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    PMP::clip(tm1, K::Plane_3(0, 0, 1, 2));
    assert(CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    PMP::clip(tm1, K::Plane_3(0, 0, 1, -2));
    assert(!CGAL::is_empty(tm1));
  }

  // clipping with identity
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(true)
                     .use_compact_clipper(true)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(num_vertices(tm1) == 8);
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(true)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(num_vertices(tm1) == 8);
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(true)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 4);
  }

  // test orientation + patch without input vertex
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;

    CGAL::make_tetrahedron(K::Point_3(0.53, -1.3, 0.2),
                           K::Point_3(0.53, 1.1, 0.2),
                           K::Point_3(0.53, -1.3, 0.4),
                           K::Point_3(0.73, -1.3, 0.2),
                           tm2);

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 6);
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;

    CGAL::make_tetrahedron(K::Point_3(0.53, -1.3, 0.2),
                           K::Point_3(0.53, 1.1, 0.2),
                           K::Point_3(0.53, -1.3, 0.4),
                           K::Point_3(0.73, -1.3, 0.2),
                           tm2);
    PMP::reverse_face_orientations(tm2);

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 6+8);
  }

  // clip meshes with intersection polyline opened
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(1, 0, 0, -2));
    assert(vertices(tm1).size() == 4);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(-1, 0, 0, 2));
    assert(vertices(tm1).size() == 3);
  }

  // test with clipper on border edge
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, 1, 0 , 0));
    assert(vertices(tm1).size() == 0);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, -1, 0 , 0));
    assert(vertices(tm1).size() == 4);
  }

  // test with clipper on border edge: full triangle
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(true));
    assert(vertices(tm1).size()!=0);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(false));
    assert(vertices(tm1).size() == 0);
  }

  // test tangencies
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(1, 0, 0, -1));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(-1, 0, 0, 1));
    assert(vertices(tm1).size() == 0);
  }

  {
    TriangleMesh tm1, tm2;
    make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1, tm2;
    make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::reverse_face_orientations(tm2);
    PMP::clip(tm1, tm2,
              params::face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 0);
  }

  // test combinaison of use_compact_clipper and clip_volume
  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    //  -> closed mesh, true/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(true));
    assert(faces(tm1).size() == 12);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, false/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(false).clip_volume(true));
    assert(faces(tm1).size() == 12);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, true/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(false));
    assert(faces(tm1).size() == 12);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, false/false
    PMP::clip(tm1, K::Plane_3(1,0,0,-1), params::use_compact_clipper(false).clip_volume(false));
    assert(faces(tm1).size() == 10);
    assert(!CGAL::is_closed(tm1));

    // -> open mesh true/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(true));
    assert(faces(tm1).size() == 10);

    // -> open mesh true/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(false));
    assert(faces(tm1).size() == 10);

    // -> open mesh false/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(false).clip_volume(false));
    assert(faces(tm1).size() == 8);

    // -> open mesh false/true
    PMP::clip(tm1, K::Plane_3(0,-1,0,0), params::use_compact_clipper(false).clip_volume(true));
    assert(faces(tm1).size() == 6);
  }

  // test special case
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-clip/tm_1.off") >> tm2;
    std::ifstream("data-clip/clipper_1.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(is_valid_polygon_mesh(tm1));
  }

  // non-manifold border vertices
  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 5 2 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 0 1 4\n3 1 2 3\n";
    ss >> tm1;
    PMP::clip(tm1, K::Plane_3(-1,0,0,2));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 7 4 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 3 5\n";
    ss >> tm1;
    CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(-1,0,0,2));
    assert(vertices(tm1).size() == 6);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
    ss >> tm1;
    for (int i=0;i<3;++i)
      CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(-1,0,0,2));
    assert(vertices(tm1).size() == 7);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
    ss >> tm1;
    for (int i=0;i<3;++i)
      CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(0,1,0,0));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
    ss >> tm1;
    for (int i=0;i<3;++i)
      CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(0,-1,0,0));
    assert(vertices(tm1).size() == 7);
  }
}

int main()
{
  std::cout << "Surface Mesh" << std::endl;
  test<Surface_mesh>();

  std::cout << "Polyhedron" << std::endl;
  test<Polyhedron>();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
