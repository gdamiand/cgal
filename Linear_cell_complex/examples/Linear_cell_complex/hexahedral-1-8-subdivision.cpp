#include <CGAL/Hexahedral_subdivision.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>

struct My_items
{
  typedef CGAL::Tag_true Use_index;
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
    typedef CGAL::Cell_attribute<LCC> Vol_attrib;
    typedef std::tuple<Vertex_attrib, void, void, Vol_attrib> Attributes;
  };
};
typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_inexact_constructions_kernel> Traits;

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, Traits, My_items> LCC;
typedef LCC::Dart_descriptor Dart_descriptor;
typedef LCC::Point           Point;

Dart_descriptor random_dart(LCC& lcc)
{
  if(lcc.is_empty())
  { return lcc.null_descriptor; }
  
  Dart_descriptor res=CGAL::get_default_random().get_int(0,lcc.upper_bound_on_dart_ids());
  while(!lcc.is_dart_used(res))
  {
    ++res;
    if(res==lcc.upper_bound_on_dart_ids()) res=0;
  }
  return res;
}

int main()
{
  LCC lcc;
  Dart_descriptor d1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,5), Point(0,0,5),
                        Point(5,0,5), Point(5,5,5));
  lcc.set_attribute<3>(d1, lcc.create_attribute<3>());
  
  Hexahedral_subdivision hs(&lcc);

  for(int i=0; i<20; ++i)
  {
    hs.subdivide(random_dart(lcc));
  }
  
  lcc.display_characteristics(std::cout)<<", valid="
                                        <<lcc.is_valid()<<std::endl;
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
