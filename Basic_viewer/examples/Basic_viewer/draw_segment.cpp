#include <CGAL/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using Point=CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;

int main()
{
  CGAL::Graphics_scene gs;
  gs.add_segment(Point(0.,0.,0.), Point(1.,1.,1.));
  CGAL::draw_graphics_scene(gs);
  return EXIT_SUCCESS;
}
