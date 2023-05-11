#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Image_3.h>
#include <CGAL/lcc_from_image3.h>

typedef float Image_word_type;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3> LCC3;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/skull_2.9.inr");
  // Load image
  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  bool simplify_vertices=false;
  bool simplify_edges=false;
  for(int i=2; i<argc; ++i)
  {
    if(argv[i]==std::string("-simplify-vertices"))
    { simplify_vertices=true; }
    else if(argv[i]==std::string("-simplify-edges"))
    { simplify_edges=true; }
    else
    { std::cout<<"Option "<<argv[i]<<" ignored."<<std::endl; }
  }
  
  LCC3 lcc;
  CGAL::lcc_from_image3(lcc, image, simplify_vertices, simplify_edges);
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
