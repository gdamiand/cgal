#define CGAL_NO_ASSERTIONS

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Polyhedron_3_to_lcc.h>

#include <thread>

//Nombre de threads qui vont parcourir la LCC en parallèle
#define NB_THREADS 64

//Combien de fois on va répéter l'itération sur 
//la LCC. Ce paramètre sert uniquement à avoir
//plus de chance de trouver une race condition
//en lançant le programme 1 fois. Si on a que
//une seule itération, même avec 10 threads, 
//on est pas sûr que le programme crash. 
//Si on répète 100 fois les itérations, on a 
//beaucoup plus de chance de faire crasher 
//le programme
#define ITERATIONS 1000

using namespace CGAL;

struct Concurrent_items : public Linear_cell_complex_min_items
{
  typedef CGAL::Tag_true Use_concurrent_container;
};

typedef Linear_cell_complex_for_combinatorial_map<3, 3, Linear_cell_complex_traits<3>, Concurrent_items> LCC;
typedef LCC::Point Point;

void import_bunny_into_lcc(LCC& lcc, int thread_id) {
    CGAL::load_off<LCC>(lcc, "../data/bunny00.off");
}

void import_hexaedron_into_lcc(LCC& lcc)
{
    lcc.make_hexahedron(Point(0, 0, 0), Point(1, 0, 0), Point(1, 0, 1), Point(0, 0, 1), Point(0, 1, 1), Point(0, 1, 0), Point(1, 1, 0), Point(1, 1, 1));
}

void add_hexaedron_to_lcc_parallel() {

    std::array<std::thread, NB_THREADS> threads;

    for(int iter = 0; iter < ITERATIONS; iter++)
    {
        LCC lcc;

        for(int i = 0; i < NB_THREADS; i++)
        {
            threads[i] = std::thread(import_hexaedron_into_lcc, std::ref(lcc));

            //threads[i].join();
        }

        //On join tous les threads à la fin de chaque itération
        for(int i = 0; i < NB_THREADS; i++)
        {
            if(threads[i].joinable())
                threads[i].join();
        }
    }

    for(int i = 0; i < NB_THREADS; i++)
    {
        if(threads[i].joinable())
            threads[i].join();
    }
}

int main()
{
    add_hexaedron_to_lcc_parallel();

    std::cout << "done" << std::endl;
}
