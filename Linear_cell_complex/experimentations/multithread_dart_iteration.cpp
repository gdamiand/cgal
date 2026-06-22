//#define CGAL_NO_ASSERTIONS

//Oui ou non activer les messages de trace du code
#define LOG 0
#define LOG2 0

#ifdef LOG
    #if LOG
        #define LOG_MESSAGE(x) x;
    #else
        #define LOG_MESSAGE(x)
    #endif
#endif

#ifdef LOG2
    #if LOG2
        #define LOG_MESSAGE2(x) x
    #else
        #define LOG_MESSAGE2(x)
    #endif
#endif


#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Polyhedron_3_to_lcc.h>

#include <thread>

//How many threads in parallel will iterate over the LCC
#define NB_THREADS 15

//How many times we're going to repeat the iteration
//of all the threads over the same LCC. This increases
//the chance of a concurrency between threads if the program
//isn't thread-safe
#define ITERATIONS 10000

using namespace CGAL;

struct ConcurrentItems : Linear_cell_complex_min_items
{
    typedef CGAL::Tag_true Use_concurrent_container;
};

//typedef Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef Linear_cell_complex_for_combinatorial_map<3, 3, Linear_cell_complex_traits<3>, ConcurrentItems> LCC;

template <unsigned int d>
using One_dart_per_cell_range = LCC::One_dart_per_cell_range<d>;

typedef LCC::Dart_range Dart_range;
typedef LCC::Dart_handle Dart_handle;
typedef LCC::Dart_descriptor Dart_descriptor;
typedef LCC::Dart Dart;
typedef LCC::Point Point;

void iterate_over_darts(LCC& lcc) {
    One_dart_per_cell_range<1> vertices_range = lcc.one_dart_per_cell<1>();

    for(typename One_dart_per_cell_range<1>::iterator
        start = vertices_range.begin();
        start.cont();
        start++)
    {
        ;
    }
}

std::atomic<unsigned int> counter;
std::array<std::stringstream, NB_THREADS * ITERATIONS> output_streams;
void iterate_over_darts_concurrent(LCC& lcc) {
    LCC::One_dart_per_cell_range dart_range = lcc.one_dart_per_cell<1>();

    unsigned int stream_number = counter++;

    for(typename LCC::One_dart_per_cell_range<1>::iterator
        start = dart_range.begin(), end = dart_range.end();
        start != end;
        ++start)
    {
        output_streams[stream_number] << lcc.point((*start).get_f(0)) << std::endl;
    }
}

template <typename LCC>
void import_bunny_into_lcc(LCC& lcc) {
    CGAL::load_off<LCC>(lcc, "../data/bunny00.off");
}

void iterator_compact_container() {
    //Tableau qui va stocker tous les handles des threads qu'on va créer
    std::array<std::thread, NB_THREADS> threads;

    LCC lcc;
    import_bunny_into_lcc(lcc);

    //En lançant deux threads qui vont itérer sur les même darts,
    //on va avoir des problèmes de concurrence sur la gestion des
    //marks en interne et ça va exploser
    for(int iter = 0; iter < ITERATIONS; iter++)
    {
        for(int i = 0; i < NB_THREADS; i++)
        {
            threads[i] = std::thread(iterate_over_darts, std::ref(lcc));

            //Si on décommente cette ligne, les deux
            //threads ne vont pas s'exécuter en parallèle
            //(puisque le main thread va attendre la fin
            //du thread 1 avant de lancer le deuxième thread)
            //et donc on ne va pas avoir de race condition

            //threads[i].join();
        }

        //On join tous les threads à la fin de chaque itération
        for(int i = 0; i < NB_THREADS; i++)
        {
            if(threads[i].joinable())
                threads[i].join();
        }
    }

    //On va attendre la fin de tous les
    //threads pour pas que le main thread
    //finisse son exécution et donc destroy
    //la LCC sur laquelle les threads sont
    //en train d'itérer
    for(int i = 0; i < NB_THREADS; i++)
    {
        if(threads[i].joinable())
            threads[i].join();
    }
}

void iterator_concurrent_compact_container() {
    //Tableau qui va stocker tous les handles des threads qu'on va créer
    std::array<std::thread, NB_THREADS> threads;

    LCC lcc;
    //import_bunny_into_lcc(lcc);
    lcc.make_hexahedron(Point(0, 0, 0), Point(1, 0, 0), Point(1, 0, 1), Point(0, 0, 1), Point(0, 1, 1), Point(0, 1, 0), Point(1, 1, 0), Point(1, 1, 1));

    for(int iter = 0; iter < ITERATIONS; iter++)
    {
        for(int i = 0; i < NB_THREADS; i++)
        {
            threads[i] = std::thread(iterate_over_darts_concurrent, std::ref(lcc));

            //threads[i].join();
        }

        //On join tous les threads à la fin de chaque itération
        for(int i = 0; i < NB_THREADS; i++)
        {
            if(threads[i].joinable())
                threads[i].join();
        }

        LOG_MESSAGE(std::cout << "\n\n\n" << std::endl);
    }

    for(int i = 0; i < NB_THREADS; i++)
    {
        if(threads[i].joinable())
            threads[i].join();
    }
}

/**
 * Checks whether or not the threads produced the same output
 * while iterating over the LCC in parallel
 */
void check_output_equality() {
    std::unordered_set<std::size_t> hash_set;

    for(int i = 0; i < NB_THREADS * ITERATIONS; i++) {
        hash_set.insert(std::hash<std::string>{}(output_streams[i].str()));
    }

    //If the hash set only contains 1 hash, this means that all the hashes
    //of the files were the same
    std::cout << (hash_set.size() == 1 ? "SUCCESS, all threads output the same result" : "FAILURE, all threads did not output the same result") << std::endl;
}

int main(int argc, char* argv[])
{
    //iterator_compact_container();
    iterator_concurrent_compact_container();
    check_output_equality();
}
