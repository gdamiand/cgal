#define CGAL_NO_ASSERTIONS

//#define LOG 0
//#define LOG2 0

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

#include <tbb/global_control.h>
#include <tbb/info.h>
#include <tbb/parallel_for.h>

#define BENCHMARK_ITERATIONS 1000

using namespace CGAL;

struct Concurrent_items : public Linear_cell_complex_min_items
{
  typedef CGAL::Tag_true Use_concurrent_container;
  //typedef CGAL::Tag_true Use_char_bitset;
};

struct Concurrent_char_items : public Linear_cell_complex_min_items
{
    typedef CGAL::Tag_true Use_concurrent_container;
    typedef CGAL::Tag_true Use_char_bitset;
};

typedef Linear_cell_complex_for_combinatorial_map<3, 3, Linear_cell_complex_traits<3>, Concurrent_items> LCC_Concurrent;
typedef Linear_cell_complex_for_combinatorial_map<3, 3, Linear_cell_complex_traits<3>, Concurrent_char_items> LCC_Concurrent_char;
typedef Linear_cell_complex_for_combinatorial_map<3> LCC;

template <unsigned int d>
using One_dart_per_cell_range_concurrent = LCC_Concurrent::One_dart_per_cell_range<d>;

template <unsigned int d>
using One_dart_per_cell_range_concurrent_char = LCC_Concurrent_char::One_dart_per_cell_range<d>;

template <unsigned int d>
using One_dart_per_cell_range = LCC::One_dart_per_cell_range<d>;

void import_bunny_into_lcc(LCC_Concurrent& lcc) {
    CGAL::load_off<LCC_Concurrent>(lcc, "../data/bunny00.off");
}

void import_bunny_into_lcc(LCC_Concurrent_char& lcc) {
    CGAL::load_off<LCC_Concurrent_char>(lcc, "../data/bunny00.off");
}

void import_bunny_into_lcc(LCC& lcc) {
    CGAL::load_off<LCC>(lcc, "../data/bunny00.off");
}

void benchmark_bunny_one_dart_per_cell_iterator_concurrent()
{
    LCC_Concurrent lcc;
    import_bunny_into_lcc(lcc);

    One_dart_per_cell_range_concurrent<1> range = lcc.one_dart_per_cell<1>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
        for(One_dart_per_cell_range_concurrent<1>::iterator start = range.begin(), end = range.end();
            start != end; start++) {
            ;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();

    double total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    double per_iteration = total / (double)BENCHMARK_ITERATIONS;

    std::cout << "Total duration: " << total << "ms" << std::endl;
    std::cout << "Average duration: " << per_iteration << "ms" << std::endl;
}


void benchmark_bunny_one_dart_per_cell_iterator_char_concurrent()
{
    LCC_Concurrent_char lcc;
    import_bunny_into_lcc(lcc);

    One_dart_per_cell_range_concurrent_char<1> range = lcc.one_dart_per_cell<1>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
        for(One_dart_per_cell_range_concurrent_char<1>::iterator start = range.begin(), end = range.end();
            start != end; start++) {
            ;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();

    double total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    double per_iteration = total / (double)BENCHMARK_ITERATIONS;

    std::cout << "Total duration: " << total << "ms" << std::endl;
    std::cout << "Average duration: " << per_iteration << "ms" << std::endl;
}


void benchmark_bunny_one_dart_per_cell_iterator()
{
    LCC lcc;
    import_bunny_into_lcc(lcc);

    One_dart_per_cell_range<1> range = lcc.one_dart_per_cell<1>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
        for(One_dart_per_cell_range<1>::iterator start = range.begin(), end = range.end();
            start != end; start++){
            ;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();

    double total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    double per_iteration = total / (double)BENCHMARK_ITERATIONS;

    std::cout << "Total duration: " << total << "ms" << std::endl;
    std::cout << "Average duration: " << per_iteration << "ms" << std::endl;
}

/**
 * Parcours tous les darts du bunny mais en utilisant plusieurs threads pour que chaque thread parcours sa partie de la map
 * et qu'on ait pontentiellement un speedup
 **/
void benchmark_bunny_map_parallel_speedup() {
    LCC lcc;
    import_bunny_into_lcc(lcc);

    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < BENCHMARK_ITERATIONS; iter++) {
        for (auto i = lcc.darts().begin(), end = lcc.darts().end(); i != end; ++i)
        {
            ;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();

    double total_sequential = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    double per_iteration_sequential = total_sequential / (double)BENCHMARK_ITERATIONS;

    std::cout << "Total duration sequential: " << total_sequential << "ms" << std::endl;
    std::cout << "Average duration sequential: " << per_iteration_sequential << "ms" << std::endl << std::endl;

    LCC_Concurrent_char lcc_c;
    import_bunny_into_lcc(lcc_c);


    start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < BENCHMARK_ITERATIONS; iter++) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, lcc_c.darts().size()), [&](const tbb::blocked_range<size_t> &r) {
            for (auto i = r.begin(), end = r.end(); i != end; ++i)
            {
                ;
            }
        });
    }
    stop = std::chrono::high_resolution_clock::now();

    double total_parallel_char = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    double per_iteration_parallel_char = total_parallel_char / (double)BENCHMARK_ITERATIONS;

    std::cout << "Total duration parallel char: " << total_parallel_char << "ms" << std::endl;
    std::cout << "Average duration parallel char: " << per_iteration_parallel_char << "ms" << std::endl << std::endl;



    start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < BENCHMARK_ITERATIONS; iter++) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, lcc.darts().size()), [&](const tbb::blocked_range<size_t> &r) {
            for (auto i = r.begin(), end = r.end(); i != end; ++i)
            {
                ;
            }
        });
    }
    stop = std::chrono::high_resolution_clock::now();

    double total_parallel = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    double per_iteration_parallel = total_parallel / (double)BENCHMARK_ITERATIONS;

    std::cout << "Total duration parallel: " << total_parallel << "ms" << std::endl;
    std::cout << "Average duration parallel: " << per_iteration_parallel << "ms" << std::endl << std::endl;



    std::cout << "Parallel over sequential speedup: " << per_iteration_sequential / per_iteration_parallel << std::endl;
}

int main()
{
    std::cout << "Concurrent container: " << std::endl;
    benchmark_bunny_one_dart_per_cell_iterator_concurrent();

    std::cout << "\n";
    std::cout << "Concurrent char container: " << std::endl;
    benchmark_bunny_one_dart_per_cell_iterator_char_concurrent();

    std::cout << "\n";

    std::cout << "Non-concurrent container: " << std::endl;
    benchmark_bunny_one_dart_per_cell_iterator();

    std::cout << "\n";
    std::cout << "Parallel speedup benchmark..." << std::endl;
    benchmark_bunny_map_parallel_speedup();
}
