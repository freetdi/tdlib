#include <fstream>
#include <string>
#include <set>
#include <boost/algorithm/string/replace.hpp>

#include "tex.cpp"
#include "helper.cpp"

#include <tdlib/combinations.hpp>
#include <tdlib/misc.hpp>

#ifdef RESULTS_HERE
#define RESULTS_DIR
#else
#define RESULTS_DIR "./Results/"
#endif

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, treedec::bag_t> TD_tree_dec_t;

double totaltime = 0;

std::stringstream convert;

void test_single(std::string filename){
    struct timeval t1, t2;
    double time1;

    TD_graph_t G;
    TD_graph_t H;

    read_dot_graph(filename, G);
    read_dot_graph(filename, H);

    TD_tree_dec_t T;

    std::cout << filename << " [" << boost::num_vertices(G) << ";" << boost::num_edges(G) << "]" << std::endl;


    gettimeofday(&t1, NULL);

    treedec::exact_decomposition_cutset(G, T);

    gettimeofday(&t2, NULL);
    time1 = time_diff(t1, t2);

    int width = treedec::get_width(T);

    if(!treedec::is_valid_treedecomposition(H, T)){
        std::cout << "invalid tree decomposition: " << filename << std::endl;
        write_dot_td("error_td.dot", T);
        exit(false);
    }

    std::string name = filename;
    boost::replace_all(name, "_", "\\_");

    convert << name << " & " << boost::num_vertices(H) << " & " << boost::num_edges(H) << " & " << width << " & " << time1 << " \\\\\n";

    totaltime += time1;

}

int main(int argc, char * const * argv){
    std::cout.setf( std::ios_base::unitbuf );

    if(argc<3){
        std::cerr << "usage\n"
                  << "    " << argv[0] << " [package] [*.dot]" << std::endl;
        exit(1);
    }

    std::string package(argv[1]);

    for(int i = 2; i < argc; i++){
        std::ifstream fin(argv[i]);
        if(!fin.is_open()){
            std::cerr << "error: could not open file " << argv[i] << std::endl;
        }
        fin.close();

        test_single(std::string(argv[i]));
    }

    std::ofstream texresults(RESULTS_DIR "7_exact_cutset.tex", std::ios::app);

    texresults << convert.str();

    texresults.close();

    std::cout << "totaltime: " << totaltime << std::endl;
}

