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
int max_width = -2;

unsigned graphs_count = 0;
double avg_width = 0.0;
unsigned diff;

void test_single(std::string filename){
    struct timeval t1, t2;
    double time1;

    TD_graph_t G, G2;
    TD_graph_t H;

    read_dot_graph(filename, G);
    read_dot_graph(filename, G2);
    read_dot_graph(filename, H);

    TD_tree_dec_t T, T2;

    std::cout << filename << " [" << boost::num_vertices(G) << ";" << boost::num_edges(G) << "]" << std::endl;

    gettimeofday(&t1, NULL);
    treedec::separator_algorithm(G, T);
    int width_sep = treedec::get_width(T);
    typename std::vector<typename boost::graph_traits<TD_graph_t>::vertex_descriptor> old_elim_ordering;
    typename std::vector<typename boost::graph_traits<TD_graph_t>::vertex_descriptor> new_elim_ordering;
    treedec::treedec_to_ordering<TD_graph_t, TD_tree_dec_t>(T, old_elim_ordering);
    treedec::minimalChordal(G, old_elim_ordering, new_elim_ordering);
    T.clear();
    treedec::ordering_to_treedec(G, new_elim_ordering, T);
    int width_tm = treedec::get_width(T);
    gettimeofday(&t2, NULL);
    time1 = time_diff(t1, t2);

    avg_width += width_tm;

    max_width = (width_tm > max_width)? width_tm : max_width;

    graphs_count++;

    if(width_sep != width_tm){
        diff++;
    }

    if(!treedec::is_valid_treedecomposition(H, T)){
        std::cout << "invalid tree decomposition: " << filename << std::endl;
        write_dot_td("error_td.dot", T);
        exit(false);
    }

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

    avg_width /= graphs_count;


    std::ofstream texresults(RESULTS_DIR "19_sep_TM.tex", std::ios::app);

    tex_tabular_entry(texresults, package, graphs_count, avg_width, max_width, diff, totaltime);

    texresults.close();

    std::cout << "totaltime: " << totaltime << std::endl;
}

