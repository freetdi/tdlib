#include <fstream>
#include <string>
#include <set>
#include <boost/algorithm/string/replace.hpp>

#include "tex.cpp"
#include "helper.cpp"

#include <tdlib/combinations.hpp>
#include <tdlib/misc.hpp>

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, treedec::bag_t> TD_tree_dec_t;


double totaltime = 0;
int max_width = -2;

unsigned graphs_count = 0;
double avg_width = 0.0;
double avg_approximation = 0.0;
unsigned max_approximation = 0;

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
    gettimeofday(&t2, NULL);
    time1 = time_diff(t1, t2);

    graphs_count++;

    int width = treedec::get_width(T);

    avg_width += width;

    max_width = (width > max_width)? width : max_width;

    treedec::PP_MD(G2, T2);

    int width2 = treedec::get_width(T2);

    float approximation = 0.0;
    if(width2 != 0){
        approximation = (float)width/(float)width2;
    }
    else{
        approximation = ((float)width+1)/((float)width2+1);
    }

    max_approximation = (approximation > max_approximation)? approximation : max_approximation;
    avg_approximation += approximation;

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

    avg_approximation /= graphs_count;

    std::ofstream texresults("./Results/14_separator.tex", std::ios::app);

    tex_tabular_entry(texresults, package, graphs_count, avg_width, avg_approximation, max_width, totaltime);

    texresults.close();

    std::cout << "totaltime: " << totaltime << std::endl;
}

