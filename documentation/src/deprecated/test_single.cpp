#include <fstream>
#include <string>
#include <set>
#include <boost/algorithm/string/replace.hpp>

#include "tex.cpp"
#include "helper.cpp"

#include <tdlib/preprocessing.hpp>
#include <tdlib/misc.hpp>

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, treedec::bag_t> TD_tree_dec_t;

double totaltime = 0;
int max_width = -2;

void test_single(std::string filename){
    struct timeval t1, t2;
    double time1, time2, time3;

    TD_graph_t G;
    TD_graph_t H;

    TD_tree_dec_t T1, T2;

    std::cout << filename << " [" << boost::num_vertices(G1) << ";" << boost::num_edges(G1) << "]" << std::endl;


    gettimeofday(&t1, NULL);
    gettimeofday(&t2, NULL);
    time1 = time_diff(t1, t2);

    if(!treedec::is_valid_treedecomposition(H, T){
        std::cout << "invalid tree decomposition: " << filename << std::endl;
        write_dot_td("error_td.dot", T);
        exit(false);
    }

    totaltime += time1;

}

int main(int argc, char * const * argv){
    std::cout.setf( std::ios_base::unitbuf );

    if(argc<2){
        std::cerr << "usage\n"
                  << "    " << argv[0] << " [*.dot]" << std::endl;
        exit(1);
    }

    for(int i = 1; i < argc; i++){
        std::ifstream fin(argv[i]);
        if(!fin.is_open()){
            std::cerr << "error: could not open file " << argv[i] << std::endl;
        }
        fin.close();

        test_single(std::string(argv[i]));
    }

    std::cout << "totaltime: " << totaltime << std::endl;
}









            if(algorithm == "PP"){
                gettimeofday(&t1, NULL);
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags, low);
                treedec::glue_bags(bags, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);

                if(boost::num_edges(G1) > 0){
                     //std::cout << "not fully preprocessable: " << *graph_name << std::endl;
                     //std::cout << "|V'| = " << boost::num_vertices(H) << std::endl;
                     write_dot_graph("./tmp/"+*graph_name+".dot", H);
                }

                if(boost::num_edges(G1) > 0)
                    packages[i].not_fully_preprocessed_count++;

                if(boost::num_edges(G1) == 0){
                    int status = treedec::is_valid_treedecomposition(G2, T);
                    if(status < 0)
                        std::cout << "invalid decomposition: " << *it_in << "[" << status << "]" << std::endl;
                }
            }
