#include <fstream>
#include <string>
#include <set>
#include <boost/algorithm/string/replace.hpp>

#include "tex.cpp"
#include "helper.cpp"

#include <tdlib/preprocessing.hpp>
#include <tdlib/lower_bounds.hpp>
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

    read_dot_graph(filename, G);

    std::cout << filename << " [" << boost::num_vertices(G) << ";" << boost::num_edges(G) << "]" << std::endl;


    gettimeofday(&t1, NULL);

    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::preprocessing(G, bags);

    gettimeofday(&t2, NULL);
    time1 = time_diff(t1, t2);

    if(boost::num_edges(G) == 0){
        return;
    }

    typedef std::vector<std::set<typename boost::graph_traits<TD_graph_t>::vertex_descriptor> > components_t;
    components_t components;
    treedec::get_components(G, components);

    typename components_t::iterator i = components.begin();
    for(; i!=components.end(); ++i){
        if(i->size() == 1){
            continue;
        }

        TD_graph_t G_;
        typename std::vector<typename boost::graph_traits<TD_graph_t>::vertex_descriptor> vdMap;
        treedec::induced_subgraph(G_, G, *i, vdMap);

        double time_delta, time_delta2, time_gamma, time_deltaD, time_delta2D, time_gammaD_left;

        int delta, delta2, gamma, deltaD, delta2D, gammaD_left;

        TD_graph_t H;

        H.clear();
        boost::copy_graph(G_, H);
        gettimeofday(&t1, NULL);
        delta = treedec::lb::delta(H);
        gettimeofday(&t2, NULL);
        time_delta = time_diff(t1, t2);

        H.clear();
        boost::copy_graph(G_, H);
        gettimeofday(&t1, NULL);
        delta2 = treedec::lb::delta2(H);
        gettimeofday(&t2, NULL);
        time_delta2 = time_diff(t1, t2);

        H.clear();
        boost::copy_graph(G_, H);
        gettimeofday(&t1, NULL);
        gamma = treedec::lb::gamma(H);
        gettimeofday(&t2, NULL);
        time_gamma = time_diff(t1, t2);

        H.clear();
        boost::copy_graph(G_, H);
        gettimeofday(&t1, NULL);
        deltaD = treedec::lb::deltaD(H);
        gettimeofday(&t2, NULL);
        time_deltaD = time_diff(t1, t2);

        H.clear();
        boost::copy_graph(G_, H);
        gettimeofday(&t1, NULL);
        delta2D = treedec::lb::delta2D(H);
        gettimeofday(&t2, NULL);
        time_delta2D = time_diff(t1, t2);

        H.clear();
        boost::copy_graph(G_, H);
        gettimeofday(&t1, NULL);
        gammaD_left = treedec::lb::gammaD_left(H);
        gettimeofday(&t2, NULL);
        time_gammaD_left= time_diff(t1, t2);

        std::string name = filename;
        boost::replace_all(name, ".dot", "");
        std::string::size_type i = name.find_last_of("/");
        name = name.substr(i+1);
        boost::replace_all(name, "_", "\\_");

        convert << name << " & " << delta << " & "  << delta2 << " & " << gamma << " & " << deltaD << " & " << delta2D 
                << " & " << gammaD_left << " \\\\" << std::endl;

        totaltime += time_delta + time_delta2 + time_delta2D + time_gamma + time_deltaD + time_delta2D + time_gammaD_left;
    }

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

    std::ofstream texresults(RESULTS_DIR "2_lb1.tex", std::ios::app);

    texresults << convert.str();

    texresults.close();

    std::cout << "totaltime: " << totaltime << std::endl;
}

