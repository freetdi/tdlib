#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <boost/algorithm/string/replace.hpp>

#include <tdlib/preprocessing.hpp>
#include <tdlib/lower_bounds.hpp>
#include <tdlib/elimination_orderings.hpp>
#include <tdlib/separator_algorithm.hpp>
#include <tdlib/combinations.hpp>
#include <tdlib/postprocessing.hpp>
#include <tdlib/misc.hpp>

#include "helper.cpp"
#include "struct_info.hpp"

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex> TD_graph_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bag> TD_tree_dec_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Thorup_graph_t;

#include "test_exact_greedy.cpp"
#include "test_exact_dynamic.cpp"
#include "Thorup.hpp"

#include "test_all_algorithms_rest_graphs.cpp"

template <typename Thorup_G_t, typename G_t>
void thorup_graph(Thorup_G_t &G2, G_t &G1){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    std::vector<typename boost::graph_traits<Thorup_G_t>::vertex_descriptor> G2vertices;
    unsigned int currentID = 0;
    while (currentID != boost::num_vertices(G1)){
        for(boost::tie(vIt, vEnd) = boost::vertices(G1); vIt != vEnd; vIt++){
            if((unsigned int) *vIt == currentID){
                G2vertices.push_back(boost::add_vertex(G2));
                currentID++;
                break;
            }
        }
    }
    for(boost::tie(vIt, vEnd) = boost::vertices(G1); vIt != vEnd; vIt++){
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G1); nIt != nEnd; nIt++){
            if(!boost::edge(G2vertices[*vIt], G2vertices[*nIt], G2).second)
                boost::add_edge(G2vertices[*vIt], G2vertices[*nIt], G2);
        }
    }
}

void testing(std::string algorithm, std::string tddir, bool cfg_graphs, bool dimacs_graphs, bool rest_graphs){
    std::vector<info> packages;
    info package;

    if(cfg_graphs){
        package.name = "stdlib"; package.indexfile = "./../Data/idx_stdlib.txt", package.graphs_path = "./../Data/TestGraphs/stdlib", 
        package.td_dest = "./Decompositions/" + tddir + "/stdlib", packages.push_back(package);
        
        package.name = "Contiki"; package.indexfile = "./../Data/idx_contiki.txt", package.graphs_path = "./../Data/TestGraphs/contiki", 
        package.td_dest = "./Decompositions/" + tddir + "/contiki", packages.push_back(package);
        
        package.name = "FUZIX"; package.indexfile = "./../Data/idx_fuzix.txt", package.graphs_path = "./../Data/TestGraphs/fuzix", 
        package.td_dest = "./Decompositions/" + tddir + "/fuzix", packages.push_back(package);
        
        package.name = "Corem."; package.indexfile = "./../Data/idx_coremark.txt", package.graphs_path = "./../Data/TestGraphs/coremark", 
        package.td_dest = "./Decompositions/" + tddir + "/coremark", packages.push_back(package);
       
        package.name = "Dhryst."; package.indexfile = "./../Data/idx_dhrystone.txt", package.graphs_path = "./../Data/TestGraphs/dhrystone", 
        package.td_dest = "./Decompositions/" + tddir + "/dhrystone", packages.push_back(package);
        
        package.name = "Whetst."; package.indexfile = "./../Data/idx_whetstone.txt", package.graphs_path = "./../Data/TestGraphs/whetstone", 
        package.td_dest = "./Decompositions/" + tddir + "/whetstone", packages.push_back(package);
    }

    if(dimacs_graphs){
        package.name = "dimacs"; package.indexfile = "./../Data/idx_dimacs.txt", package.graphs_path = "./../Data/TestGraphs/dimacs", 
        package.td_dest = "./Decompositions/" + tddir + "/dimacs", packages.push_back(package);
    }

    if(rest_graphs){
        package.name = "rest_graphs"; package.indexfile = "./../Data/idx_rest_graphs.txt"; package.graphs_path = "./../Data/TestGraphs/rest_graphs";
        package.td_dest = "./Decompositions/" + tddir + "/rest_graphs", packages.push_back(package);
    }

    std::ofstream texresults;

    if(algorithm == "PP"){
        texresults.open("./Results/1_PP.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}\n";
        texresults << "\\hline" << std::endl; 
        texresults << "package & \\#CFGs & avg width & max width & \\#not fully prep. & total time[ms] \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    if(algorithm == "PP_selection"){
        texresults.open("./Results/2_PP_selection.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|}\n";
        texresults << "\\hline" << std::endl; 
        texresults << "package & name & $|V|$ & $|E|$ & $|V'|$ & $|E'|$ & tw & avg bag size & time[ms] & fully reduced \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    if(algorithm == "LB_degree_based1"){
        texresults.open("./Results/3_LB_degree_based1.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|l|}" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "name & delta & delta2 & gamma & deltaD & delta2D & gammaD-left \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    if(algorithm == "LB_degree_based2"){
        texresults.open("./Results/4_LB_degree_based2.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "name & gammaD-right & gammaD-min-e & deltaC-min-d & deltaC-max-d & deltaC-leastC \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    if(algorithm == "LB_improved_graphs1"){
        texresults.open("./Results/5_LB_improved_graphs1.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "name & LBNdeltaD & LBNdeltaC & LBNCdeltaD & LBNCdeltaC \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    if(algorithm == "LB_improved_graphs2"){
        texresults.open("./Results/6_LB_improved_graphs2.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "name & LBPdeltaD & LBPdeltaC & LBPCdeltaD & LBPCdeltaC \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    if(algorithm == "LB_MCS"){
        texresults.open("./Results/7_LB_MCS.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "name & random & all & min-deg & last-mcs & max-mcs \\\\" << std::endl;
        texresults << "\\hline" << std::endl; 
        texresults << "\\hline" << std::endl; 
    }

    //8-9: exact algorithms

    if(algorithm == "MD"){
        texresults.open("./Results/10_MD.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "PP+MD"){
        texresults.open("./Results/11_PP_MD.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "FI"){
        texresults.open("./Results/12_FI.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "PP+FI"){
        texresults.open("./Results/13_PP_FI.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Thorup"){
        texresults.open("./Results/14_thorup.tex");
        
        texresults << "\\begin{tabular}{|l|l|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "source & \\#CFGs & avg width & max width & time[ms] & \\#$width_{T} > widht_{PP+MD}$ & max difference \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Sep_packages"){
        texresults.open("./Results/16_separator_packages.tex");
        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & avg appr. & max width & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "MSVS_trivial"){
        texresults.open("./Results/17_MSVS_trivial.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Sep+MSVS"){
        texresults.open("./Results/18_Sep_MSVS.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & \\#different & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "FI+TM"){
        texresults.open("./Results/19_FI_TM.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & \\#different & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "PP+FI+TM"){
        texresults.open("./Results/20_PP_FI_TM.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & \\#different & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Sep+TM"){
        texresults.open("./Results/21_Sep_TM.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & \\#different & total time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: Sep"){
        texresults.open("./Results/S_Sep.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: MSVS_trivial"){
        texresults.open("./Results/S_MSVS_trivial.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: Sep+MSVS"){
        texresults.open("./Results/S_Sep_MSVS.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: Sep+TM"){
        texresults.open("./Results/S_Sep_TM.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: MD"){
        texresults.open("./Results/S_MD.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: PP+MD"){
        texresults.open("./Results/S_PP_MD.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: FI"){
        texresults.open("./Results/S_FI.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: FI+TM"){
        texresults.open("./Results/S_FI_TM.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: PP+FI+TM"){
        texresults.open("./Results/S_PP_FI_TM.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    if(algorithm == "Summary: Thorup"){
        texresults.open("./Results/S_Thorup.tex");

        texresults << "\\begin{tabular}{|l|l|l|l|l|}\n";
        texresults << "\\hline\n";
        texresults << "package & \\#CFGs & avg width & max width & time[ms] \\\\\n";
        texresults << "\\hline\n";
    }

    std::cout << packages.size() << " packages: " << "";

    std::vector<std::string> data_PP_selection1, data_PP_selection2;
    std::vector<std::string> texdata;

    std::vector<std::string>::iterator it_in, it_out, graph_name;

    struct timeval t1, t2;
    TD_graph_t G1, G2, H;
    TD_tree_dec_t T, T1, T2; 

    std::vector<double> runtime_vertices(20);
    std::vector<unsigned int> count_vertices(20);
    for(unsigned int u = 0; u < 20; u++){
        runtime_vertices[u] = 0;
        count_vertices[u] = 0;
    }

    for(unsigned int i = 0; i < packages.size(); i++){
        std::cout << packages[i].name << " " << "";
        std::vector<std::string> filenames_in;
        std::vector<std::string> filenames_out;
        std::vector<std::string> graph_names;

        collect_testgraphs(packages[i].indexfile, packages[i].graphs_path, packages[i].td_dest, "", filenames_in, filenames_out, graph_names);

        it_in = filenames_in.begin();
        it_out = filenames_out.begin();
        graph_name = graph_names.begin();

        packages[i].max_width = 0;
        packages[i].avg_width = 0;
        packages[i].not_fully_preprocessed_count = 0;
        packages[i].total_time = 0;
        packages[i].graphs_count = filenames_in.size();
        packages[i].max_approximation = 0;
        packages[i].avg_approximation = 0;
        packages[i].diff_width = 0;
        packages[i].max_diff = 0;
        packages[i].diffs_count = 0;

        double running_time = 0;

        for(; it_in != filenames_in.end();){
            G1.clear();
            G2.clear();
            H.clear();
            T.clear();
            T1.clear();
            T2.clear();

            int low = -1;

            read_dot_graph(*it_in, G1);
            read_dot_graph(*it_in, G2);

            std::string name = *graph_name;
            boost::replace_all(name, "stdlib_", "");
            boost::replace_all(name, "contiki_", "");
            boost::replace_all(name, "fuzix_", "");
            boost::replace_all(name, "_", "\\_");

                    if(name == "psock\\_psock\\_readbuf_len")
                        name = "psock\\_readbuf\\_len";
                    if(name == "psock\\_psock\\_send")
                        name = "psock\\_send";
                    if(name == "shell-irc\\_process\\_thread\\_shell\\_irc\\_process")
                        name = "irc\\_process";
                    if(name == "uip\\_uip\\_process")
                        name = "uip\\_process";
                    if(name == "vfprintf\\_vfprintf")
                        name = "vfprintf";
                    if(name == "vfscanf\\_vfscanf")
                        name = "vfscanf";
                    if(name == "about\\_process\\_thread\\_about\\_process")
                        name = "thread\\_about\\_process";
                    if(name == "calc\\_process\\_thread\\_calc\\_process")
                        name = "calc\\_process";
                    if(name == "collect\\_node_packet\\_received")
                        name = "packet\\_received";
                    if(name == "ctk\\_process\\_thread\\_ctk\\_process")
                        name = "ctk\\_process";
                    if(name == "elfloader\\_elfloader\\_load")
                        name = "elfloader\\_load";
                    if(name == "shell-netperf\\_process\\_thread\\_shell\\_netperf\\_process")
                        name = "netperf\\_process";
                    if(name == "tcpdump\\_tcpdump\\_format")
                        name = "tcpdump\\_format";

            int n = boost::num_vertices(G1);
            int e = boost::num_edges(G1);

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

            if(algorithm == "PP_selection"){
                int low =-1;

                gettimeofday(&t1, NULL);
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags, low);
                treedec::glue_bags(bags, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);

                int n_ = boost::num_vertices(H);
                int e_ = boost::num_edges(H);

                float avg = 0;
                typename boost::graph_traits<TD_tree_dec_t>::vertex_iterator tIt, tEnd;
                for(boost::tie(tIt, tEnd) = boost::vertices(T);  tIt != tEnd; tIt++)
                    avg += (float)T[*tIt].bag.size();
                avg /= boost::num_vertices(T);

                std::stringstream convert;

                if(*it_in == "./../Data/TestGraphs/stdlib/stdlib_atanf.dot"
                 ||*it_in == "./../Data/TestGraphs/stdlib/stdlib_print_format.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_about_process_thread_about_process.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_calc_process_thread_calc_process.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_chameleon-raw_input.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_collect_node_packet_received.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_ctk_process_thread_ctk_process.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_cxmac_send_packet.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_elfloader_elfloader_load.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_shell-netperf_process_thread_shell_netperf_process.dot"
                 ||*it_in == "./../Data/TestGraphs/contiki/contiki_tcpdump_tcpdump_format.dot"){

                    convert << packages[i].name << " & " << name << " & " << n << " & " << e << " & " << n_ << " & " << e_ 
                            << " & " << low << " & " << avg << " & " << running_time << " & " << "yes\\\\\n";
                    data_PP_selection1.push_back(convert.str());
                }
                else if (n_ > 0){
                    convert << packages[i].name << " & " << name << " & " << n << " & " << e << " & " << n_ << " & " << e_ 
                            << " & $ \\geq $ " << low << " & " << " $>$ " << avg << " & " << running_time << " & " << "no\\\\\n";
                    data_PP_selection2.push_back(convert.str());
                }
            }

            if(algorithm == "LB_degree_based1"){
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags);

                double time_delta, time_delta2, time_gamma, time_deltaD, time_delta2D, time_gammaD_left;

                int delta, delta2, gamma, deltaD, delta2D, gammaD_left;

                gettimeofday(&t1, NULL);
                delta = treedec::lb::delta(G1);
                gettimeofday(&t2, NULL);
                time_delta = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                delta2 = treedec::lb::delta2(G1);
                gettimeofday(&t2, NULL);
                time_delta2 = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                gamma = treedec::lb::gamma(G1);
                gettimeofday(&t2, NULL);
                time_gamma = time_diff(t1, t2);
            
                gettimeofday(&t1, NULL);
                deltaD = treedec::lb::deltaD(G1);
                gettimeofday(&t2, NULL);
                time_deltaD = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                delta2D = treedec::lb::delta2D(G1);
                gettimeofday(&t2, NULL);
                time_delta2D = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                gammaD_left = treedec::lb::gammaD_left(G1);
                gettimeofday(&t2, NULL);
                time_gammaD_left= time_diff(t1, t2);

                std::stringstream convert;
                convert << name << " & " << delta << " & "  << delta2 << " & " << gamma << " & " << deltaD << " & " << delta2D 
                        << " & " << gammaD_left << " \\\\" << std::endl;
                texdata.push_back(convert.str());
            }

            if(algorithm == "LB_degree_based2"){
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags);

                double time_gammaD_right, time_gammaD_min_e, time_deltaC_min_d, time_deltaC_max_d, time_deltaC_least_c;

                int gammaD_right, gammaD_min_e, deltaC_min_d, deltaC_max_d, deltaC_least_c;

                gettimeofday(&t1, NULL);
                gammaD_right = treedec::lb::gammaD_right(G1);
                gettimeofday(&t2, NULL);
                time_gammaD_right= time_diff(t1, t2);
                
                gettimeofday(&t1, NULL);
                gammaD_min_e = treedec::lb::gammaD_min_e(G1);
                gettimeofday(&t2, NULL);
                time_gammaD_min_e= time_diff(t1, t2);
                
                gettimeofday(&t1, NULL);
                deltaC_min_d = treedec::lb::deltaC_min_d(G1);
                gettimeofday(&t2, NULL);
                time_deltaC_min_d = time_diff(t1, t2);
                
                gettimeofday(&t1, NULL);
                deltaC_max_d = treedec::lb::deltaC_max_d(G1);
                gettimeofday(&t2, NULL);
                time_deltaC_max_d = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                deltaC_least_c = treedec::lb::deltaC_least_c(G1);
                gettimeofday(&t2, NULL);
                time_deltaC_least_c = time_diff(t1, t2);

                std::stringstream convert;
                convert << name << " & " << gammaD_right << " & " << gammaD_min_e << " & " << deltaC_min_d << " & " 
                        << deltaC_max_d << " & " << deltaC_least_c << " \\\\" << std::endl;

                texdata.push_back(convert.str());
            }

            if(algorithm == "LB_improved_graphs1"){
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags);

                double time_LBNdeltaD, time_LBNdeltaC, time_LBNCdeltaD, time_LBNCdeltaC;

                int LBNdeltaD, LBNdeltaC, LBNCdeltaD, LBNCdeltaC;

                gettimeofday(&t1, NULL);
                LBNdeltaD = treedec::lb::LBN_deltaD(G1);
                gettimeofday(&t2, NULL);
                time_LBNdeltaD = time_diff(t1, t2);
                
                gettimeofday(&t1, NULL);
                LBNdeltaC = treedec::lb::LBN_deltaC(G1);
                gettimeofday(&t2, NULL);
                time_LBNdeltaC = time_diff(t1, t2);;
                
                gettimeofday(&t1, NULL);
                LBNCdeltaD = treedec::lb::LBNC_deltaD(G1);
                gettimeofday(&t2, NULL);
                time_LBNCdeltaD = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                LBNCdeltaC = treedec::lb::LBNC_deltaC(G1);
                gettimeofday(&t2, NULL);
                time_LBNCdeltaC = time_diff(t1, t2);

                std::stringstream convert;
                convert << name << " & " << LBNdeltaD << " & " << LBNdeltaC << " & "  << LBNCdeltaD << " & " << LBNCdeltaC
                        << "\\\\" << std::endl;

                texdata.push_back(convert.str());
            }

            if(algorithm == "LB_improved_graphs2"){
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags);

                double time_LBPdeltaD, time_LBPdeltaC, time_LBPCdeltaD, time_LBPCdeltaC;

                int LBPdeltaD, LBPdeltaC, LBPCdeltaD, LBPCdeltaC;

                gettimeofday(&t1, NULL);
                LBPdeltaD = treedec::lb::LBP_deltaD(G1);
                gettimeofday(&t2, NULL);
                time_LBPdeltaD = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                LBPdeltaC = treedec::lb::LBP_deltaC(G1);
                gettimeofday(&t2, NULL);
                time_LBPdeltaC = time_diff(t1, t2); 

                gettimeofday(&t1, NULL);
                LBPCdeltaD = treedec::lb::LBPC_deltaD(G1);
                gettimeofday(&t2, NULL);
                time_LBPCdeltaD = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                LBPCdeltaC = treedec::lb::LBPC_deltaC(G1);
                gettimeofday(&t2, NULL);
                time_LBPCdeltaC = time_diff(t1, t2); 

                std::stringstream convert;
                convert << name << " & " << LBPdeltaD << " & " << LBPdeltaC << " & " << LBPCdeltaD << " & " << LBPCdeltaC << "\\\\" << std::endl;

                texdata.push_back(convert.str());
            }

            if(algorithm == "LB_MCS"){
                std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
                treedec::preprocessing(G1, bags);

                double time_MCS, time_MCSC;

                int MCS, MCSC;

                MCS = -1;
                MCSC = -1;

                gettimeofday(&t1, NULL);
                MCS = treedec::lb::MCS(H);
                gettimeofday(&t2, NULL);
                time_MCS = time_diff(t1, t2);

                gettimeofday(&t1, NULL);
                MCSC = treedec::lb::MCSC(G1);
                gettimeofday(&t2, NULL);
                time_MCSC = time_diff(t1, t2);;

                std::stringstream convert;
                convert << name << " & " << MCS << " & "  << MCSC << "\\\\" << std::endl;

                texdata.push_back(convert.str());
            }

            if(algorithm == "MD"){
                gettimeofday(&t1, NULL);
                treedec::minDegree_decomp(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "PP+MD"){
                gettimeofday(&t1, NULL);
                treedec::PP_MD(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "FI"){
                gettimeofday(&t1, NULL);
                treedec::fillIn_decomp(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "PP+FI"){
                gettimeofday(&t1, NULL);
                treedec::PP_FI(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Thorup"){
                Thorup_graph_t TG;
                thorup_graph(TG, G1);

                gettimeofday(&t1, NULL);
                thorup_tree_decomposition(T, G1);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);

                treedec::minDegree_decomp(G2, T2);

                int max1 = treedec::get_width(T);
                int max2 = treedec::get_width(T2);

                if(max1-max2 > packages[i].max_diff){
                    packages[i].max_diff = max1-max2;
                }
            
                if(max1 != max2)
                    packages[i].diffs_count++;
            }

            if(algorithm == "Sep_packages"){
                if(*it_in == "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
                    it_in++;
                    it_out++;
                    graph_name++;
                    continue;
                }

                gettimeofday(&t1, NULL);
                treedec::separator_algorithm(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);

                treedec::PP_MD(G2, T2);

                int max1 = treedec::get_width(T);
                int max2 = treedec::get_width(T2);
                float approximation = 0.0;

                if(max2 != 0)
                    approximation = (float)max1/(float)max2;
                else
                    approximation = ((float)max1+1)/((float)max2+1);

                packages[i].max_approximation = (approximation > packages[i].max_approximation)? approximation : packages[i].max_approximation;
                packages[i].avg_approximation += approximation;
            }

            if(algorithm == "MSVS_trivial"){
                gettimeofday(&t1, NULL);
                treedec::trivial_decomposition(G1, T);
                treedec::MSVS(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Sep+MSVS"){
                if(*it_in == "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
                    it_in++;
                    it_out++;
                    graph_name++;
                    continue;
                }

                gettimeofday(&t1, NULL);
                treedec::separator_algorithm(G1, T);
                int width_sep = treedec::get_width(T);
                treedec::MSVS(G1, T);
                int width_msvs = treedec::get_width(T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);


                if(width_sep != width_msvs){
                    packages[i].diff_width++;
                }
            }

            if(algorithm == "FI+TM"){
                treedec::fillIn_decomp(G1, T1);

                gettimeofday(&t1, NULL);
                treedec::FI_TM(G2, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);


                if(treedec::get_width(T1) != treedec::get_width(T)){
                    //std::cout << *it_in << std::endl;
                    packages[i].diff_width++;
                }
            }

            if(algorithm == "PP+FI+TM"){
                treedec::PP_FI(G1, T1);

                gettimeofday(&t1, NULL);
                treedec::PP_FI_TM(G2, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);


                if(treedec::get_width(T1) != treedec::get_width(T)){
                    //std::cout << *it_in << std::endl;
                    packages[i].diff_width++;
                }
            }

            if(algorithm == "Sep+TM"){
                if(*it_in == "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
                    it_in++;
                    it_out++;
                    graph_name++;
                    continue;
                }

                gettimeofday(&t1, NULL);
                treedec::separator_algorithm(G1, T);
                int width_sep = treedec::get_width(T);
                typename std::vector<typename boost::graph_traits<TD_graph_t>::vertex_descriptor> old_elim_ordering;
                typename std::vector<typename boost::graph_traits<TD_graph_t>::vertex_descriptor> new_elim_ordering;
                treedec::treedec_to_ordering<TD_graph_t, TD_tree_dec_t>(T, old_elim_ordering);
                treedec::minimalChordal(G1, old_elim_ordering, new_elim_ordering);
                T.clear();
                treedec::ordering_to_treedec(G1, new_elim_ordering, T);
                int width_tm = treedec::get_width(T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);


                if(width_sep != width_tm){
                    //std::cout << *it_in << std::endl;
                    packages[i].diff_width++;
                }
            }

            if(algorithm == "Summary: Sep"){
                if(*it_in == "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
                    it_in++;
                    it_out++;
                    graph_name++;
                    continue;
                }
                gettimeofday(&t1, NULL);
                treedec::separator_algorithm(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: MSVS_trivial"){
                gettimeofday(&t1, NULL);
                treedec::MSVS_trivial(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: Sep+MSVS"){
                if(*it_in == "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
                    it_in++;
                    it_out++;
                    graph_name++;
                    continue;
                }
                gettimeofday(&t1, NULL);
                treedec::separator_algorithm_MSVS(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: Sep+TM"){
                if(*it_in == "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
                    it_in++;
                    it_out++;
                    graph_name++;
                    continue;
                }
                gettimeofday(&t1, NULL);
                treedec::separator_algorithm_TM(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: MD"){
                gettimeofday(&t1, NULL);
                treedec::minDegree_decomp(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: PP+MD"){
                gettimeofday(&t1, NULL);
                treedec::PP_MD(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: FI"){
                gettimeofday(&t1, NULL);
                treedec::fillIn_decomp(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: FI+TM"){
                gettimeofday(&t1, NULL);
                treedec::FI_TM(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: PP+FI+TM"){
                gettimeofday(&t1, NULL);
                treedec::PP_FI_TM(G1, T);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(algorithm == "Summary: Thorup"){
                Thorup_graph_t TG;
                thorup_graph(TG, G1);

                gettimeofday(&t1, NULL);
                thorup_tree_decomposition(T, TG);
                gettimeofday(&t2, NULL);
                running_time = time_diff(t1, t2);
            }

            if(tddir != ""){
                write_dot_td(*it_out, T);
            }

            packages[i].total_time += running_time;

            int width = treedec::get_width(T);

            if(width > 0){
                if(width > packages[i].max_width)
                    packages[i].max_width = width;
            }

            if(width > 0){
                packages[i].avg_width += width;
            }

            it_in++;
            it_out++;
            graph_name++;
        }
        packages[i].avg_approximation /= packages[i].graphs_count;
        packages[i].avg_width /= packages[i].graphs_count;
    }

    if(algorithm == "PP"){ 
        for(unsigned int i = 0; i < packages.size(); i++){
            texresults << packages[i].name << " & " << packages[i].graphs_count << " & " << packages[i].avg_width << " & " 
                       << packages[i].max_width << " & " << packages[i].not_fully_preprocessed_count << " & " 
                       <<  packages[i].total_time << " \\\\\n";
        }

        texresults << "\\hline" << std::endl; 
        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }

    if(algorithm == "PP_selection"){
        for(unsigned int i = 0; i < data_PP_selection1.size(); i++)
            texresults << data_PP_selection1[i];

        texresults << "\\hline" << std::endl;

        for(unsigned int i = 0; i < data_PP_selection2.size(); i++)
            texresults << data_PP_selection2[i];

        texresults << "\\hline" << std::endl; 
        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }

    if(algorithm == "LB_degree_based1" || algorithm == "LB_degree_based2" 
    || algorithm == "LB_improved_graphs1" || algorithm == "LB_improved_graphs2" 
    || algorithm == "LB_MCS"){
        for(unsigned int i = 0; i < texdata.size(); i++)
            texresults << texdata[i];

        texresults << "\\hline" << std::endl;

        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }

    if(algorithm == "MD" || algorithm == "PP+MD" || algorithm == "FI" || algorithm == "PP+FI" || algorithm == "MSVS_trivial"
    || algorithm == "Summary: Sep" || algorithm == "Summary: MSVS_trivial" || algorithm == "Summary: Sep+MSVS" || algorithm == "Summary: Sep+TM"
    || algorithm == "Summary: MD" || algorithm == "Summary: PP+MD" || algorithm == "Summary: FI" || algorithm == "Summary: FI+TM"
    || algorithm == "Summary: PP+FI+TM" || algorithm == "Summary: Thorup"){
        for(unsigned int i = 0; i < packages.size(); i++){
            texresults << packages[i].name << " & " << packages[i].graphs_count << " & " << packages[i].avg_width << " & " 
                       << packages[i].max_width  << " & " <<  packages[i].total_time << " \\\\\n";
        }

        texresults << "\\hline" << std::endl; 
        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }

    if(algorithm == "Thorup"){
        for(unsigned int i = 0; i < packages.size(); i++){
            texresults << packages[i].name << " & " << packages[i].graphs_count << " & " << packages[i].avg_width << " & " << packages[i].max_width
                       << " & " << packages[i].total_time << " & " << packages[i].diffs_count << " & " 
                       << packages[i].max_diff << " \\\\\n";
        }

        texresults << "\\hline" << std::endl; 
        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }

    if(algorithm == "Sep_packages"){
        for(unsigned int i = 0; i < packages.size(); i++){
            texresults << packages[i].name << " & " << packages[i].graphs_count << " & " << packages[i].avg_width << " & "
                       << packages[i].avg_approximation << " & " << packages[i].max_width << " & " << packages[i].total_time << " \\\\\n";
        }
        texresults << "\\hline" << std::endl; 
        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }

    if(algorithm == "Sep+MSVS" || algorithm == "FI+TM" || algorithm == "PP+FI+TM" || algorithm == "Sep+TM"){ 
        for(unsigned int i = 0; i < packages.size(); i++){
            texresults << packages[i].name << " & " << packages[i].graphs_count << " & " << packages[i].avg_width << " & " 
                       << packages[i].max_width << " & " << packages[i].diff_width << " & " 
                       <<  packages[i].total_time << " \\\\\n";
        }

        texresults << "\\hline" << std::endl; 
        texresults << "\\end{tabular}" << std::endl;
        texresults.close();
    }
}

int main(){
    std::cout.setf( std::ios_base::unitbuf );

    std::vector<bool> apply(30);
    apply[ 0] = false; //none
    apply[ 1] = false; //PP
    apply[ 2] = false; //PP_selection
    apply[ 3] = false; //LB_degree_based1
    apply[ 4] = false; //LB_degree_based2
    apply[ 5] = false; //LB_improved_graphs1
    apply[ 6] = false; //LB_improved_graphs2
    apply[ 7] = false; //LB_MCS
    apply[ 8] = false; //exact_greedy
    apply[ 9] = false; //exact_dynamic
    apply[10] = false; //MD
    apply[11] = false; //PP+MD
    apply[12] = false; //FI
    apply[13] = false; //PP+FI
    apply[14] = false; //Thorup
    apply[15] = true; //Sep_runtime
    apply[16] = true; //Sep_packages
    apply[17] = false; //MSVS_trivial
    apply[18] = true; //Sep+MSVS
    apply[19] = false; //FI+TM
    apply[20] = false; //PP+FI+TM
    apply[21] = true; //Sep+TM
    apply[22] = false; //Summary
    apply[23] = false; //all_algorithms_rest_graphs
    apply[24] = false; //maxClique
    apply[25] = false; //maxIndependentSet
    apply[26] = false; //minVertexCover
    apply[27] = false; //minDominatingSet
    apply[28] = false; //minColoring

    std::cout << "( 1/25) PP, " << "";
    if(apply[1]){ testing("PP", "PP", true, true, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 2/25) PP_selection, " << "";
    if(apply[2]){ testing("PP_selection", "", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 3/25) LB_degree_based1, " << "";
    if(apply[3]){ testing("LB_degree_based1", "", false, false, true); }
    else { std::cout << "skipped"; }
    std::cout << std::endl;

    std::cout << "( 4/25) LB_degree_based2, " << "";
    if(apply[4]){ testing("LB_degree_based2", "", false, false, true); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 5/25) LB_improved_graphs1, " << "";
    if(apply[5]){ testing("LB_improved_graphs1", "", false, false, true); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 6/25) LB_improved_graphs2, " << "";
    if(apply[6]){ testing("LB_improved_graphs2", "", false, false, true); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 7/25) LB_MCS, " << "";
    if(apply[7]){ testing("LB_MCS", "", false, false, true); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 8/25) exact_greedy, ";
    if(apply[8]){ test_exact_greedy(); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "( 9/25) exact_dynamic, ";
    if(apply[9]){ //test_exact_dynamic();
    }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(10/25) MD, " << "";
    if(apply[10]){ testing("MD", "MD", true, true, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(11/25) PP+MD, " << "";
    if(apply[11]){ testing("PP+MD", "PP_MD", true, true, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(12/25) FI, " << "";
    if(apply[12]){ testing("FI", "FI", true, true, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(13/25) PP+FI, " << "";
    if(apply[13]){ testing("PP+FI", "PP_FI", true, true, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(14/25) Thorup " << "";
    if(apply[14]){ testing("Thorup", "Thorup", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(15/25) Sep_runtime, ";
    if(apply[15]){ testing("Sep_runtime", "", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(16/25) Sep_packages, ";
    if(apply[16]){ testing("Sep_packages", "separator_algorithm", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(17/25) MSVS_trivial, ";
    if(apply[17]){ testing("MSVS_trivial", "MSVS_trivial", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(18/25) Sep+MSVS, ";
    if(apply[18]){ testing("Sep+MSVS", "Sep_MSVS", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(19/25) FI+TM, ";
    if(apply[19]){ testing("FI+TM", "FI_TM", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(20/25) PP+FI+TM, ";
    if(apply[20]){ testing("PP+FI+TM", "PP_FI_TM", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(21/25) Sep+TM, ";
    if(apply[21]){ testing("Sep+TM", "Sep_TM", true, false, false); }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(22/25) Summary" << std::endl;
    if(apply[22]){
        std::cout << "    ( 1/10) Sep, ";
        testing("Summary: Sep", "", true, false, false); 
        std::cout << std::endl;
        std::cout << "    ( 2/10) MSVS_trivial, ";
        testing("Summary: MSVS_trivial", "", true, false, false); 
        std::cout << std::endl;
        std::cout << "    ( 3/10) Sep+MSVS, ";
        testing("Summary: Sep+MSVS", "", true, false, false); 
        std::cout << std::endl;
        std::cout << "    ( 4/10) Sep+TM, ";
        testing("Summary: Sep+TM", "", true, false, false);
        std::cout << std::endl;
        std::cout << "    ( 5/10) MD, ";
        testing("Summary: MD", "", true, false, false);
        std::cout << std::endl;
        std::cout << "    ( 6/10) PP+MD, ";
        testing("Summary: PP+MD", "", true, false, false);
        std::cout << std::endl;
        std::cout << "    ( 7/10) FI, ";
        testing("Summary: FI", "", true, false, false);
        std::cout << std::endl;
        std::cout << "    ( 8/10) FI+TM, ";
        testing("Summary: FI+TM", "", true, false, false);
        std::cout << std::endl;
        std::cout << "    ( 9/10) PP+FI+TM, ";
        testing("Summary: PP+FI+TM", "", true, false, false);
        std::cout << std::endl;
        std::cout << "    (10/10) Thorup, ";
        testing("Summary: Thorup", "", true, false, false);
        std::cout << std::endl;
    }
    else { std::cout << "skipped"; }
    std::cout << "" << std::endl;

    std::cout << "(25/28) all_algorithms_rest_graphs ";
    if(apply[25]){ test_all_algorithms_rest_graphs(); }
    else { std::cout << "skipped"; }

    std::cout << "" << std::endl;

    std::cout << "(24/28) max_clique, " << "not implemented yet";
    std::cout << "" << std::endl;

    std::cout << "(25/28) max_independent_set, " << "not implemented yet";
    std::cout << "" << std::endl;

    std::cout << "(26/28) min_vertex_cover, " << "not implemented yet";
    std::cout << "" << std::endl;

    std::cout << "(27/28) min_dominating_set, " << "not implemented yet";
    std::cout << "" << std::endl;

    std::cout << "(28/28) min_coloring, " << "not implemented yet";
    std::cout << "" << std::endl;
}


