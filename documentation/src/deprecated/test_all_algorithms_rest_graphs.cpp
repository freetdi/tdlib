
void test_all_algorithms_rest_graphs(){
    std::vector<std::string> graphs;
    std::vector<int> tws;

    graphs.push_back("./../Data/TestGraphs/contiki/contiki_contikimac_powercycle.dot"); tws.push_back(5);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_dhcpc_handle_dhcp.dot"); tws.push_back(6);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_ircc_handle_connection.dot");tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_ircc_handle_input.dot"); tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_lpp_dutycycle.dot"); tws.push_back(5);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_psock_psock_readbuf_len.dot");  tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_psock_psock_send.dot"); tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_shell-irc_process_thread_shell_irc_process.dot"); tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot");  tws.push_back(7);
    graphs.push_back("./../Data/TestGraphs/fuzix/fuzix_vfprintf_vfprintf.dot"); tws.push_back(6);
    graphs.push_back("./../Data/TestGraphs/fuzix/fuzix_vfscanf_vfscanf.dot"); tws.push_back(6);

    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_contikimac_powercycle.dot"); tws.push_back(5);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_dhcpc_handle_dhcp.dot"); tws.push_back(6);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_ircc_handle_connection.dot"); tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_ircc_handle_input.dot");tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_lpp_dutycycle.dot"); tws.push_back(5);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_psock_psock_readbuf_len.dot"); tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_psock_psock_send.dot");  tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_shell-irc_process_thread_shell_irc_process.dot");  tws.push_back(4);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/contiki_uip_uip_process.dot"); tws.push_back(7);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/fuzix_vfprintf_vfprintf.dot");  tws.push_back(6);
    graphs.push_back("./../Data/TestGraphs/rest_graphs/fuzix_vfscanf_vfscanf.dot"); tws.push_back(6);

    struct timeval t1, t2;  

    std::ofstream texresults1("./Results/25_all_algorithms_rest_graphs_1.tex");
    std::ofstream texresults2("./Results/25_all_algorithms_rest_graphs_2.tex");

    texresults1 << "\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|}\n";
    texresults1 << "\\hline\n";
    texresults1 << "name & tw & MD & PP+MD & FI & FI+TM & MSVS & Sep. & Sep.+MSVS & Sep.+TM \\\\\n";
    texresults1 << "\\hline\n";

    texresults2 << "\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|}\n";
    texresults2 << "\\hline\n";
    texresults2 << "name & tw & MD & PP+MD & FI & FI+TM & MSVS & Sep. & Sep.+MSVS & Sep.+TM \\\\\n";
    texresults2 << "\\hline\n";

    TD_graph_t G1,G2,G3,G4,G5,G6,G7,G8,G9;
    boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> H;
    TD_tree_dec_t T1,T2,T3,T4,T5,T6,T7,T8,T9;

    for(unsigned int i = 0; i < graphs.size(); i++){
        H.clear();
        T1.clear();
        T2.clear();
        T3.clear();
        T4.clear();
        T5.clear();
        T6.clear();
        T7.clear();
        T8.clear();
        T9.clear();

        G1.clear();
        G2.clear();
        G3.clear();
        G4.clear();
        G5.clear();
        G6.clear();
        G7.clear();
        G8.clear();
        G9.clear();

        read_dot_graph(graphs[i], G1);

        std::cout << graphs[i] << std::endl;

        if(i >= 11){
            int low = -1;
            std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
            treedec::preprocessing(G1, bags, low);
            bags.clear();
        }

        boost::copy_graph(G1, G2);
        boost::copy_graph(G1, G3);
        boost::copy_graph(G1, G4);
        boost::copy_graph(G1, G5);
        boost::copy_graph(G1, G6);
        boost::copy_graph(G1, G7);
        boost::copy_graph(G1, G8);
        boost::copy_graph(G1, G9);


        int n = boost::num_vertices(G1);
        int e = boost::num_edges(G1);

        int width_MD = -2;
        int width_PP_MD = -2; 
        int width_FI = -2;
        int width_FI_TM = -2;
        int width_MSVS = -2;
        int width_Sep = -2;
        int width_Sep_MSVS = -2;
        int width_Sep_TM = -2;

        //MD
        treedec::minDegree_decomp(G1, T1);
        width_MD = treedec::get_width(T1);

        //PP+MD
        treedec::PP_MD(G2, T2);
        width_PP_MD = treedec::get_width(T2);

        //FI
        treedec::fillIn_decomp(G3, T3);
        width_FI = treedec::get_width(T3);

        //FI+TM
        treedec::FI_TM(G4, T4);
        width_FI_TM = treedec::get_width(T4);

        //MSVS trivial
        treedec::MSVS_trivial(G6, T6);
        width_MSVS = treedec::get_width(T6);

        if(graphs[i] != "./../Data/TestGraphs/contiki/contiki_uip_uip_process.dot"){
            //Separator
            treedec::separator_algorithm(G7, T7);
            width_Sep = treedec::get_width(T7);

            //Separator+MSVS
            treedec::separator_algorithm_MSVS(G8, T8);
            width_Sep_MSVS = treedec::get_width(T8);

            //Separator+TM
            treedec::separator_algorithm_TM(G9, T9);
            width_Sep_TM = treedec::get_width(T9);
        }

        boost::replace_all(graphs[i], "./../Data/TestGraphs/contiki/contiki_", "");
        boost::replace_all(graphs[i], "./../Data/TestGraphs/fuzix/fuzix_", "");
        boost::replace_all(graphs[i], "./../Data/TestGraphs/rest_graphs/contiki_", "");
        boost::replace_all(graphs[i], "./../Data/TestGraphs/rest_graphs/fuzix_", "");
        boost::replace_all(graphs[i], ".dot", "");
        boost::replace_all(graphs[i], "_", "\\_");

        if(i < 11){
            texresults1 << graphs[i] << " & " << tws[i] << " & ";

            if(width_MD == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_MD << "} & ";
            else
                texresults1 << width_MD << " & ";

            if(width_PP_MD == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_PP_MD << "} & ";
            else
                texresults1 << width_PP_MD << " & ";

            if(width_FI == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_FI << "} & ";
            else
                texresults1 << width_FI << " & ";

            if(width_FI_TM == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_FI_TM << "} & ";
            else
                texresults1 << width_FI_TM << " & ";

            if(width_MSVS == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_MSVS << "} & ";
            else
                texresults1 << width_MSVS << " & ";

            if(width_Sep == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_Sep << "} & ";
            else if(width_Sep != -2)
                texresults1 << width_Sep << " & ";
            else
                texresults1 << " - " << " & ";

            if(width_Sep_MSVS == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_Sep_MSVS << "} & ";
            else if(width_Sep_MSVS != -2)
                texresults1 << width_Sep_MSVS << " & ";
            else
                texresults1 << " - " << " & ";

            if(width_Sep_TM == tws[i]) 
                texresults1 << "\\textcolor{cgreen}{" << width_Sep_TM << "}" << " \\\\\n"; 
            else if(width_Sep_TM != -2)
                texresults1 << width_Sep_TM << " \\\\\n"; 
            else
                texresults1 << " - " << " \\\\\n";
        }
        else{
            texresults2 << graphs[i] << " & " << tws[i] << " & ";

            if(width_MD == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_MD << "} & ";
            else
                texresults2 << width_MD << " & ";

            if(width_PP_MD == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_PP_MD << "} & ";
            else
                texresults2 << width_PP_MD << " & ";

            if(width_FI == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_FI << "} & ";
            else
                texresults2 << width_FI << " & ";

            if(width_FI_TM == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_FI_TM << "} & ";
            else
                texresults2 << width_FI_TM << " & ";

            if(width_MSVS == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_MSVS << "} & ";
            else
                texresults2 << width_MSVS << " & ";

            if(width_Sep == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_Sep << "} & ";
            else if(width_Sep != -2)
                texresults2 << width_Sep << " & ";

            if(width_Sep_MSVS == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_Sep_MSVS << "} & ";
            else if(width_Sep_MSVS != -2)
                texresults2 << width_Sep_MSVS << " & ";

            if(width_Sep_TM == tws[i]) 
                texresults2 << "\\textcolor{cgreen}{" << width_Sep_TM << "}" << " \\\\\n"; 
            else if(width_Sep_TM != -2)
                texresults2 << width_Sep_TM << " \\\\\n"; 
        }
    }

    texresults1 << "\\hline\n";
    texresults1 << "\\end{tabular}\n";

    texresults2 << "\\hline\n";
    texresults2 << "\\end{tabular}\n";

    texresults1.close();
    texresults2.close();
}
