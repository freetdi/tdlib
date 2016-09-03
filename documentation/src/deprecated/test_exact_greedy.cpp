#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <boost/thread.hpp>

#include <tdlib/combinations.hpp>
#include <tdlib/misc.hpp>

#include "helper.cpp"

#ifndef MUTEX
#define MUTEX
boost::mutex m;
#endif

void test_parallel_instance_greedy(std::string fin, std::string fout, std::string name, std::vector<bool> &job_finished, std::vector<bool> &threads_finished, unsigned int idx1, unsigned int idx2){
    struct timeval t1, t2;
    TD_graph_t G,H;
    TD_tree_dec_t T;

    read_dot_graph(fin, G);
    read_dot_graph(fin, H);

    unsigned int n = boost::num_vertices(G);
    unsigned int e = boost::num_edges(G);

    gettimeofday(&t1, NULL);
    treedec::exact_decomposition_cutset(G, T);
    gettimeofday(&t2, NULL);
    double running_time = time_diff(t1, t2);

    int status = treedec::is_valid_treedecomposition(H, T);
    if(status < 0)
        std::cout << "[!] [" << name << "] invalid decomposition: " << status << std::endl;
      
    std::cout << "[" << name << "] width  : " << treedec::get_width(T) << std::endl;
    std::cout << "[" << name << "] time[s]: " << running_time/1000 << std::endl;

    m.lock();
    write_dot_td(fout, T);

    std::ofstream texresults("./Results/8_exact_greedy.tex", std::ios::app);
    texresults << name << " & " << n << " & " << e << " & " << treedec::get_width(T) << " & " << running_time << " \\\\\n";
    texresults.close();

    job_finished[idx1] = true;
    threads_finished[idx2] = true;
    m.unlock();

    std::cout << "-[finishing thread: " << "slot " <<  idx2 << ", task " << name << "]" << std::endl;
}

void test_exact_greedy(){
    std::vector<info> packages;
    info package;

    std::ofstream texresults("./Results/8_exact_greedy.tex", std::ios::app);
    texresults << "\\begin{tabular}{|l|l|l|l|l|l|l|}\n";
    texresults << "\\hline" << std::endl; 
    texresults << "name & $|V|$ & $|E|$ & width & time_{computation}[ms] \\\\" << std::endl;
    texresults << "\\hline" << std::endl; 
    texresults.close();
    
    package.name = "rest_graphs"; package.indexfile = "./../Data/idx_rest_graphs.txt"; package.graphs_path = "./../Data/TestGraphs/rest_graphs";
    package.td_dest = "./Decompositions/exact_greedy/"; packages.push_back(package);
    
    std::vector<std::string>::iterator it_in, it_out, it_name;
    
    for(int i = 0; i < packages.size(); i++){
        std::vector<std::string> filenames_in;
        std::vector<std::string> filenames_out;
        std::vector<std::string> graph_names;
        
        collect_testgraphs(packages[i].indexfile, packages[i].graphs_path, packages[i].td_dest, "", filenames_in, filenames_out, graph_names);
        
        it_in = filenames_in.begin();
        it_out = filenames_out.begin();
        it_name = graph_names.begin();

        unsigned int cores = boost::thread::hardware_concurrency();

        unsigned int tasks_count = filenames_in.size();

        cores = (tasks_count < cores)? tasks_count : cores;

        std::vector<bool> threads_finished(cores, true);
        std::vector<bool> tasks_finished(tasks_count, false);
        bool all_scheduled = false;
        unsigned int job = 0;

        while(std::find(tasks_finished.begin(), tasks_finished.end(), false) != tasks_finished.end()){
            std::vector<bool>::iterator it = std::find(threads_finished.begin(), threads_finished.end(), true);
            if(it == threads_finished.end() || all_scheduled){
                //all slots are occupied
                usleep(1000);
                continue;
            }
            else{
                unsigned int pos = it - threads_finished.begin();
                std::cout << "+[starting thread: " << "slot " <<  pos << ", task " << *it_in << "]" << std::endl;

                m.lock();
                threads_finished[pos] = false;
                m.unlock();

                boost::thread t(test_parallel_instance_greedy, *it_in, *it_out, *it_name, boost::ref(tasks_finished), boost::ref(threads_finished), job, pos);
                it_in++;
                it_out++;
                it_name++;
                job++;
                if(job == tasks_count)
                    all_scheduled = true;
            }
        }
    }

    texresults.open("./Results/8_exact_greedy.tex", std::ios::app);
    texresults << "\\end{tabular}\n"; 
}
