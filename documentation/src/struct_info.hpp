#ifndef STRUCT_INFO
#define STRUCT_INFO

struct info{
    std::string name;
    std::string indexfile;
    std::string graphs_path;
    std::string td_dest;
    
    unsigned int graphs_count;
    int max_width;
    float avg_width;
    unsigned int not_fully_preprocessed_count;
    double total_time;

    float max_approximation;
    float avg_approximation;
    unsigned int diff_width;
    unsigned int max_diff;
    unsigned int diffs_count;
    unsigned int width_improvement_count;
};

#endif
