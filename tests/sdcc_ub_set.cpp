
#define GDIR boost::undirectedS
#define TDIR boost::bidirectionalS

#define BACKEND_GRAPH_TYPE \
  boost::adjacency_list<boost::vecS, boost::vecS, GDIR, cfg_node>

#define INPUT_GRAPH_TYPE \
  boost::adjacency_list<boost::setS, boost::vecS, GDIR, cfg_node>

// TODO
#define CONSTRUCTOR_WIP

#include "sdcc_.cpp"
