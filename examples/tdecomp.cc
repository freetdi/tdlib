// (c) 2016, 2017 Felix Salfelder
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//
// compute gr to td using various algorithms (in parallel)
//

#include "config.h"
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>

#ifdef USE_GALA
// #define USE_RANDOM_MD // 1
#define USE_FI
// #define USE_MSVS_TRIVIAL // 4
// #define USE_BMD // 8
// #define USE_SOME // maybe
// #define USE_EX // 32
#define USE_EX17 // 64
// pace17 heuristics
//
#endif

// these should work without gala
//#define USE_P17 // 128
#define USE_THORUP almost
//#define USE_FIPPTM


typedef enum{
	f_GR,
	f_DOT
} FF_t;

FF_t fformat=f_GR;

#include "twh.h"

struct Vertex{
    unsigned int id;
};

void graphviz(){

	boost::dynamic_properties dp(boost::ignore_other_properties);
	typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
			  Vertex> dotGt;
	dotGt p;
	// std::ifstream dot_graph(src.c_str());
	try{
		read_graphviz(std::cin, p, dp);
		// done=true;
	}catch(...){ untested();
		// BUG catch what?!
		std::cerr << "error parsing header\n";
		exit(2);
	}

	global_result = boost::num_vertices(p);
	std::cout << "c status " << global_result << " " << status_ms() << " initial\n";

	m = choose_m(m, p);

//	boost::print_graph(p);
	twh(p, m, mask);
}

void gr(){
    PARSE* p;

    try{
        p = new PARSE(std::cin); // , oUPPER);
    }catch(...){ untested();
        // BUG catch what?!
        std::cerr << "error parsing header\n";
        exit(2);
    }
    global_result = p->num_vertices();
    std::cout << "c status " << p->num_vertices() << " " << status_ms() << " initial\n";

	 m = choose_m(m, *p);

    twh(*p, m, mask);
}


int main(int argc, char * const * argv)
{
    trace1("main", getpid());
    global_result = -1u;
    // threads_running = 0;
    finished = 0;
    trace=false;

	 setup_handlers();
	 parseargs(argc, argv);

	 switch(fformat){
		 case f_GR:
			 gr();
			 break;
			case f_DOT:
			 graphviz();
	 }
    return(0);
}
