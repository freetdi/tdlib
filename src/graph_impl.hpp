// Felix Salfelder, 2015 - 2016
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
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
#ifndef TD_GRAPH_IMPL_H
#define TD_GRAPH_IMPL_H
#include "graph_traits.hpp"

namespace treedec{

namespace impl{

template<typename B, typename E, typename G_t>
void make_clique(B nIt1, E nEnd, G_t &G, typename treedec::graph_callback<G_t>* cb=NULL)
{
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    B nIt2;
    for(; nIt1!=nEnd; ++nIt1){
#if 1
        if(cb){
            (*cb)(*nIt1); // hmm
        }
#endif
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){
            std::pair<edge_descriptor, bool> ep=boost::add_edge(*nIt1, *nIt2, G);

            if(!ep.second){ untested();
				}else if(cb){ untested();
					(*cb)(*nIt1, *nIt2);
            }else{ untested();
				}

            if(!boost::is_directed(G)){ untested();
            }else if( boost::edge(*nIt2, *nIt1, G).second ){ untested();
                // change later.
            }else{ untested();
                boost::add_edge(*nIt2, *nIt1, G);
            }
        }
    }
}

} // impl

}

#endif // guard
