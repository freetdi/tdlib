// Felix Salfelder, 2017, 2021
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

// python interface common types.

#ifndef GRAPH_PY_HPP
#define GRAPH_PY_HPP

#include "config.h"
#include <boost/graph/graph_traits.hpp>
#include "graph_traits.hpp"
#include "treedec_traits.hpp"
#include "boost_compat.h"

#ifdef HAVE_GALA_GRAPH_H
#include "gala_graph_py.hpp"
#endif

#include <boost/python.hpp>

#define setS boost::setS
#define vecS boost::vecS
#define directedS boost::directedS
#define undirectedS boost::undirectedS

typedef boost::adjacency_list<vecS, vecS, directedS>   _balvvd;
typedef boost::adjacency_list<setS, vecS, directedS>   _balsvd;
typedef boost::adjacency_list<vecS, vecS, undirectedS> _balvvu;
typedef boost::adjacency_list<setS, vecS, undirectedS> _balsvu;

typedef typename treedec::graph_traits<_balvvd>::treedec_type _balvvd_treedec;
typedef typename treedec::graph_traits<_balsvd>::treedec_type _balsvd_treedec;
typedef typename treedec::graph_traits<_balvvu>::treedec_type _balvvu_treedec;
typedef typename treedec::graph_traits<_balsvu>::treedec_type _balsvu_treedec;

// graphs used in python-tdlib
typedef boost::adjacency_list<setS, vecS, undirectedS> TD_graph_t; //type 0
typedef boost::adjacency_list<vecS, vecS, undirectedS> TD_graph_vec_t; //type 1
typedef boost::adjacency_list<setS, vecS, directedS> TD_graph_directed_t; //type 2
typedef boost::adjacency_list<vecS, vecS, directedS> TD_graph_directed_vec_t; //type 3

#undef setS
#undef vecS
#undef directedS
#undef undirectedS

#endif
