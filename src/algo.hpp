/* Copyright (C) 2017 Felix Salfelder
 * Author: Felix Salfelder
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * graph decomposition algorithms, traits.
 */

#ifndef TREEDEC_ALGO_HPP
#define TREEDEC_ALGO_HPP

#include <limits>
#include "timer.hpp"
#include "graph_traits.hpp"

#define TREEDEC_ALGO_TC template <typename GGG, template<class GG_, class...> class CFGG> class

namespace treedec{

namespace algo{

struct config_base{
    static void interruption_point(){
        // compile-time disabled signalling
    }
    static void commit_lb(unsigned, std::string=""){ untested();
    }
    static void commit_ub(unsigned, std::string=""){ untested();
    }
    // "C" style seems to be most practical.
    static void message(unsigned badness, const char* fmt, ...){
        // compile-time disabled messaging
        (void) badness; (void) fmt;
    }
};


template<class GraphType, class ... rest>
struct default_config : config_base {
    using vst=typename boost::graph_traits<GraphType>::vertices_size_type;
    static constexpr unsigned max_vertex_index=std::numeric_limits<vst>::max();
};

namespace draft{

template<class O>
O* clone(O const* o)
{
    if(o){ untested();
        return o->clone();
    }else{
        return NULL;
    }
}

class algo1{
protected:
	algo1(const algo1& o)
#ifdef TIMER
		: _timer(clone(o._timer))
#endif
		{(void)o;}
public:
	algo1(std::string label)
#ifdef TIMER
	   : _timer(NULL) {
			_timer=new DOUBLE_TIMER("raw" + label);
		}
#else
	{(void)label;}
#endif
	virtual ~algo1(){
#ifdef TIMER
		if(_timer){
		//	std::cout << _label << ": ";
			std::cout << *_timer << "\n";
		}
#endif
	}

	virtual void do_it() = 0;

public: // results. to be reimplemnted.
    unsigned lower_bound_bagsize() const{
        incomplete();
        return 0;
    }
    unsigned bagsize() const{
        incomplete();
        return 0;
    }
    template<class O>
    void get_elimination_ordering(O&) const{
        unreachable();
    }
//for now
public:
	double get_runtime(){
#ifdef TIMER
            if(_timer){
                return _timer->total();
            }
#endif
            return -1;
        }

protected:
	void timer_on(){
#ifdef TIMER
		if(_timer){
			_timer->start();
		}
#endif
	}
	void timer_off(){
#ifdef TIMER
		if(_timer){
			_timer->stop();
		}
#endif
	}

private:
//	std::string _label; incomplete.
#ifdef TIMER
	TIMER_BASE* _timer;
#endif
};

} // draft

} // algo1

} // treedec

#endif // guard

// vim:ts=8:sw=4:et
