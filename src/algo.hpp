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

#ifndef TD_ALGO_HPP
#define TD_ALGO_HPP

#include <limits>
#include "timer.hpp"

namespace treedec{
namespace algo{

template<class GraphType>
struct default_config{
	typedef typename boost::graph_traits<GraphType>::vertices_size_type vst;
	static constexpr unsigned max_vertex_index=std::numeric_limits<vst>::max();
};

namespace draft{

template<class O>
O* clone(O const* o)
{
    if(o){
        return o->clone();
    }
    else{
        return NULL;
    }
}

class algo1{
protected:
	algo1(const algo1& o)
		: _label(o._label), _timer(clone(o._timer)){}
public:
	algo1(std::string label)
	   : _label(label), _timer(NULL) {
#ifdef TIMER
			_timer=new DOUBLE_TIMER("raw");
#endif
		}
	virtual ~algo1(){
		if(_timer){
			std::cout << _label << ": ";
			std::cout << *_timer << "\n";
		}
	}

	virtual void do_it() = 0;

//for now
public:
	double get_runtime(){
            if(_timer){
                return _timer->total();
            }
            return -1;
        }

protected:
	void timer_on(){
		if(_timer){
			_timer->start();
		}
	}
	void timer_off(){
		if(_timer){
			_timer->stop();
		}
	}

private:
	std::string _label;
	TIMER_BASE* _timer;
};

} // draft

} // algo1

} // treedec

#endif // guard
