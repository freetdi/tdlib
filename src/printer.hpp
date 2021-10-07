// Felix Salfelder, 2016-2017
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; version 3, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//

#ifndef TREEDEC_PRINTER_HPP
#define TREEDEC_PRINTER_HPP

#include <ostream>
#include <string>
#include "treedec.hpp"
#include "container.hpp"

namespace treedec{

namespace draft {

	// print a tree decomposition in td format.
	template<class G>
	class printer {
	public:
		typedef void vertex_bundled;
		typedef void vertex_property_type;
		typedef void edge_property_type;
	public:
		explicit printer(size_t nv) : _num_vertices(nv) {}
		explicit printer(std::ostream& s, G const& g, std::string reason="td")
		  : _num_vertices(0),
		    _nva(0),
		    // _bagsize??,
		    _offset(1),
		    _s(&s),
		    _reason(reason) {
			_ngv = boost::num_vertices(g);
		}
		~printer() {
			if(_s){
			  *_s << "\n";
			}else{
			}
		}

	public:
		printer& operator=(printer const& o){
			trace1("printer assign", o._num_vertices);
			_num_vertices = o._num_vertices;
			return *this;
		}	

		virtual void head(size_t numbags=0, size_t bagsize=0, size_t numvert=0) {
			(void)numvert;
			if(numbags==0 && bagsize==0){
			}else{
				// obsolete?
				_num_vertices = numbags;
				_bagsize = bagsize;
			}
			s() << "s " << _reason << " " << _num_vertices
			          << " " << _bagsize << " " << _ngv;
		}
		void set_bagsize(size_t b) {
			_bagsize = b;
			head();
		}
		size_t add_vertex() {
			return _nva;
		}
		void edge(size_t x, size_t y) { itested();
			s() << "\n" << x+_offset << " " << y+_offset;
		}
		void announce_bag(size_t x) {
			if(_nva==x){
			}else{ untested();
				incomplete(); //?
			}
			s() << "\nb " << ++_nva;
		}
		void push_back(size_t x) { untested();
			s() << " " << x+_offset;
		}
		// kludge: just tell something has been added.
		// (which is not completely wrong)
		std::pair<bool, bool> insert(size_t x) {
			s() << " " << x+_offset;
			return std::make_pair(true, true);
		}

		printer& bag(unsigned i){
			announce_bag(i);
			return *this;
		}
	private:
		std::ostream& s(){assert(_s); return *_s;}

	public:
		size_t num_vertices() const{
			return _num_vertices;
		}

	private:
		size_t _num_vertices{0};
		unsigned _nva{0};
		size_t _bagsize{0};
		unsigned _offset{0};
		std::ostream* _s{nullptr};
		size_t _ngv{0};
//		G const& _g; // needed? maybe vertex name map?
		std::string _reason;
	}; // printer

// fill bags. one at a time.
template<class G>
printer<G>& bag(size_t i, printer<G>& g) {
	return g.bag(i);
}

} // treedec::draft

template<class T>
struct treedec_traits<treedec::draft::printer<T> >{
	// typedef size_t vertex_descriptor;
	// typedef size_t edge_descriptor;
};

template <class T>
void set_bagsize(
		treedec::draft::printer<T>& t, size_t s)
{
	t.set_bagsize(s);
}

template<class T>
using grtdprinter=draft::printer<T>;

namespace detail{

#if 0
template<class G>
struct container_modify<treedec::draft::printer<G> >{

	// push, insert new item
	template<class E>
	static void push(C& c, E e) { untested();
		// incomplete(); // not supported by all containers.
		// assert(!container_inspect<C>::contains(c, e));
		c.insert(e);
	}
	template<class E>
	static void insert(C& c, E e) { untested();
		c.insert(e);
	}
	// sort. no-op for some containers
	static void insert(C& c) { untested();
		incomplete();
	}

	template<class C>
	static void sort(C&){
		// no-op.
	}
};
#endif

} // detail
} // treedec

namespace boost{

	template<class G>
	treedec::grtdprinter<G>& get(treedec::bag_t,
			treedec::grtdprinter<G>& g, unsigned i) {
		return g.bag(i);
	}

	// pass vertex_descriptors. for now
	template<class G>
	size_t add_vertex( treedec::grtdprinter<G>& p)
	{
		return p.add_vertex();
	}

	template<class G>
	size_t num_vertices( treedec::grtdprinter<G> const& p)
	{
		return p.num_vertices();
	}

	template<class G>
	std::pair<std::pair<size_t, size_t>, bool>
	add_edge(size_t x, size_t y, treedec::grtdprinter<G>& p)
	{
		p.edge(x,y);
		return std::make_pair( std::make_pair(x, y), true );
	}

	template<class T>
	struct graph_traits<treedec::grtdprinter<T> >{
		typedef size_t vertex_descriptor;
		typedef std::pair<size_t, size_t> edge_descriptor;
	};


	template<class G>
	struct property_map<treedec::grtdprinter<G>, vertex_all_t> {
		typedef treedec::grtdprinter<G> graph_type;
		typedef property_map<treedec::grtdprinter<G>, vertex_all_t> type;
		size_t& operator[](size_t& n) const { return n; }
	public:
		property_map( graph_type & g) : _g(g) {
		}

	public: // for now
		treedec::grtdprinter<G>& _g;
	};

	template<class G>
	struct property_map<treedec::grtdprinter<G>, vertex_all_t, treedec::bag_t> {
		typedef property_map<treedec::grtdprinter<G>, vertex_all_t, treedec::bag_t> type;
	public:
		property_map( treedec::grtdprinter<G>& g) : _g(g) { untested();
		}

	private:
		treedec::grtdprinter<G>& _g;
	};
	template<class G>
	struct property_map<treedec::grtdprinter<G>, edge_all_t> {
		typedef property_map<treedec::grtdprinter<G>, edge_all_t> type;
		size_t& operator[](std::pair<size_t, size_t>& n) const { return n.first; }
	};

	template<class G>
	property_map<treedec::grtdprinter<G>, vertex_all_t>
	get(vertex_all_t, treedec::grtdprinter<G>& g) {
		return property_map<treedec::grtdprinter<G>, vertex_all_t>(g);
	}
	template<class G>
	property_map<treedec::grtdprinter<G>, edge_all_t>
	get(edge_all_t, treedec::grtdprinter<G>&) {
		return property_map<treedec::grtdprinter<G>, edge_all_t>();
	}

	template<class G>
	void put(property_map<treedec::grtdprinter<G>, edge_all_t> const&,
	         std::pair<size_t, size_t>&, const no_property)
	{
	}

	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t, void> const&,
	         size_t, const no_property)
	{ untested();
	}
	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t, void> const&,
	         size_t, const unsigned)
	{ untested();
	}
	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t>
	         const&, size_t&)
	{ incomplete();
	}
	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t>
	         const&, size_t&, const treedec::bag_t)
	{ untested();
		unreachable();
	}
	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t,
	         treedec::bag_t> const&, size_t&, const treedec::bag_t)
	{ incomplete();
	}
	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t> const&,
			std::set<unsigned>&, const treedec::bag_t)
	{ incomplete();
	}

	template<class G>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t> const&,
			std::vector<unsigned>&, const treedec::bag_t)
	{ incomplete();
	}

	template<class G, class U>
	void put(boost::property_map<treedec::grtdprinter<G>, vertex_all_t>& m,
		 size_t v, property<treedec::bag_t, std::vector<U> > const& p)
	{
		auto& b=bag(v, m._g);
		for(auto i : p.m_value){
			treedec::push(b, i);
		}
	}
	template<class G, class U>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t, void>& m,
	         size_t v, const property<treedec::bag_t, std::set<U> >& p)
	{
		auto& b=bag(v, m._g);
		for(auto i : p.m_value){
			treedec::push(b, i);
		}
	}
	template<class G, class U>
	void put(property_map<treedec::grtdprinter<G>, vertex_all_t, void>& m,
	         size_t v, std::set<U> const& p)
	{ untested();
		auto& b=bag(v, m._g);
		for(auto i : p){ untested();
			treedec::push(b, i);
		}
	}

} // boost

#endif
