// Felix Salfelder, 2016
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//

/* parse .gr files. undirected graphs
 *
 * c A comment ...
 * c first noncommment:
 * c "p tw" followed by number of vertices, number of edges.
 * p tw 16 32
 * c A comment ...
 * 1 2
 * 1 14
 * c these are edges. counting starts at 1.
 * 1 8
 * 1 9
 * ...
 * c line number determined by header and comments
 * c #comments + #vertices + #edges + 1.
 * 16 15
 */

#include <istream>
#include <fstream>
#include <iostream>
#include "assert.h"

static const unsigned BUFSIZE=256;

inline void myatoi(const char*& line, unsigned& v, char end=0)
{
	v = (*line)-'0';
	++line;
	while(*line!=end){ itested();
		v *= 10;
		v += (*line)-'0';
		++line;
	}
}

inline void myas2is(const char* line, unsigned& v, unsigned& w)
#if 0 // fallback.
{
	sscanf(line, "%d %d", &v, &w);
}
#else
{ // faster
	myatoi(line, v, ' ');
	++line;
	myatoi(line, w);
}
#endif
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
typedef enum{
	oLOWER = -1,
	oNONE = 0,
	oUPPER =1
} pORIENT_t;
/*--------------------------------------------------------------------------*/
class PARSE;
namespace detail {
	template<class G, class C=typename boost::graph_traits<G>::vertex_iterator::iterator_category>
	struct parse_backend{ //
		static void do_it(G& g, typename boost::graph_traits<G>::vertex_iterator base, PARSE& P);
	};
} // detail
/*--------------------------------------------------------------------------*/
class PARSE{
public:
	PARSE(const std::string& filename, pORIENT_t orient=oNONE)
		: _is(new std::ifstream(filename)), _orient(orient)
	{
		init(); // hmm delegate constructor?
	}
	PARSE(std::istream& is, pORIENT_t orient=oNONE)
		: _is(&is), _orient(orient)
	{
		init(); // hmm delegate constructor?
	}
#if 0
	GRPARSE(int fd){
		FILE* f=fdopen(fd);
		if(!f) throw "cannot open fd";
		init(f);
	}
#endif

	size_t num_vertices() const { return _v;}
	size_t num_edges() const { return _e; }

	template<class G>
	void parse_boost(G&);
private:
	void init()
	{
		_line=(char*) malloc(BUFSIZE);
		get_line();
		parse_header();
	}
	char* get_line(){
		_is->getline(_line, BUFSIZE);
		skip_comments();
		if(_is->eof()){
			return NULL;
		}else{
			return _line;
		}
	}
	void skip_comments()
	{
		while(*_line=='c'){
			get_line();
		}
	}
	void parse_header()
	{
		if(*_line!='p'){ untested();
			throw "there's no header";
		}
		char* seek=_line+2;

		// BUG, should match "tw" (or a fixed string?!)
		while(*seek!=' '){
			++seek;
		}
		++seek;
		myas2is(seek, _v, _e);
		if(errno){ untested();
			throw "p line?\n";
		}else{
		}
	}

	private:
	std::istream* _is;
	char* _line;
	unsigned _v;
	unsigned _e;

public:
	class iterator{
	public:
		friend class PARSE; //??
	public:
		typedef std::pair<unsigned, unsigned> value_type;
		iterator(PARSE& p, bool end=false) : _p(p)
		{
			if(end){
				_value.first = -1;
			}else{
				_p.next(*this);
			}
		}
		bool operator==(const iterator& o) const
		{ itested();
			assert(_value.first==-1u || o._value.first==-1u);
			return(_value.first==-1u && o._value.first==-1u);
		}
		bool operator!=(const iterator& o) const
		{
			return !operator==(o);
		}
		iterator& operator++()
		{
			_p.next(*this);
			return *this;
		}
		value_type const& operator*() const
		{ itested();
			return _value;
		}
		value_type const* operator->() const
		{ itested();
			return &_value;
		}
	private:
		PARSE& _p;
		value_type _value;
	}; // iterator

	iterator begin()
	{
		return iterator(*this);
	}
	iterator end()
	{
		return iterator(*this, true);
	}
	template<class P>
	void orient(P& p)
	{
		switch(_orient){ untested();
			case oNONE:
				break;
			case oUPPER:
				if(p.first < p.second){
					std::swap(p.first, p.second);
				}else{ itested();
				}
				break;
			case oLOWER: untested();
				if(p.first > p.second){ untested();
					std::swap(p.first, p.second);
				}else{ untested();
				}
				break;
		}
	}

	void next(iterator& i)
	{ itested();
		if(char* line=get_line()){
			myas2is(line, i._value.first, i._value.second);
			assert(i._value.second);
			assert(i._value.first);
			orient(i._value);
			--(i._value.second);
			--(i._value.first);
		}else{
			i._value.first=-1;
		}
	}
	template<class G, class C>
	friend struct detail::parse_backend;
	pORIENT_t _orient;
}; // PARSE
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
template<class G>
void PARSE::parse_boost(G&g)
{
	// g.resize(_v)?
	assert(boost::num_vertices(g) == _v);
	assert(boost::num_edges(g) == 0);
	typename boost::graph_traits<G>::vertex_iterator base=boost::vertices(g).first;
	detail::parse_backend<G, typename boost::graph_traits<G>::vertex_iterator::iterator_category>::
		do_it(g, base, *this);
}
/*--------------------------------------------------------------------------*/
namespace detail{
	template<class G, class C>
	void parse_backend<G,C>::do_it(G& g,
		                typename boost::graph_traits<G>::vertex_iterator base, PARSE& P)
	{
		incomplete(); // slow.
//			--base; // .gr starts at 1
		while(const char* line = P.get_line()){ itested();
			unsigned /*?*/ v, w;
			myas2is(line, v, w);
			assert(v);
			assert(w);
			assert(v<=P._v);
			assert(w<=P._v);
			assert(!errno); //incomplete.
			typename boost::graph_traits<G>::vertex_descriptor V, W;
			V = *(base+(v-1));
			W = *(base+(w-1));
			assert(V!=W);
			boost::add_edge(V, W, g);
		}
	}

	template<class G>
	struct parse_backend<G, std::random_access_iterator_tag>{ //
		static void do_it(G& g, typename boost::graph_traits<G>::vertex_iterator base, PARSE& P)
		{
			--base; // .gr starts at 1
			size_t E = boost::num_edges(g);
			assert(!E);
			// use iterator?!
			while(char* line = P.get_line()){ itested();
				unsigned /*?*/ v, w;
				myas2is(line, v, w);
				assert(v);
				assert(w);
				assert(v<=P._v);
				assert(w<=P._v);
				assert(!errno); //incomplete.
				typename boost::graph_traits<G>::vertex_descriptor V, W;
				V = *(base+v);
				W = *(base+w);
				assert(V!=W);
				boost::add_edge(V, W, g);
				++E;
				assert(E==boost::num_edges(g));
			}
			(void) E;
		}
	}; // parse_backend

}
