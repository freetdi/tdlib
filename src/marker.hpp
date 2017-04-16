#pragma once

// inspired by boost/graph/md/marker
// could just include if there was a seperate header...

namespace treedec{

namespace draft{

template<class U, class key_type /* idmap... */>
class sMARKER{
	BOOST_STATIC_ASSERT( std::numeric_limits<U>::is_integer
	                 && !std::numeric_limits<U>::is_signed);
	BOOST_STATIC_ASSERT( std::numeric_limits<key_type>::is_integer
	                 && !std::numeric_limits<key_type>::is_signed);
private: // hide
	sMARKER(){ unreachable(); }
	sMARKER(const sMARKER&){ unreachable(); }
public:
	sMARKER(size_t howmany)
		: _tide(1 /*so we can unmark*/),
		  _tags(howmany)
	{
		assert(sizeof(U));

	}
	~sMARKER(){
	}
	void clear(){
		if(_tide==std::numeric_limits<U>::max()){ untested();
			reset();
		}else{
			++_tide;
		}
	}
	void mark(key_type x){
		assert(x<_tags.size());
		_tags[x] = _tide;
	}
	void unmark(key_type x){
		_tags[x] = _tide-1;
	}
	bool is_marked(key_type x) const{
		return(_tags[x] == _tide);
	}
private:
	void reset(){ untested();
		std::fill(_tags.begin(), _tags.end(), 0);
		_tide = 1;
	}

private:
	U _tide;
	std::vector<U> _tags;
};

} // draft

// mark neighbours of v up to v.
template<class M, typename V, class G>
void mark_smaller_neighbours(M& marker, V v, G const& g)
{ untested();
    auto pp=adjacent_vertices(v, g);
    for(; pp.first!=pp.second; ++pp.first){
		 if(*pp.first>=v){ untested();
			 // break; // need sorting...
			 continue;
		 }else{ untested();
			 marker.mark(*pp.first);
		 }
    }
}

} // treedec
