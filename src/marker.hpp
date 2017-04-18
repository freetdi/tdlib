// Felix Salfelder, 2016-2017
//
// (c) 2014-2016 Goethe-Universität Frankfurt
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


} // treedec
