/* (c) Felix Salfelder 2021
 * License: GPLv3+
 *
 * Derived from pace_meiji (c) 2017, Hisao Tamaki, MIT license
 */

#include "tr_bs.hpp"

template<class key_type, class value_t, unsigned MAX_CHILDREN_SIZE=512>
class LayeredSieve {
//  int _n;
  unsigned _tbs;

	typedef key_type cfg_myset;
  typedef BlockSieve<cfg_myset, value_t, MAX_CHILDREN_SIZE> bs_t;
  std::vector<bs_t> _sieves;

public:
  typedef bs_t::set_out_hack set_out_hack;
  
public:
   LayeredSieve(){
	}
	
//  explicit LayeredSieve(int n, int targetWidth)
	 // : _n(n)
 void init(int n, int targetWidth) {
//    this.targetWidth = targetWidth;
	 _tbs = targetWidth + 1;
    
    int k = 33 - std::countl_zero(_tbs - 1);
	 trace1("alloc sieves", k);
    _sieves.resize(0);
    _sieves.resize(k); //  = new BlockSieve[k];
    for (int i = 0; i < k; i++) {
//      int margin = (1 << i) - 1;
		// sieve[0] margin 0
		// sieve[1] margin 1
		// sieve[2] margin 3
		// sieve[3] margin 7
		 _sieves[i].set_n(n);
//      _sieves[i].settargetWidth(targetWidth);
//      _sieves[i].setMargin(margin);
    }
  }
  
public:
#if 1
  value_t& put(cfg_myset const& vertices, value_t const& neighbors) {
    int ns = BlockSieve<cfg_myset,value_t, MAX_CHILDREN_SIZE>::get_size(neighbors);
//	 trace3("put", vertices, ns, neighbors);
	 return put(vertices, ns, neighbors);
//     int margin = targetWidth + 1 - ns;
//     int i = 32 - Integer.numberOfLeadingZeros(margin);
//     sieves[i].put(vertices, neighbors);
  }
private:
#endif
  value_t& put(cfg_myset const& vertices, int neighborSize, value_t const& value) {
    int margin = _tbs - neighborSize;
//	 assert(_tbs>=neighborSize); //?
//    margin=0 -> i = 0
//    margin=1 -> i = 1
//    margin=2 -> i = 2
//    margin=3 -> i = 2
//    margin=4 -> i = 3
//    margin=5 -> i = 3
//    margin=6 -> i = 3
//    margin=7 -> i = 3
//    margin=8 -> i = 4
//
//    int i ~ ceil(log2(margin+1))?
    // int i = 32 - Integer.numberOfLeadingZeros(margin);
	 int i = 32 - __builtin_clz(margin);
	 assert(i>=0);
//	 trace3("put", margin, i, _sieves.size());
	 assert(i<int(_sieves.size()));
    return _sieves[i][vertices] = value;
  }

public:
	void collectSuperblocks(cfg_myset const& component, cfg_myset const& neighbors,
        set_out_hack& list) {

		int i = 0;
		int mm = 1;
		for (auto& sieve : _sieves) {
			int margin = (1 << i) - 1;
			assert(margin == mm-1);
			sieve.collectSuperblocks(_tbs, mm - 1, component, neighbors, list);
			++i;
			mm = mm << 1;
		}
	}
  
	std::ostream& stats(std::ostream& o) const {
		o << "[";
		std::string comma = "";
		for (size_t i=0; i < _sieves.size(); i++) {
			o << comma << _sieves[i].size();
			comma = ", ";
		}
		o << "]";
		return o;
	}
	std::vector<int> getSizes() {
		std::vector<int> sizes(_sieves.size());
    for (int i = 0; i < _sieves.size(); i++) {
      sizes[i] = _sieves[i].size();
    }
    return sizes;
  }
};

