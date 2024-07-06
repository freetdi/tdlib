/* (c) Felix Salfelder 2021
 * License: GPLv3+
 *
 * Derived from pace_meiji (c) 2017, Hisao Tamaki, MIT license
 */

#ifndef TR_BAG_H
#define TR_BAG_H

#include "tr_sep.h"

#include <ostream>


template<class G>
class Separator;

static const unsigned INVALID_NODE = -1;

template<class G>
class IODecomposer;

// is this (just) a graph?
template<class G>
class Bag{
public: // types
	typedef Separator<G> sep_t;
	typedef Bag<G> bag_t;

	friend class IODecomposer<G>;

private: // data
	Bag const* _parent{nullptr};
	myset _vertices; // vertices of the paren graph?  // map?
	G const* _graph{nullptr};
	unsigned _size{-1u};
	std::vector<int> _conv;
	std::vector<int> _inv;
	std::vector<Bag const*> _nestedBags;
	std::vector<sep_t*> _separators;
	std::vector<Separator<G> const* > _incidentSeparators;
  int _width{0}; //bs?
  int _bagsize{0}; //bs?
  int _separatorWidth;
  int _lowerBound;
  int _inheritedLowerBound;
  bool _optimal;
  // SafeSeparator _ss;
  //
public:
	int size() const{return _size;}

	explicit Bag(G const& g)
	  : _graph(&g),
	    _size(boost::num_vertices(g)) {
		auto vp = boost::vertices(g);
		unsigned k = 0;
		for(; vp.first!=vp.second; ++vp.first){
			auto v=*vp.first;
			assert(v==k);
			_vertices.insert(v);
			++k;
		}
  }
  G const& graph() const{assert(_graph); return *_graph; }
  myset const& vertices(){ return _vertices; }
  
	explicit Bag(Bag const* parent, myset const& vertices)
	  : _parent(parent),
	    _vertices(vertices) {
		_size = _vertices.size();
	}
  
  void initializeForDecomposition() {
    if (_graph){
	 }else{
		 assert(_parent && "graph missing?!");
		 makeLocalGraph();
    }
    assert(_nestedBags.empty());
    assert(_separators.empty());
	 assert (!_width);
    _width = 0;
    _separatorWidth = 0;
  }
  
	void attachSeparator(Separator<G> const& separator) {
		// refcnt??
		_incidentSeparators.add(&separator);
	}

  void makeRefinable() {
    makeLocalGraph();
    _nestedBags.clear();
    _separators.clear();
  }
  
	int maxNestedBagSize() const {
		if (_nestedBags.size()){
			int max = 0;
			for (auto bag : _nestedBags) {
				if (bag->size() > max) {
					max = bag->size();
				}
			}
			return max;
		}else{
			return -1;
		}
	}
  
  Bag const& addNestedBag(myset const& _vertices) {
    bag_t const* bag = new Bag(this, _vertices);
    _nestedBags.push_back(bag);
    return *bag;
  }

  Separator<G>* addSeparator(myset const& _vertices) {
    auto s = new Separator<G>(this, _vertices);
    _separators.push_back(s);
    return s;
  }
  
  void addIncidentSeparator(Separator<G> const* separator) {
    _incidentSeparators.push_back(separator);
  }

  // import edges from parent graph.
	void makeLocalGraph() {
		// when is size!=_graph.size?
		_graph = new G(size);
		_conv.resize(_parent.size(), INVALID_NODE);
		_inv.resize(size); // too big? vertices_size?
    
		int k = 0;
		for (int v : _vertices){
			_conv[v] = k;
			_inv[k++] = v; // _vertices array.
		}

//		_graph.inheritEdges(parent.graph, conv, inv);
		auto vv = boost::vertices(_graph);
		for (; vv.first!=vv.second; ++vv.first) {
			auto v = *vv.first;
			int x = _inv[v];
			auto bb = boost::adjacent_vertices(_parent->_graph, x);
			for (; bb.first!=bb.second; ++bb.first) {
				int u = _conv[*bb.first];
				assert (u != v);
				if (u > v) {
					boost::add_edge(_graph, u,  v);
				}else{
				}
			}
		}

		// import separators from parent?
    for (Separator separator: _incidentSeparators) {
//      System.out.println("filling " + separator);
			myset convd = convert(separator._vertices, _conv);
//			graph.fill(convd);
			for(auto i=convd.begin(); i!=convd.end();){
				auto j = i;
				++j;
				auto i_ = j;
				for(; j!=convd.end(); ++j){
					boost::add_edge(_graph, *i, *j);
				}
				i = i_;
			}
    }
  }
  
   int getWidth() {
    if (_nestedBags.empty()) {
      return size - 1;
    }else{
	 }
    int max = 0;
    for (Bag bag: _nestedBags) {
      int w = bag.getWidth();
      if (w > max) {
        max = w;
      }
    }
    for (Separator separator: _separators) {
      int w = separator.size(); // sic.
      if (w > max) {
        max = w;
      }
    }
    return max;

  }
  
	// set?!
	void computeBagSize() {
		// assumes that the bag is flat
		// assert(is_flat); // what does it mean??

		//    System.out.println("setWidth for " + this._vertices);
		//    System.out.println("nestedBags = " + nestedBags);
		_bagsize = 0;
		_separatorWidth = 0;

		if (_nestedBags.empty()){
			_bagsize = size;
			return;
		}else{
			for (Bag bag: _nestedBags) {
				if (bag.size > _bagsize) {
					_bagsize = bag.size;
				}else{
				}
			}

			for (Separator separator: _separators) {
				if (separator.size > _separatorWidth) {
					_separatorWidth = separator.size;
				}
			}

			if (_separatorWidth + 1 > _bagsize ) {
				_bagsize = _separatorWidth + 1;
			}
		}

		_width = _bagsize - 1;
	}
  
   void flatten() {
    if (_nestedBags.empty()){
      return;
    }else{
	 }
    
    validate();
    for (Bag bag: _nestedBags) {
		 bag.flatten();
    }
    validate();
	 std::vector<sep_t*> newSeparatorList;

    for (sep_t* separator: _separators) {
//      System.out.println(separator.incidentBags.size() + " incident bags of " + 
//          separator);
		 std::vector<bag_t*> newIncidentBags;
		 for (bag_t* bag: separator._incidentBags) {
			 if (bag.parent != this){
				 newIncidentBags.add(bag);
			 }else if (bag.nestedBags.empty()) {
				 newIncidentBags.add(bag);
			 }else{
				 myset cs = convert(separator._vertices, bag->_conv);
				 bag_t* nested = bag->findNestedBagContaining(cs);

				 if (nested) {
				 }else{
					 // error??
					 bag->dump();
					 std::cerr <<" does not have a bag containing " 
						 << cs << " which is originally " << separator._vertices << "\n";
					 dump();
				 }

				 newIncidentBags.add(nested);
				 nested.addIncidentSeparator(separator);
			 }
      }
      if (!newIncidentBags.isEmpty()) {
        separator.incidentBags = newIncidentBags;
        newSeparatorList.add(separator);
      }
//      System.out.println("processed separator :" + separator);
    }
    _separators = newSeparatorList;
    
	 std::vector<bag_t*> temp;
	 std::swap(temp, _nestedBags);
    for (bag_t* bag: temp) {
		 assert(bag);
		 if (bag->_nestedBags.empty()) {
			 //        System.out.println("adding original bag " + bag._vertices);
			 _nestedBags.push_back(bag);
		 }else{
			 for (bag_t* nested: bag._nestedBags) {
				 //          System.out.println("inverting " + nested);
				 nested->invert();
				 _nestedBags.add(nested);
				 //          System.out.println("inverted " + nested);
			 }
			 for (sep_t* s : bag->_separators) {
				 //          System.out.println("inverting sep " + separator);
				 s->invert();
				 _separators.push_back(s);
				 //          System.out.println("inverted sep " + separator);
			 }
		 }
	 }
    computeBagSize();
//    System.out.println("bag of size " + size + " flattened into " + nestedBags.size() + " bags and width " +
//        width);
//    for (Bag bag: nestedBags) {
//      System.out.println("incident separators of " + bag._vertices);
//      for (Separator s: bag.incidentSeparators) {
//        System.out.println("  " + s._vertices);
//        for (Bag b: s.incidentBags) {
//          System.out.println("        " + b._vertices);
//        }
//      }
//    }
  }
  
  bag_t const* findNestedBagContaining(myset _vertices) {
    for (bag_t* bag: _nestedBags) {
      if (_vertices.is_subset_of(bag._vertices)) {
        return bag;
      }else{
		}
    }
    return nullptr;
  }

   void invert() {
		assert(_parent);
		_vertices = convert(_vertices, _parent->_inv);
		_parent = _parent->_parent;
	}
  
	void convert() {
		_vertices = convert(_vertices, _parent->_conv);
	}
  
   myset convert(myset const& s) const {
    return convert(s, _conv);
  }
  
	template<class MAP>
  myset convert(myset const& s, MAP const& conv) {
    if (conv.size() < s.size()) {
		 assert(false);
    }
    myset result(conv.size());
    for (auto v : s){
      result.insert(conv[v]);
    }
    return result;
  }
  
	void* toTreeDecomposition() {
		incomplete();
		return nullptr;
#if 0
    computeBagSize();
    TreeDecomposition td = new TreeDecomposition(0, width, _graph);
    for (Bag bag: nestedBags) {
      td.addBag(bag._vertices.toArray());
    }
    
    for (Separator separator: separators) {
       myset const& vs = separator.vertices();
      Bag full = null;
      for (Bag bag: separator.incidentBags) {
        if (vs.isSubset(bag._vertices)) {
          full = bag;
          break;
        }
      }
 
      if (full != null) {
        int j = nestedBags.indexOf(full) + 1;
        for (Bag bag: separator.incidentBags) {

          if (bag != full) {
            td.addEdge(j, nestedBags.indexOf(bag) + 1);
          }
        }
      }
      else {
        int j = td.addBag(separator._vertices.toArray());
        for (Bag bag: separator.incidentBags) {
          td.addEdge(j, nestedBags.indexOf(bag) + 1);
        }
      }
    }
    
    return td;
#endif
  }
  
#if 0
  void detectSafeSeparators() {
    ss = new SafeSeparator(graph);
    for (Separator separator: separators) {
//      separator.figureOutSafetyBySPT();
      separator.figureOutSafety(ss);
    }
  }
#endif
  
  void pack() {
	  std::vector<Bag*> nestedBags_uc;
	  for (bag_t* bag: _nestedBags) {
		  if (bag->_parent == this) {
			  std::vector<bag_t*> bagsToPack;
			  bag.collectBagsToPack(bagsToPack, nullptr);
			  //        System.out.println("bags to pack: " + bagsToPack);
			  if (bagsToPack.size() >= 2) {
				  myset vertices(boost::num_vertices(*_graph));
				  for (bag_t* toPack: bagsToPack) {
					  vertices.merge(toPack->_vertices);
				  }
				  Bag packed = new Bag(this, _vertices);
				  packed.initializeForDecomposition();
				  packed.nestedBags = bagsToPack;
				  for (Bag toPack: bagsToPack) {
					  toPack.parent = packed;
					  toPack.convert();
				  }
				  nestedBags_uc.push_back(packed);
			  }
			  else {
				  nestedBags_uc.push_back(bag);
			  }
		  }
	  }
	  _nestedBags = nestedBags_uc;

	  // ArrayList<Separator> newSeparatorList = new ArrayList<>();
	  std::vector<sep_t*> seplist_uc;

	  for (Separator separator: _separators) {
		  bool internal = true;
		  bag_t* parent = nullptr;
		  for (Bag b: separator._incidentBags) {
			  if (b._parent == this) {
				  internal = false;
				  break;
			  } else if (!parent) {
				  parent = b._parent;
			  } else if (b._parent != parent) {
				  internal = false;
				  break;
			  }
		  }
		  if (internal) {
			  separator.parent = parent;
			  separator.convert();
			  parent.separators.add(separator);
		  }else{
			  std::vector<bag_t*> newIncidentBags;
			  for (Bag b: separator._incidentBags) {
				  if (b._parent == this) {
					  newIncidentBags.push_back(b);
				  } else {
					  newIncidentBags.push_back(b._parent);
					  b._parent._incidentSeparators.push_back(separator);
					  b._incidentSeparators.erase(separator);
				  }
			  }
			  separator.incidentBags = newIncidentBags;
			  seplist_uc.push_back(separator);    
		  }
	  }

	  _separators = seplist_uc;

	  for (bag_t* bag: _nestedBags) {
		  bag->computeBagSize();
	  }
	  computeBagSize();
  }
  
  // Bag::
  void collectBagsToPack(std::vector<bag_t*> list, sep_t* from) {
    list.push_back(this);
    for (Separator separator: _incidentSeparators) {
//      System.out.println(" safe = " + separator.safe);
      if (separator == from){
		}else if( separator.safe){
		}else if( separator.wall){
      }else{
			separator.collectBagsToPack(list,  this);
		}
    }
  }
  
  int countSafeSeparators() {
    int count = 0;
    for (Separator separator: _separators) {
      if (separator.is_safe()) {
        ++count;
      }else{
		}
    }
    return count;
  }
  
  void validate() {
    if (!_nestedBags.empty()){
//      assert !nestedBags.isEmpty() : "no nested bags " + this; 
      for (bag_t* b: _nestedBags) {
        b->validate();
        assert(!b._vertices.empty() && "empty bag");
        assert(b._parent == this && "parent");
//            "\n which is " + b.parent +
//            "\n is supposed to be " + this;
      }
      for (Separator s: _separators) {
        assert(!s._vertices.isEmpty() && "empty seprator ");
        assert(s._parent == this);
//            "\n which is " + s.parent +
//            "\n is supposed to be " + this;
      }
      for (bag_t* b: _nestedBags) {
        for (Separator s: b._incidentSeparators) {
          assert(!s._vertices.isEmpty());
          assert(s._parent == this);
			// : "parent of " + s + 
         //     "\n which is " + s.parent +
         //     "\n is supposed to be " + this + 
         //     "\n where the separator is incident to bag " + b;
          assert(s._vertices.is_subset_of(b._vertices));
			//  : "separator vertex set " + s._vertices + 
         //  "\n is not a subset of the bag vertex set " + b._vertices;
        }
      }
      for (sep_t* separator: _separators) {
        for (bag_t* b : separator._incidentBags) {
          assert(b);
          assert(b._parent == this);
			// : "parent of " + b + 
         //     "\n which is " + b.parent +
         //     "\n is supposed to be " + this + 
         //     "\n where the bag is incident to separator " + separator;
          assert(separator._vertices.is_subset_of(b._vertices));
			// : "separator vertex set " + 
         //     separator._vertices + 
         // "\n is not a subset of the bag vertex set " + b._vertices;
        }
      }
    }
  }
  
	void dump(std::ostream& o, std::string indent = "") {
		o << indent << "bag:" << _vertices
		  << indent << "width = " << _width << ", conv = "
		  //<< _conv
		  ;
    //if (_nestedBags != null) {
	 //}
	 {
    //  System.out.println(indent + nestedBags.size() + " subbags:"); 
      for (bag_t* bag: _nestedBags) {
        bag->dump(indent + "  ");
      }
      for (auto separator: _separators) {
        separator->dump(indent + "  ");
      }
    }
  }
  
	void canonicalize() {
		bool moving = true;
		while (moving = true) {
			moving = false;
			for (auto bag: _nestedBags) {
				if (bag.trySplit()) {
					moving = true;
				}else{
				}
			}
			if (moving) {
				flatten();
			}else{
			}
		}
	}

	bool trySplit() {
		return false;
	}

	std::ostream& print(std::ostream& o) const {
		if (_parent){
			o << "bag " << _parent->indexOf(this) << ": ";
		} else {
			o << "root bag : ";
		}
		o << _vertices;
		return o;
	}

  int compare(Bag const& b) const {
    if (size() != b.size()) {
      return b.size() - size();
    }
    return b._vertices.compare_int(_vertices);
  }
};

#endif
