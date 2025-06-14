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
	
	typedef G::edge_container_type cfg_myset;
	typedef Separator<G> sep_t;
	typedef Bag<G> bag_t;
	
	
	
	friend class IODecomposer<G>;

private: // data
	Bag const* _parent{nullptr};
	cfg_myset _vertices; // vertices of the paren graph?  // map?
	G const* _graph{nullptr};
	unsigned _size{-1u};
	std::vector<int> _conv;
	std::vector<int> _inv;
	std::vector<Bag const*> _nestedBags;
	std::vector<sep_t*> _separators;
	mutable std::vector<Separator<G> const* > _incidentSeparators;
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
  cfg_myset const& vertices(){ return _vertices; }
  
	explicit Bag(Bag const* parent, cfg_myset const& vertices)
	  : _parent(parent),
	    _vertices(vertices) {
		_size = _vertices.size();
		_graph = &parent->graph();
	}
  
  void initializeForDecomposition() { untested();
    if (_graph){ untested();
	 }else{ untested();
		 assert(_parent && "graph missing?!");
		 makeLocalGraph();
    }
    assert(_nestedBags.empty());
    assert(_separators.empty());
	 assert (!_width);
    _width = 0;
    _separatorWidth = 0;
  }
  
	void attachSeparator(Separator<G> const& separator) { untested();
		// refcnt??
		_incidentSeparators.add(&separator);
	}

  void makeRefinable() { untested();
    makeLocalGraph();
    _nestedBags.clear();
    _separators.clear();
  }
  
	int maxNestedBagSize() const { untested();
		if (_nestedBags.size()){ untested();
			int max = 0;
			for (auto bag : _nestedBags) { untested();
				if (bag->size() > max) { untested();
					max = bag->size();
				}
			}
			return max;
		}else{ untested();
			return -1;
		}
	}
  
  Bag& addNestedBag(cfg_myset const& _vertices) {
    bag_t* bag = new Bag(this, _vertices);
    _nestedBags.push_back(bag);
    return *bag;
  }

  Separator<G>* addSeparator(cfg_myset const& _vertices) {
    auto s = new Separator<G>(this, _vertices);
    _separators.push_back(s);
    return s;
  }
  
  void addIncidentSeparator(const Separator<G> * separator) const {
    _incidentSeparators.push_back(separator);
  }

  // import edges from parent graph.
	void makeLocalGraph() { untested();
		// when is size!=_graph.size?
		_graph = new G(size);
		_conv.resize(_parent.size(), INVALID_NODE);
		_inv.resize(size); // too big? vertices_size?
    
		int k = 0;
		for (int v : _vertices){ untested();
			_conv[v] = k;
			_inv[k++] = v; // _vertices array.
		}

//		_graph.inheritEdges(parent.graph, conv, inv);
		auto vv = boost::vertices(_graph);
		for (; vv.first!=vv.second; ++vv.first) { untested();
			auto v = *vv.first;
			int x = _inv[v];
			auto bb = boost::adjacent_vertices(_parent->_graph, x);
			for (; bb.first!=bb.second; ++bb.first) { untested();
				int u = _conv[*bb.first];
				assert (u != v);
				if (u > v) { untested();
					boost::add_edge(_graph, u,  v);
				}else{ untested();
				}
			}
		}

		// import separators from parent?
    for (Separator separator: _incidentSeparators) { untested();
//      System.out.println("filling " + separator);
			cfg_myset convd = convert(separator._vertices, _conv);
//			graph.fill(convd);
			for(auto i=convd.begin(); i!=convd.end();){ untested();
				auto j = i;
				++j;
				auto i_ = j;
				for(; j!=convd.end(); ++j){ untested();
					boost::add_edge(_graph, *i, *j);
				}
				i = i_;
			}
    }
  }
  
   int getWidth() { untested();
    if (_nestedBags.empty()) { untested();
      return size - 1;
    }else{ untested();
	 }
    int max = 0;
    for (Bag bag: _nestedBags) { untested();
      int w = bag.getWidth();
      if (w > max) { untested();
        max = w;
      }
    }
    for (Separator separator: _separators) { untested();
      int w = separator.size(); // sic.
      if (w > max) { untested();
        max = w;
      }
    }
    return max;

  }
  
	// set?!
	void computeBagSize() { untested();
		// assumes that the bag is flat
		// assert(is_flat); // what does it mean??

		//    System.out.println("setWidth for " + this._vertices);
		//    System.out.println("nestedBags = " + nestedBags);
		_bagsize = 0;
		_separatorWidth = 0;

		if (_nestedBags.empty()){ untested();
			_bagsize = size;
			return;
		}else{ untested();
			for (Bag bag: _nestedBags) { untested();
				if (bag.size > _bagsize) { untested();
					_bagsize = bag.size;
				}else{ untested();
				}
			}

			for (Separator separator: _separators) { untested();
				if (separator.size > _separatorWidth) { untested();
					_separatorWidth = separator.size;
				}
			}

			if (_separatorWidth + 1 > _bagsize ) { untested();
				_bagsize = _separatorWidth + 1;
			}
		}

		_width = _bagsize - 1;
	}
  
   void flatten() { untested();
    if (_nestedBags.empty()){ untested();
      return;
    }else{ untested();
	 }
    
    validate();
    for (Bag bag: _nestedBags) { untested();
		 bag.flatten();
    }
    validate();
	 std::vector<sep_t*> newSeparatorList;

    for (sep_t* separator: _separators) { untested();
//      System.out.println(separator.incidentBags.size() + " incident bags of " + 
//          separator);
		 std::vector<bag_t*> newIncidentBags;
		 for (bag_t* bag: separator._incidentBags) { untested();
			 if (bag.parent != this){ untested();
				 newIncidentBags.add(bag);
			 }else if (bag.nestedBags.empty()) { untested();
				 newIncidentBags.add(bag);
			 }else{ untested();
				 cfg_myset cs = convert(separator._vertices, bag->_conv);
				 bag_t* nested = bag->findNestedBagContaining(cs);

				 if (nested) { untested();
				 }else{ untested();
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
      if (!newIncidentBags.isEmpty()) { untested();
        separator.incidentBags = newIncidentBags;
        newSeparatorList.add(separator);
      }
//      System.out.println("processed separator :" + separator);
    }
    _separators = newSeparatorList;
    
	 std::vector<bag_t*> temp;
	 std::swap(temp, _nestedBags);
    for (bag_t* bag: temp) { untested();
		 assert(bag);
		 if (bag->_nestedBags.empty()) { untested();
			 //        System.out.println("adding original bag " + bag._vertices);
			 _nestedBags.push_back(bag);
		 }else{ untested();
			 for (bag_t* nested: bag._nestedBags) { untested();
				 //          System.out.println("inverting " + nested);
				 nested->invert();
				 _nestedBags.add(nested);
				 //          System.out.println("inverted " + nested);
			 }
			 for (sep_t* s : bag->_separators) { untested();
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
//    for (Bag bag: nestedBags) { untested();
//      System.out.println("incident separators of " + bag._vertices);
//      for (Separator s: bag.incidentSeparators) { untested();
//        System.out.println("  " + s._vertices);
//        for (Bag b: s.incidentBags) { untested();
//          System.out.println("        " + b._vertices);
//        }
//      }
//    }
  }
  
  bag_t const* findNestedBagContaining(cfg_myset _vertices) { untested();
    for (bag_t* bag: _nestedBags) { untested();
      if (_vertices.is_subset_of(bag._vertices)) { untested();
        return bag;
      }else{ untested();
		}
    }
    return nullptr;
  }

   void invert() { untested();
		assert(_parent);
		_vertices = convert(_vertices, _parent->_inv);
		_parent = _parent->_parent;
	}
  
	void convert() { untested();
		_vertices = convert(_vertices, _parent->_conv);
	}
  
   cfg_myset convert(cfg_myset const& s) const { untested();
    return convert(s, _conv);
  }
  
	template<class MAP>
  cfg_myset convert(cfg_myset const& s, MAP const& conv) { untested();
    if (conv.size() < s.size()) { untested();
		 assert(false);
    }
    cfg_myset result(conv.size());
    for (auto v : s){ untested();
      result.insert(conv[v]);
    }
    return result;
  }
  
	void* toTreeDecomposition() { untested();
		incomplete();
		return nullptr;
#if 0
    computeBagSize();
    TreeDecomposition td = new TreeDecomposition(0, width, _graph);
    for (Bag bag: nestedBags) { untested();
      td.addBag(bag._vertices.toArray());
    }
    
    for (Separator separator: separators) { untested();
       cfg_myset const& vs = separator.vertices();
      Bag full = null;
      for (Bag bag: separator.incidentBags) { untested();
        if (vs.isSubset(bag._vertices)) { untested();
          full = bag;
          break;
        }
      }
 
      if (full != null) { untested();
        int j = nestedBags.indexOf(full) + 1;
        for (Bag bag: separator.incidentBags) { untested();

          if (bag != full) { untested();
            td.addEdge(j, nestedBags.indexOf(bag) + 1);
          }
        }
      }
      else { untested();
        int j = td.addBag(separator._vertices.toArray());
        for (Bag bag: separator.incidentBags) { untested();
          td.addEdge(j, nestedBags.indexOf(bag) + 1);
        }
      }
    }
    
    return td;
#endif
  }
  
#if 0
  void detectSafeSeparators() { untested();
    ss = new SafeSeparator(graph);
    for (Separator separator: separators) { untested();
//      separator.figureOutSafetyBySPT();
      separator.figureOutSafety(ss);
    }
  }
#endif
  
  void pack() { untested();
	  std::vector<Bag*> nestedBags_uc;
	  for (bag_t* bag: _nestedBags) { untested();
		  if (bag->_parent == this) { untested();
			  std::vector<bag_t*> bagsToPack;
			  bag.collectBagsToPack(bagsToPack, nullptr);
			  //        System.out.println("bags to pack: " + bagsToPack);
			  if (bagsToPack.size() >= 2) { untested();
				  cfg_myset vertices(boost::num_vertices(*_graph));
				  for (bag_t* toPack: bagsToPack) { untested();
					  vertices.merge(toPack->_vertices);
				  }
				  Bag packed = new Bag(this, _vertices);
				  packed.initializeForDecomposition();
				  packed.nestedBags = bagsToPack;
				  for (Bag toPack: bagsToPack) { untested();
					  toPack.parent = packed;
					  toPack.convert();
				  }
				  nestedBags_uc.push_back(packed);
			  }
			  else { untested();
				  nestedBags_uc.push_back(bag);
			  }
		  }
	  }
	  _nestedBags = nestedBags_uc;

	  // ArrayList<Separator> newSeparatorList = new ArrayList<>();
	  std::vector<sep_t*> seplist_uc;

	  for (Separator separator: _separators) { untested();
		  bool internal = true;
		  bag_t* parent = nullptr;
		  for (Bag b: separator._incidentBags) { untested();
			  if (b._parent == this) { untested();
				  internal = false;
				  break;
			  } else if (!parent) { untested();
				  parent = b._parent;
			  } else if (b._parent != parent) { untested();
				  internal = false;
				  break;
			  }
		  }
		  if (internal) { untested();
			  separator.parent = parent;
			  separator.convert();
			  parent.separators.add(separator);
		  }else{ untested();
			  std::vector<bag_t*> newIncidentBags;
			  for (Bag b: separator._incidentBags) { untested();
				  if (b._parent == this) { untested();
					  newIncidentBags.push_back(b);
				  } else { untested();
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

	  for (bag_t* bag: _nestedBags) { untested();
		  bag->computeBagSize();
	  }
	  computeBagSize();
  }
  
  // Bag::
  void collectBagsToPack(std::vector<bag_t*> list, sep_t* from) { untested();
    list.push_back(this);
    for (Separator separator: _incidentSeparators) { untested();
//      System.out.println(" safe = " + separator.safe);
      if (separator == from){ untested();
		}else if( separator.safe){ untested();
		}else if( separator.wall){ untested();
      }else{ untested();
			separator.collectBagsToPack(list,  this);
		}
    }
  }
  
  int countSafeSeparators() { untested();
    int count = 0;
    for (Separator separator: _separators) { untested();
      if (separator.is_safe()) { untested();
        ++count;
      }else{ untested();
		}
    }
    return count;
  }
  
  void validate() { untested();
    if (!_nestedBags.empty()){ untested();
//      assert !nestedBags.isEmpty() : "no nested bags " + this; 
      for (bag_t* b: _nestedBags) { untested();
        b->validate();
        assert(!b._vertices.empty() && "empty bag");
        assert(b._parent == this && "parent");
//            "\n which is " + b.parent +
//            "\n is supposed to be " + this;
      }
      for (Separator s: _separators) { untested();
        assert(!s._vertices.isEmpty() && "empty seprator ");
        assert(s._parent == this);
//            "\n which is " + s.parent +
//            "\n is supposed to be " + this;
      }
      for (bag_t* b: _nestedBags) { untested();
        for (Separator s: b._incidentSeparators) { untested();
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
      for (sep_t* separator: _separators) { untested();
        for (bag_t* b : separator._incidentBags) { untested();
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
  
	void dump(std::ostream& o, std::string indent = "") { untested();
		o << indent << "bag:" << _vertices
		  << indent << "width = " << _width << ", conv = "
		  //<< _conv
		  ;
    //if (_nestedBags != null) { untested();
	 //}
	 { untested();
    //  System.out.println(indent + nestedBags.size() + " subbags:"); 
      for (bag_t* bag: _nestedBags) { untested();
        bag->dump(indent + "  ");
      }
      for (auto separator: _separators) { untested();
        separator->dump(indent + "  ");
      }
    }
  }
  
	void canonicalize() { untested();
		bool moving = true;
		while (moving) { untested();
			moving = false;
			for (auto bag: _nestedBags) { untested();
				if (bag.trySplit()) { untested();
					moving = true;
				}else{ untested();
				}
			}
			if (moving) { untested();
				flatten();
			}else{ untested();
			}
		}
	}

	bool trySplit() { untested();
		return false;
	}

	std::ostream& print(std::ostream& o) const { untested();
		if (_parent){ untested();
			o << "bag " << _parent->indexOf(this) << ": ";
		} else { untested();
			o << "root bag : ";
		}
		o << _vertices;
		return o;
	}

  int compare(Bag const& b) const { untested();
    if (size() != b.size()) { untested();
      return b.size() - size();
    }
    return b._vertices.compare_int(_vertices);
  }
};

#endif

