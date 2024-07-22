/* (c) Felix Salfelder 2021
 * License: GPLv3+
 *
 * Derived from pace_meiji (c) 2017, Hisao Tamaki, MIT license
 */

#ifndef TR_SEP_H
#define TR_SEP_H

#include "tr_bag.h"
#include "tr_ssep.h"

#include "bits/bool.hpp"
#include "graph_util.hpp"
#include "iter.hpp"

template<class G>
class Separator { // public Bag?
  typedef G::edge_container_type cfg_myset;
  typedef Bag<G> bag_t;
  private:
  Bag<G> const* _parent{nullptr};
  G const* _graph;
  cfg_myset _vertices;
  size_t _size;
  std::vector<bag_t const*> _incidentBags;
  public:
  bool safe{false};
  bool unsafe{false};
  bool wall{false};
  
  std::vector<unsigned> parentVertex;

  explicit Separator(bag_t const* parent)
	  :_parent(parent),
    _graph ( parent->_graph) {
   // incidentBags = new ArrayList<>();
  }

 explicit Separator(bag_t const* parent, cfg_myset const& vertexSet) 
	  :_parent(parent),
    _graph (&parent->graph()) {
    _vertices = vertexSet;
    _size = vertexSet.size();
  }
  
   void addIncidentBag(bag_t const* bag) {
    _incidentBags.push_back(bag);
  }
  
   void removeVertex(int v) {
    if (_vertices.contains(v)) {
       --_size;
		 _vertices.erase(v);
    }else{
	 }
  }
  
	void invert() {
		assert(_parent);
		_vertices = convert(_vertices, _parent->_inv);
		_parent = _parent->_parent;
	}

	void convert() {
		assert(_parent);
		_vertices = convert(_vertices, _parent->_conv);
	}
  
  template<class C>
  cfg_myset convert(cfg_myset s, C const& conv) {
    cfg_myset result;
    for (auto v : s) {
      assert(conv[v] >= 0);
      result.insert(conv[v]);
    }
    return result;
  }

	template<class T>
	void collectBagsToPack(T list, bag_t const* from) {
		for (bag_t* bag : _incidentBags) {
			if (bag != from) {
				bag->collectBagsToPack(list,  this);
			}else{
			}
		}
	}

#if 0
	void figureOutSafety(SafeSeparator const& ss) {
		if (!safe && !unsafe) { untested();
			safe = ss.isSafeSeparator(_vertices);
			unsafe = !safe;
		}else{
		}
	}
#endif
  
	void figureOutSafetyBySPT() {
    if (!safe && !unsafe) { untested();
      safe = isSafe();
      unsafe = !safe;
    }else{ untested();
	 }
  }
  
   bool isSafe() {
    return isSafeBySPT();
  }
  
  bool isSafeBySPT() {
	  assert(_graph);
    parentVertex.resize(boost::num_vertices(*_graph));
    //ArrayList<XBitSet> components = _graph.getComponents(vertexSet);

	 auto mask = _vertices; //  treedec::util::make_incidence_mask(_vertices);
	 auto vert = boost::vertices(_graph);
    auto cmps_range = treedec::make_components_range(vert.first, vert.second, _graph, mask);
	 for(; cmps_range.first!=cmps_range.second; ++cmps_range.first){
		 auto comp_range=*cmps_range.first;
		 cfg_myset compo;
		 assert(compo.empty());
		 for(; comp_range.first!=comp_range.second; ++comp_range.first){
			 compo.add(*comp_range.first);
		 }
      if (!isSafeComponentBySPT(compo)) {
        return false;
      }else{
		}
    }
    return true;
  }
  
	bool isSafeComponentBySPT(cfg_myset const& component) const {
		cfg_myset neighborSet = _graph.neighborSet(component);
		cfg_myset rest = _graph.all.subtract(neighborSet).subtract(component);

		for (auto v : neighborSet) {
			cfg_myset missing = neighborSet;
			missing = subtract(_graph.neighborSet[v]);

			for (auto w : missing) {
				missing.erase(w);
				if(w>v){
					break;
				}else{
				}
			}

			if (missing.empty()) {
			}else{
				cfg_myset spt = shortestPathTree(v, missing, rest);
				if (spt.empty()) {
					return false;
				}else{
					rest.subtract(spt);
				}
			}
		}
		return true;
	}

  cfg_myset shortestPathTree(int v, cfg_myset const& targets, cfg_myset const& available) {
    cfg_myset aunion = available;
	 aunion.merge(targets);
    
    cfg_myset reached;
    reached.insert(v);
    cfg_myset leaves(reached);
    while (!targets.is_subset_of(reached) && !leaves.empty()) {
      cfg_myset newLeaves;
      for (auto u : leaves){
        cfg_myset children = _graph.neighborSet[u].intersectWith(aunion).subtract(reached);
        for (auto w : children) {
          reached.insert(w);
          parentVertex[w] = u;
          if (available.contains(w)) {
            newLeaves.insert(w);
          }
        }
      }
      leaves = newLeaves;
    }
    
   cfg_myset spt;

	if (targets.is_subset_of(reached)) {
		for (auto u : targets){
			int w = parentVertex[u];
			while (w != v) {
				spt.insert(w);
				w = parentVertex[w];
			}
		}
	}else{
	 }
    return spt;
  }


  // void dump(String indent) {
  //   System.out.println(indent + "sep:" + toString());
  // }
  
  std::ostream& print(std::ostream& sb) {
    sb << _vertices;
    sb << "<";
    for (bag_t const* bag: _incidentBags){
      if (!bag) {
        sb << ("null bag ");
      } else {
        sb << _parent.nestedBags.indexOf(bag) << ":" << bag._vertices;
        sb << (" ");
      }
    }
    sb << ")";
    
    return sb;
  }
  
};

#endif

