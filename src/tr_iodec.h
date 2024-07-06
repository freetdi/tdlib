/* (c) Felix Salfelder 2021
 * License: GPLv3+
 *
 * Derived from pace_meiji (c) 2017, Hisao Tamaki, MIT license
 */

#include "tr_sieve.hpp"
#include "tr_sep.h"
#include <queue>

//static bool VERBOSE = true;
//static bool DEBUG = false;

template<class G>
class Bag;

struct firstcmp {
	bool operator()(myset const& a, myset const& b) const {
		assert(a.size());
		assert(b.size());

		return *a.begin() < *b.begin();
	}
};
struct mysetless {
	bool operator()(myset const& a, myset const& b) const {
			return a.compare_int(b) < 0;
		int l1 = a.size();
		int l2 = b.size();
		if (l1 != l2) { itested();
			return l1 < l2;
		}else{ itested();
			return a.compare_int(b) < 0;
		}
	}
}setless;
#include "bits/cache.hpp"


// std::ostream& operator<<(std::ostream& o, myset const& s)
// {
// 	o<< "set " << s;
// 	return o;
// }

template<class G>
class IODecomposer;

template<class G>
myset neighborSet(myset const& set, G const& g) 
{
    myset result;
    for (auto v : set){
      result.merge(g.out_edges(v));
    }
    result.subtract(set);
    return result;
}

// BAG decomposer
template<class G>
class IODecomposer {
	typedef Bag<G> bag_t;
	typedef LayeredSieve<myset, myset const*, 256> blocksieve_t;

	class Block {
	public: // needed for std::map
	//	Block() { // incomplete();
	//	}
	// 	explicit Block( Block const& o )
	// 	  :_component(o._component),
	// 	   _outbound(o._outbound),
	// 	   _separator(o._separator)
	// 	  	  { untested();
	// 			 untested();
	// 		  }
	public:
		template<class A>
		explicit Block(myset const& component, G const& g, A const& all)
		   : _component(component) {
			_separator = neighborSet(component, g);

			myset rest = all;
			rest.subtract(component);
			rest.subtract(_separator); // counter?

			assert(_component.size());
			int minCompo = *component.begin();

			// the scanning order ensures that the first full component
			// encountered is the outbound one
			for (auto v : rest) {
				myset c = g.out_edges(v);
				myset toBeScanned = c;
				toBeScanned.subtract(_separator);
				assert(!c.contains(v));
				c.insert(v);
				while (!toBeScanned.empty()) {
					myset save = c;
					for (auto w : toBeScanned){
						c.merge(g.out_edges(w));
					}
					toBeScanned = c;
					toBeScanned.subtract(save);
					toBeScanned.subtract(_separator);
				}
				if (_separator.is_subset_of(c)) {
					// full block other than "component" found
					if (int(v) < minCompo) {
						_outbound = c;
						_outbound.subtract(_separator);
					} else {
						// v > minCompo
						_outbound = _component;
					}
					return;
				}else{
					rest.subtract(c);
				}
			}
		}

		size_t separator_size() const{
			return _separator.size();
		}

		bool isOutbound() const {
			return _outbound == _component;
		}

		bool ofMinimalSeparator() const {
			return !_outbound.empty();
		}

		std::ostream& print(std::ostream& sb) const { untested();
			if (_outbound == _component) { untested();
				sb << "o";
			} else { untested();
//				if (iBlockCache.get(component) != null) { untested();
//					sb << "f";
//				} else { untested();
//					sb << "i";
//				}
			}
			sb << _component << "(" << _separator << ")";
			return sb;
		}

		int compareTo(Block const& b) const{ untested();
			assert(_component.size());
			assert(_component.size() == _component.recount());
			assert(b._component.size());
			assert(b._component.size() == b._component.recount());
			return int(_component.front()) - int(b._component.front());
		}

		myset const& key() const{
			return _component;
		}
		size_t size() const{
			return _component.size();
		}
		myset const& component() const{
			return _component;
		}
		myset const& outbound() const{
			return _outbound;
		}
		myset const& separator() const{
			return _separator;
		}
	private:
		myset _component; // const?
		myset _outbound; // const?
		myset _separator; // const?
	}; // Block

	/// PMC := max clique in some chordal completion of G
	//   <=> is cliquish and has no-full component
	class PMC {
		PMC(const PMC&) = delete;
	public:
		template<class BL>
		explicit PMC(myset const& v, BL const& blockList, G const& g) :
			_vertices(v) {
//			trace2("DEBUG PMC", v.size(), blockList.size());
			if (_vertices.empty()) { untested();
				return;
			}else{
			}
			for (Block const* block: blockList) {
				assert(block);
//				trace2("DEBUG PMC", block->component(), block->isOutbound());

				if (!block->isOutbound()){
				}else if(_outbound_block == nullptr){
					_outbound_block = block;
				}else if(_outbound_block->separator().is_subset_of(block->separator())){
					_outbound_block = block;
				}else{ untested();
				}
			}

			_inbounds.resize(0); // dropping refs.
			if (_outbound_block == nullptr) {
				_inbounds = blockList;
			}else{
				// std::cerr << "OUT " << _outbound_block->size() << "\n";
			//	_inbounds = new Block[blockList.size()];
				_inbounds.reserve(blockList.size());
//				int k = 0;
				for (Block const* block: blockList) {
//					assert(block->separator().size());
					if (!block->separator().is_subset_of(_outbound_block->separator())) {
						assert(block);
						_inbounds.push_back( block);
					}else{
					}
				}
			}
			checkValidity(g);

#if 0
				trace1("DEBUG: pmc created from", v);
				trace1("DEBUG:", _vertices);
				trace1("DEBUG:", _inbounds.size());
				//System.out.println("PMC created:");
				//System.out.println(this);
			trace0("DEBUG PMC");
#ifdef DO_TRACE
			print(std::cerr);
#endif
#endif
		} // ::PMC

		// PMC::
		void checkValidity(G const& g) {
			for (Block const* b: _inbounds) {
				assert(b);
				if (!b->ofMinimalSeparator()) {
					_isValid = false;
					return;
				}else{
				}
			}

			for (auto v : _vertices) {
				myset rest = _vertices;
				rest.subtract(g.out_edges(v));
				rest.erase(v);
				if (!_outbound_block){
				}else if( _outbound_block->separator().contains(v)) {
					rest.subtract(_outbound_block->separator());
				}
				for (Block const* b : _inbounds) {
					assert(b);
					if (b->separator().contains(v)) {
						rest.subtract(b->separator());
					}else{
					}
				}
				if (!rest.empty()) {
					_isValid = false;
					return;
				}else{
				}
			}
			_isValid = true;
		}

		// PMC::
#if 0
		bool isReady() const{ untested();
			for (auto ii : _inbounds) { untested();
				if (_iBlockCache.find(ii.component) == _iBlockCache.end()) { untested();
					return false;
				}else{ untested();
				}
			}
			return true;
		}
#endif

		// PMC::
#if 1
		myset getTarget() const {
			assert(_outbound_block);
			{
				myset combined = _vertices;
				combined.subtract(_outbound_block->separator());
				for (Block const* b: _inbounds) {
					combined.merge(b->component());
				}
				return combined;
			}
		}
#endif



		// PMC::
		void carryOutDecomposition(bag_t const& b) const {
			incomplete();
			(void)b;
#if 0
			auto const* bag = &b;
			trace1("carryOutDecomposition:", this);

			for (Block const* inbound: _inbounds) { untested();
				trace1("inbound  = ", inbound);
				auto ibi = _iBlockCache.find(inbound.component);
				if (ibi == _iBlockCache.end()) { untested();
					trace1("inbound iBlock is null, block = ", inbound);
					continue;
				}else{ untested();
					IBlock* iBlock = *ibi;
					assert(iBlock);

					Bag subBag = _rootBag.addNestedBag( iBlock.endorser._vertices);
					Separator separator = _rootBag.addSeparator(inbound.separator);

					separator.addIncidentBag(bag);
					separator.addIncidentBag(subBag);

					bag.addIncidentSeparator(separator);
					subBag.addIncidentSeparator(separator);

					iBlock.endorser.carryOutDecomposition(subBag);
				}
			}
#endif
		}

		myset inletsInduced() const { untested();
			myset result;
			for (Block b : _inbounds) { untested();
				result.merge(b.separator);
			}
			return result;
		}

		// PMC::
		std::ostream& print(std::ostream& sb) const{ itested();
			sb << "PMC";
			if (_isValid) { itested();
				sb << "(valid):\n";
			}else{ itested();
				sb << "(invalid):\n";
			}
			sb << "  sep     : " << _vertices << "\n";
			sb << "  outbound: " << _outbound_block << "\n";

			for (Block const* b : _inbounds) { itested();
				sb << "  inbound : " << b << "\n";
			}
			return sb;
		}

		bool is_valid() const{ return _isValid; }
		myset const& vertices() const{ return _vertices; }
		std::vector<Block const*> const& inbounds() const{
			return _inbounds;
		}
		Block const* outbound_block() const{
		  return _outbound_block;
		}

		size_t size() const{ return _vertices.size(); }

	private:
		myset _vertices;
		std::vector<Block const*> _inbounds;
		Block const* _outbound_block{nullptr};
		bool _isValid{false};
	}; // PMC

	void endorse(PMC const* p) {
		assert(p);
//		trace1("DEBUG, endorsing ", p);

		if (!p->outbound_block()) {
			trace0("DEBUG solution found in endorse()");
			_solution = p;
		} else {
//			trace1("DEBUG", p->outbound_block()->component());
//			trace1("DEBUG", p->outbound_block()->separator());
			endorse(p, p->getTarget());
		}
	}

	void endorse(PMC const* p, myset const& target) {
		// std::cerr << "END " << p->size() << " " << target.size() << "\n";
		// if (separator.equals(bs1)) { untested();
		// System.err.println("endorsed = " + endorsed +
		// ", " + endorserMap.get(endorsed));
		// }
		//

		auto x = _iBlockCache.cache(target, this, p);
		bool inserted = x.second;
		if (inserted) {
			// std::cerr << "RQP " << *x.first->block().component().begin() << "\n";
			_readyQueue.push(x.first);
		}else{
			// std::cerr << "RQP miss " << *target.begin() << "\n";
		}
	}



	class IBlock {
	public: // needed for std::map
		// IBlock() {
		//   	incomplete();
		// }
	public:
		IBlock(IBlock const& o)
		  : _block(o._block), _endorser(o._endorser) { untested();
		  }
		explicit IBlock(Block const* block, PMC const* endorser)
			: _block(block),
			_endorser(endorser){ untested();
		}
		explicit IBlock(myset const& target, IODecomposer* d, PMC const* endorser)
			: _block(&d->getBlock(target)),
			_endorser(endorser){
		}

		std::ostream& print(std::ostream& sb) { untested();
			sb << "IBlock:" + _block.separator + "\n";
			sb << "  in  :" + _block.component + "\n";
			sb << "  out :" + _block.outbound + "\n";
			return sb;
		}

		Block const& block() const{
			assert(_block);
			return *_block;
		}
		myset const& key() const{
			assert(_block);
			return _block->key();
		}

	private:
		Block const* _block{nullptr};
		PMC const* _endorser{nullptr};
	};

private:

	bool isReady(PMC const& p) const{
		for (auto ii : p.inbounds()) {
			if (!_iBlockCache.has(ii->component())) {
				return false;
			}else{
			}
		}
		return true;
	}
  
//  SafeSeparator ss;

  // bag lives outside of Decomposer...
public:
	explicit IODecomposer(Bag<G>& bag, int lowerBound=0, unsigned upperBound=-1)
	  : _graph(bag.graph()),
       _rootBag(bag),
	    _lowerBound(lowerBound), _upperBound(upperBound)
	{
		assert(boost::num_vertices(_graph));
		// _rootBag = bag; // root bag?
#ifndef NDEBUG
		auto r = treedec::make_components_range(_graph);
		assert(r.first!=r.second); // empty?
		++r.first;
		assert(r.first==r.second);
#endif
		// not used in exact?
		//    ss = new SafeSeparator(g);
		_all_vertices = bag.vertices();
		assert(boost::num_vertices(_graph) == _all_vertices.size());
	}

	void process(IBlock const* ibl);

	void makeSimpleTBlock(IBlock const& ibl) {
		// assert(ibl);
		auto const& bl = ibl.block();
		auto pp = _oBlockCache.cache(bl.separator(), bl.outbound(), nullptr);
		if(pp.second){
//			auto& x = 
				_oBlockSieve.put(bl.outbound(), &pp.first->separator());
			// x = pp.first->separator();
			crown(*pp.first);
		}
	}

	// take ownership
	void process(PMC const* p) {
		assert(p);
//		trace1("processing ", p);
		if (isReady(*p)) { untested();
			trace1("endorsing ", p);
			endorse(p);
		}else{
			_pendingEndorsers.push_back(p);
		}
	}

  
	void decompose() {
		assert(_blockCache.empty());
		assert(_iBlockCache.empty());
		assert(_pendingEndorsers.empty());
		//    assert(_pmcCache.empty());

//		_targetWidth = 2;
		// target bag size.
		unsigned tbs = _targetWidth + 1;

		trace3("...", _targetWidth, _upperBound, _lowerBound);
		while (unsigned(_targetWidth) <= unsigned(_upperBound)) {
			trace2("...", tbs, _rootBag.size());

			if (unsigned(_rootBag.size()) <= tbs) {
				_rootBag._nestedBags.clear(); // BUG: memory leak
				_rootBag._separators.clear(); // BUG: memory leak
				return;
			}else{
			}

			// endorserMap = new HashMap<>();

			_oBlockSieve.init(_rootBag.size(), _targetWidth);
			_oBlockCache.clear(); //  = new HashMap<>();

			std::queue<IBlock const*> empty;
			std::swap( _readyQueue, empty );
				// std::cerr << "RQFs " << _iBlockCache.size() <<  "\n";
			for(auto i : _iBlockCache.s()){
				//_readyQueue.push(&i.second);

				// std::cerr << "RQF " << i->block().component().size() << " " <<
				// *i->block().component().begin() << " " <<
				// *i->block().component().rbegin() << "\n";
				_readyQueue.push(i);
			}
				// std::cerr << "RQF.\n";

			for (auto v_ : _all_vertices){
				unsigned v = v_;
				myset cnb(_graph.out_edges(v));

//				trace2("loop", v, cnb.size());
				assert(!cnb.contains(v));

				if (unsigned(cnb.size()) > _targetWidth) {
					continue;
				}else{
				}
				cnb.insert(v);
				assert(cnb.contains(v));

				auto bl = getBlocks(cnb);
//				trace2("loop", bl.size(), cnb.size());
//
				//      if (!pmcCache.contains(cnb)) { untested();
				PMC* pmc = new PMC(cnb, bl, _graph);
				if (!pmc->is_valid()) {
					//          pmcCache.add(cnb);
				}else if (isReady(*pmc)) {
					endorse(pmc);
//						delete pmc; // ?
				}else{
					_pendingEndorsers.push_back(pmc);
				}
					//        }
			}

			trace1("entering loop", _readyQueue.size());
			trace2("PES", _targetWidth, _pendingEndorsers.size());

			while (true) {
				while (!_readyQueue.empty()) {

					IBlock const* ready = _readyQueue.front();
					assert(ready);
					_readyQueue.pop();

					process(ready); // (_graph);
					// delete ready?

					if (_solution) {
						stats(std::cout, "solution1: ") << "\n";
						bag_t const& bag = _rootBag.addNestedBag(_solution->vertices());
						_solution->carryOutDecomposition(bag);
						return;
					}else{
					}
				}

				trace1("PES", _pendingEndorsers.size() );

				std::vector<PMC const*> endorsers;
				std::swap(endorsers, _pendingEndorsers);

				for (auto endorser : endorsers) {
					//endorser->process(*this);
					process(endorser);
					// incomplete();
					// delete endorser?!
					if (_solution) { untested();
						stats(std::cout, "solution2: ") << "\n";
						bag_t const& bag = _rootBag.addNestedBag(_solution->vertices());
						_solution->carryOutDecomposition(bag);
						return; // LEAK
					}else{
					}
				}

				if (_readyQueue.empty()) {
					break;
				}else{ untested();
				}
			}

			stats(std::cout, "failed: ") << "\n";

			_targetWidth++;
			tbs++;
		}
		return;
	}
private:

#if 0
  bool crossesOrSubsumes(XBitSet separator1, XBitSet endorsed, XBitSet separator2) { untested();
    ArrayList<XBitSet> components = g.getComponents(separator1);
    for (XBitSet compo: components) { untested();
      if (endorsed.is_subset_of(compo)) { untested();
        // subsumes
        return true;
      }
    }
    // test crossing
    myset diff = separator2;
	 diff.subtract(separator1);
    for (XBitSet compo: components) { untested();
      if (diff.is_subset_of(compo)) { untested();
        return false;
      }
    }
    return true;
  }
#endif

	Block const& getBlock(myset const& component) {
		auto p = _blockCache.cache(component, _graph, _all_vertices);
		return *p.first;
	}

	bool isFullComponent(myset const& component, myset const& sep) { untested();
		for (auto v: sep){ untested();
			if (!component.intersects(_graph.out_edges(v))) { untested();
				return false;
			}else{ untested();
			}
		}
		return true;
	}

	// not const. calls getBlock...
	std::vector<Block const*> getBlocks(myset const& separator) {
		std::vector<Block const*> result;
		myset rest = _all_vertices; // bag has them..?
		rest.subtract(separator);
		while ( !rest.empty() ){
			auto v = *rest.begin();
			trace2("loop", v, rest.size());
			myset c(_graph.out_edges(v));
			c.subtract(separator);
			myset toBeScanned(c);
			c.insert(v);
			while (!toBeScanned.empty()) {
				myset save(c);
				for (auto w : toBeScanned) { itested();
					c.merge(_graph.out_edges(w));
				}
				c.subtract(separator);
				//        toBeScanned = c - save;
				toBeScanned = c;
				toBeScanned.subtract(save);
			}

			Block const* block = &getBlock(c);
			result.push_back(block);
			assert(c.contains(v));
			rest.subtract(c);
			assert(!rest.contains(v));
		}
		return result;
	}



	class OBlock {
	public:
		explicit OBlock(myset const& separator, myset const& openComponent, void*unused=nullptr) 
			: _openComponent(&openComponent)
		{
			(void)unused;
			// check. needed?
			_separator = separator;
		}

		// OBlock

		std::ostream& print(std::ostream& sb) { untested();
			sb << "TBlock:\n"
				<< "  sep :" << _separator << "\n"
				<< "  open:" << _openComponent << "\n";
			return sb;
		}

	public:
		myset const& key() const{
			return _separator;
		}
		myset const& separator() const{
			return _separator;
		}
		myset const& openComponent() const{
			return *_openComponent;
		}
	private:
		myset _separator; // ??
		myset const* _openComponent; // pointer?
	}; // OBlock

	void crown(OBlock const& o) {
		for (auto v : o.separator()){
			myset i = _graph.out_edges(v);
			i.intersect(o.openComponent());

			myset newsep = o.separator();
			newsep.merge(i);

			int tbs = _targetWidth + 1;
			if (int(newsep.size()) <= tbs) {
				//          if (!pmcCache.contains(newsep)) {  
				auto pmc = new PMC(newsep, getBlocks(newsep), _graph); // ??
				if (!pmc->is_valid()) {
					//              pmcCache.add(newsep);
					delete pmc;
				}else if (isReady(*pmc)) {
					endorse(pmc);
					// delete?
				}else{
					_pendingEndorsers.push_back(pmc);
				}
				//            }
			}
		}
	}

	void plugin(OBlock const& o, IBlock const& ibl) {
		unsigned tbs = _targetWidth + 1;
		// std::cerr<< "PLG " << ibl.block().separator().size()  << " " << ibl.block().component().size() 
		  //                  << " " << o.separator().size()  << " " << o.openComponent().size()  << "\n";

		myset newsep = o.separator();
		newsep.merge(ibl.block().separator());

		if (newsep.size() > tbs){ itested();
			return;
		}else{
		}

		auto blockList = getBlocks(newsep);

		Block const* fullBlock = nullptr;
		unsigned nSep = newsep.size();

		for (Block const* block : blockList) {
			if (block->separator_size() == nSep) {
				if (fullBlock) {
					//             minimal separator: treated elsewhere
					return;
				}else{
					fullBlock = block;
				}
			}else{
			}
		}

		if (fullBlock == nullptr) {
			//        if (!pmcCache.contains(newsep)) {  
			auto pmc = new PMC(newsep, blockList, _graph);
			// std::cerr<< "PMC " << newsep.size()  << " " << blockList.size() <<  " " << pmc->inbounds().size() << "\n";

			if (!pmc->is_valid()) {
				//            pmcCache.add(newsep);
					delete pmc;
			}else if (isReady(*pmc)) {
				endorse(pmc);
					// delete?
			} else {
				_pendingEndorsers.push_back(pmc);
			}
				//          }
		}else{
			if (newsep.size() > _targetWidth) {
				return;
			}else{
			}

#if 1
			auto pp = _oBlockCache.cache(newsep, fullBlock->component(), nullptr);
			if(pp.second){
//				auto& x = 
					_oBlockSieve.put(fullBlock->component(), &pp.first->separator());
				// x = pp.first->separator();
				crown(*pp.first);
			}else{
			}
#else
			auto oBlockf = _oBlockCache.find(newsep);
			if (oBlockf == _oBlockCache.end()) { untested();
				OBlock ob(newsep, fullBlock->component());
				_oBlockCache[newsep] = ob;
				_oBlockSieve.put(fullBlock->component(), newsep);
				auto const& n = _oBlockCache[newsep];
				crown(n);
			}else{
			}
#endif
		}
	}

	int numberOfEnabledBlocks() { untested();
		return _iBlockCache.size();
	}

  void dumpPendings() { untested();
    trace0("pending endorsers\n");
    for (auto endorser : _pendingEndorsers) { untested();
      trace1("...", endorser);
    }
  }

public:
	std::ostream& stats(std::ostream& o, std::string h="") const {
		trace1("DEBUG", _targetWidth);
		o << h
		  << "n = " << _all_vertices.size()
		  << " width = " << _targetWidth
		  << ", oBlocks = " << _oBlockCache.size()
		  << " ";
		  _oBlockSieve.stats(o);
		o << " , endorsed = " << _iBlockCache.size();
		o << ", pendings = " << _pendingEndorsers.size();
		o << ", blocks = " << _blockCache.size();
		return o;
	}

public:
  void set_min_bs(unsigned u){ untested();
	  _lowerBound = u-1;
  }
  void set_max_bs(unsigned u){
	  _upperBound = u-1;
  }
private:
	G const& _graph;
	myset _all_vertices;
	Bag<G>& _rootBag; // rootbag?
//	LayeredSieve<myset, myset, 128> _oBlockSieve; // store oblocks here?
//	LayeredSieve<myset, myset, 256> _oBlockSieve; // store oblocks here?
	blocksieve_t _oBlockSieve; // store oblocks here?
	std::queue<IBlock const*> _readyQueue;
	std::vector<PMC const*> _pendingEndorsers;

	BlockCacheMap<myset, OBlock, mysetless> _oBlockCache; // only need myset keys?
	BlockCacheMap<myset, Block, mysetless> _blockCache;

//	BlockCacheHash<myset, IBlock> _iBlockCache;
	BlockCacheMap<myset, IBlock, mysetless> _iBlockCache; // faster than hash? ordering?

//	Set<myset> pmcCache; not used;
  
  unsigned _lowerBound{0};
  unsigned _upperBound{-1u}; // max possible tw ?
  
  unsigned _targetWidth{0};

  PMC const* _solution{nullptr};
}; // IODecomposer

template<class G>
class TR {
	typedef Bag<G> bag_t;

public:
	TR(G const& g) : _g(g) {

	}

	void do_it(){
		int mindegree = 1; // for now.

		Bag<G> bag(_g);
		assert(int(bag.size()) == int(boost::num_vertices(_g)));
      IODecomposer<G> mtd(bag, mindegree);
		mtd.set_max_bs(boost::num_vertices(_g));

		mtd.decompose();
	}

private:
	G const& _g;
}; // IODecomposer

template<class G>
inline void IODecomposer<G>::process(IBlock const* ibl)
{
	// std::cerr<< "PIB " << ibl->block().component().size() << " " << ibl->block().separator().size() << "\n";
	// std::cerr<< "=== " << *ibl->block().component().begin() << " " << *ibl->block().separator().begin() << "\n";
	// assert(ibl);
	makeSimpleTBlock(*ibl);

	Block const& bl = ibl->block();

	blocksieve_t::set_out_hack buf;
	_oBlockSieve.collectSuperblocks( bl.component(), bl.separator(), buf);

	// std::cerr<< "CS " << bl.separator().size()  << " " << bl.component().size() << " " << buf.size() << "\n";

	for (myset const* tsep : buf) {
		// sort them and compare better?
		assert(tsep);
		// why not store the OBlock in the sieve?
		auto i = _oBlockCache.s().find(*tsep);
		assert(i!=_oBlockCache.s().end());
		//			i->second.plugin(ibl, tbs);
		plugin(**i, *ibl);
	}
}

