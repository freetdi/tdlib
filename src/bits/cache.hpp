
#ifndef TREEDEC_CACHE_H
#define TREEDEC_CACHE_H
template<class Key, class Block, class compare>
class BlockCacheMap{
	typedef Key myset;
private:
	class blockcompare{
	public:
		using is_transparent = void;
		bool operator()(Block const* a, Block const* b) const {
			return a->key()==b->key();
		}
		bool operator()(myset const& a, myset const& b) const {
			return a == b;
		}
		bool operator()(myset const& a, Block const* b) const {
			return a == b->key();
		}
		bool operator()(Block const* a, myset const& b) const {
			return a->key() == b;
		}
	};
	class blockless{
	public:
		using is_transparent = void;
		bool operator()(Block const* a, Block const* b) const {
			return compare()(a->key(), b->key());
		}
		bool operator()(myset const& a, myset const& b) const {
			return compare()(a, b);
		}
		bool operator()(myset const& a, Block const* b) const {
			return compare()(a, b->key());
		}
		bool operator()(Block const* a, myset const& b) const {
			return compare()(a->key(), b);
		}

	};
	typedef std::set<Block const*, blockless> map_t;
public:
	BlockCacheMap(BlockCacheMap const&) = delete;
	explicit BlockCacheMap(){}

public:

#if 0
	template<class A, class B>
	std::pair<Block const*, bool> cache(myset const& key, A const& a, B const& b){
		auto it = _s.find(key);
		if (it != _s.end() && (*it)->key() == key){
//		  	assert((*it)->key() == key);
			return std::make_pair(*it, false);
		} else {
			auto x = new Block(key, a, b); // TODO placement new
			Block const* X = *_s.emplace_hint(it, x);
			assert(X->key() == key);
			return std::make_pair(X, true);
		}

	}
#endif
	template<class A, class B>
	std::pair<Block const*, bool> cache(myset const& key, A const& a, B const& b){
		auto it = _s.lower_bound(key);
		if (it != _s.end() && (*it)->key() == key){
			return std::make_pair(*it, false);
		} else {
			auto x = new Block(key, a, b); // TODO placement new
			Block const* X = *_s.emplace_hint(it, x);
			assert(X->key() == key);
			return std::make_pair(X, true);
		}

	}

public:
	bool has(Key const& k) const{
		auto ii = _s.find(k);
		return ii != _s.end();
	}
	bool empty() const{
		return _s.empty();
	}
	void clear(){
		for(auto i : _s){
			delete i;
		}
		_s.clear();
	}
	size_t size() const{
		return _s.size();
	}

	
	Block * find(myset const* key) {
		// BUG //
		incomplete();
		return const_cast<Block*>(*_s.find(key));
	}
	Block * find(myset const& key) {
		// BUG //
		return const_cast<Block*>(*_s.find(key));
	}

	auto end(){
		return *_s.end();
	}

public:
	map_t const& s() const{return _s;}
private:
	map_t _s;
};

template<class Key, class Block>
class BlockCacheHash{
	typedef Key myset;
private:
	class blockcompare{
	public:
		using is_transparent = void;
		bool operator()(Block const* a, Block const* b) const {
			return a->key()==b->key();
		}
		bool operator()(myset const& a, myset const& b) const {
			return a == b;
		}
		bool operator()(myset const& a, Block const* b) const {
			return a == b->key();
		}
		bool operator()(Block const* a, myset const& b) const {
			return a->key() == b;
		}
	};
	class blockhash{
	public:
		using is_transparent = void;
		using transparent_key_equal = blockcompare;
		size_t operator()(Block const* const& a) const{
			return operator()(a->key());
		}
		size_t operator()(myset const& a) const{
			size_t hash=0;
			for(unsigned k = 0; k< a.howmany(); ++k){
				hash ^= a.chunk(k);
			}

//			unsigned numbytes = (myset::max_element + 1) / 8;
//			char* data = (char*) &a.chunk(0);
//			for(unsigned k = 0; k< numbytes; ++k){
//				hash ^= data[k];
//			}

			return hash;
		}
	};
	class blockless{
	public:
		using is_transparent = void;
		bool operator()(Block const* a, Block const* b) const {
			return setless(a->key(), b->key());
		}
		bool operator()(myset const& a, myset const& b) const {
			return setless(a, b);
		}
		bool operator()(myset const& a, Block const* b) const {
			return setless(a, b->key());
		}
		bool operator()(Block const* a, myset const& b) const {
			return setless(a->key(), b);
		}

	};
	// typedef std::set<Block const*, blockless> map_t;
	typedef std::unordered_set<Block const*, blockhash, blockcompare> map_t;
public:
	BlockCacheHash(BlockCacheHash const&) = delete;
	explicit BlockCacheHash(){}

public:
	template<class A, class B>
	std::pair<Block const*, bool> cache(myset const& key, A const& a, B const& b){
		auto it = _s.find(key);
		if (it != _s.end() && (*it)->key() == key){
//		  	assert((*it)->key() == key);
			return std::make_pair(*it, false);
		} else {
			auto x = new Block(key, a, b); // TODO placement new
			Block const* X = *_s.emplace_hint(it, x);
			assert(X->key() == key);
			return std::make_pair(X, true);
		}

	}

#if 0
	template<class A, class B>
	std::pair<Block const*, bool> cache_set(myset const& key, A const& a, B const& b){
		auto it = _s.lower_bound(key);
		if (it != _s.end() && (*it)->key() == key){
//		  	assert((*it)->key() == key);
			return std::make_pair(*it, false);
		} else {
			auto x = new Block(key, a, b); // TODO placement new
			Block const* X = *_s.emplace_hint(it, x);
			assert(X->key() == key);
			return std::make_pair(X, true);
		}

	}
#endif

public:
	void clear(){
		for(auto i : _s){
			delete i;
		}
		_s.clear();
	}
	bool has(Key const& k) const{
		auto ii = _s.find(k);
		return ii != _s.end();
	}
	bool empty() const{
		return _s.empty();
	}
	size_t size() const{
		return _s.size();
	}

public:
	map_t const& s() const{return _s;}
private:
	map_t _s;
};
#endif
