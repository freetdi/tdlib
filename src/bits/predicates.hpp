
#pragma once

#include <gala/cbset.h>
namespace util{ // BUG: "draft"
/*--------------------------------------------------------------------------*/
template<class S>
struct not_in_set{
	not_in_set(S const& s) : _s(s){}
public:
	bool operator()(typename S::value_type const& x) const{
		return !cbset::contains(_s, x);
	}
private:
	S const& _s;
};
/*--------------------------------------------------------------------------*/
template<class S>
not_in_set<S> make_not_in_set(const S& s)
{
	return not_in_set<S>(s);
}
/*--------------------------------------------------------------------------*/

template<class A, class B>
struct conj{
	conj(A const& a, B const& b)
		: _a(a), _b(b) {}
public:
	template<class X>
	bool operator()(X const& x) const{
		return _a(x) && _b(x);
	}
private:
	A const& _a;
	B const& _b;
};

template<class A, class B>
conj<A, B>const make_conj(const A& a, const B& b)
{
	return conj<A, B>(a, b);
}

template<class S>
struct lt{
	lt(S const& val) : _s(val){}
public:
	template<class T>
	conj<lt, T> operator&&(T const& t) const{ untested();
		return conj<lt, T>(*this, t);
	}
public:
	template<class X>
	bool operator()(X const& x) const{
		return x<_s;
	}
private:
	S const& _s;
};

template<class S>
lt<S> make_lt(const S& s)
{
	return lt<S>(s);
}

} // util

