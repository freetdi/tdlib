#include <iostream>

#ifndef incomplete
#define incomplete() \
	std::cout << "incomplete " << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n"
#endif

#ifndef unreachable
#define unreachable() \
	std::cerr << "unreachable " << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n"
#endif

// interactively tested
#undef itested
#ifdef TRACE_ITESTED
#include <iostream>
#define itested() \
	std::cerr << "@@# itested \n@@@:" << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n"
#else
#define itested()
#endif

#undef untested
#ifdef TRACE_UNTESTED
#include <iostream>
#define untested() \
	std::cerr << "@@# untested \n@@@:" << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n"
#define untested1(x) \
	std::cerr << "@@# untested \n###:" << x << "\n";
#else
#define untested()
#define untested1(x)
#endif
#define USE(x) (1)?(void)(0):(void)(x)

#undef trace0
#undef trace1
#undef trace2
#undef trace3
#undef trace4
#undef trace5

#ifdef DO_TRACE
#define trace0(s) ( std::cerr << "@#@" << (s) << "\n")
#define trace1(s,x) ( \
		std::cerr <<  "@#@" << (s) << "  " << #x << "=" << (x)  \
		     << std::endl )
#define trace2(s,x,y) ( \
		std::cerr <<  "@#@" << (s) << "  " << #x << "=" << (x)  \
		     << "  " << #y << "=" << (y)  \
		     << std::endl )
#define trace3(s,x,y,z) ( \
		std::cerr <<  "@#@" << (s) << "  " << #x << "=" << (x)  \
		     << "  " << #y << "=" << (y)  \
		     << "  " << #z << "=" << (z)  \
		     << std::endl )
#define trace4(s,x,y,z,w) ( \
		std::cerr <<  "@#@" << (s) << "  " << #x << "=" << (x)  \
		     << "  " << #y << "=" << (y)  \
		     << "  " << #z << "=" << (z)  \
		     << "  " << #w << "=" << (w)  \
		     << std::endl )
#define trace5(s,x,y,z,w,u) ( \
		std::cerr <<  "@#@" << (s) << "  " << #x << "=" << (x)  \
		     << "  " << #y << "=" << (y)  \
		     << "  " << #z << "=" << (z)  \
		     << "  " << #w << "=" << (w)  \
		     << "  " << #u << "=" << (u)  \
		     << std::endl )
#else
#define trace0(s) (USE(s))
#define trace1(s,x) (USE(s),USE(x))
#define trace2(s,x,y) USE(s);USE(x);USE(y)
#define trace3(s,x,y,z) USE(s);USE(x);USE(y);USE(z)
#define trace4(s,x,y,z,w) USE(s);USE(x);USE(y);USE(z);USE(w)
#define trace5(s,x,y,z,w,u) USE(s);USE(x);USE(y);USE(z);USE(w);USE(u)
#endif

#ifndef TRACE_H
#define TRACE_H
# ifndef NDEBUG
#include <vector>
template <class C>
inline bool contains(C const& c, typename C::value_type const& v)
{
	return(c.find(v) != c.end());
}
template <class X>
inline bool contains(std::vector<X> const& c, X const& v)
{
	typedef typename std::vector<X>::const_iterator it;

	for(it i=c.begin(); i!=c.end(); ++i){
		if(v == *i) return true;
	}
	return false;
}
# endif
#endif
