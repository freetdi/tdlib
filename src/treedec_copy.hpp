#pragma once

namespace treedec {

namespace hack{

struct forgetprop
{
  template <class G, class H>
  void operator()(G, H) const
  {
  }
};

}

// a hack.
//
// boost::copy is (somehow) unable to convert bags.
// this is required to hold some old code together.
//
// we are making all possible assumptions here.
template<class S, class T>
void obsolete_copy_treedec(S const& src, T& tgt)
{

	  boost::copy_graph(src, tgt,
				 boost::vertex_copy(hack::forgetprop()).
				 edge_copy(hack::forgetprop()));

	auto p=boost::vertices(src);
	auto q=boost::vertices(src);
	for(;p.first!=p.second; ){

		auto const& bs=boost::get(bag_t(), src, *p.first);
		auto& bt=boost::get(bag_t(), tgt, *q.first);

		push(bt, bs.begin(), bs.end());

		++p.first;
		++q.first;
	}

}

}
