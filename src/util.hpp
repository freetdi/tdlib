#ifndef TREEDEC_UTIL_HPP
#define TREEDEC_UTIL_HPP

namespace treedec{

template <typename S>
bool is_permutation(S const &s)
{ untested();

	std::set<typename S::value_type> X(s.begin(), s.end());

	return (s.size()==X.size());
}

// check if s is a permutation of t
template <typename S, typename T>
bool is_permutation(S const &s, T const &t)
{ untested();

	if(s.size()!=t.size()){ untested();
		return false;
	}else{ untested();

		std::set<typename T::value_type> X, Y;

		for( auto x : s){ untested();
			X.insert(x);
		}
		if(s.size()!=X.size()){
			return false;
		}

		for( auto y : t){ untested();
			Y.insert(y);
		}

		return X == Y;
	}
}

}

#endif
