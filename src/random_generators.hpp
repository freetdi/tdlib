#ifndef TD_RANDOM_GENERATORS
#define TD_RANDOM_GENERATORS

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>
#include "trace.hpp"


namespace treedec{

namespace random{

static bool coin()
{
    static boost::random::mt11213b _rnd_gen; //fastest according to boost reference.
    static boost::random::uniform_int_distribution<> _dist;
    static unsigned _rnd, _which;

    if(!_which){
        _rnd=_dist(_rnd_gen);
        _which=1;
    }
    bool c = _rnd & _which;
    _which <<= 1;

    return c;
}

inline unsigned randint(){ untested();
    static boost::random::mt11213b _rnd_gen; //fastest according to boost reference.
    static boost::random::uniform_int_distribution<> _dist;
    static unsigned _rnd;

    _rnd=_dist(_rnd_gen);

    return _rnd;
}

}

}

#endif
