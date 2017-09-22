// Lukas Larisch, 2014 - 2016
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//

#ifndef TD_RANDOM_GENERATORS
#define TD_RANDOM_GENERATORS

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>
#include "trace.hpp"


namespace treedec{

namespace random{

static bool coin(){
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

inline unsigned randint(){
    static boost::random::mt11213b _rnd_gen; //fastest according to boost reference.
    static boost::random::uniform_int_distribution<> _dist;
    static unsigned _rnd;

    _rnd=_dist(_rnd_gen);

    return _rnd;
}

}

}

#endif
