// Felix Salfelder, 2021
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
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


// hacks to avoid boost incompatibilities.
// perhaps this should be version dependent.

#ifndef BOOST_COMPAT_H
#define BOOST_COMPAT_H

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;
#define BOOST_BIND_HPP_INCLUDED // avoid deprecated boost/bind.hpp include

#include <iterator>
#define ITERATOR_DWA122600_HPP_ // avoid deprecated header

namespace boost{
	namespace detail{
		using std::iterator_traits;
	}
}

#endif
