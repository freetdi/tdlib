// Felix Salfelder, 2021
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//

#ifndef TREEDEC_EXCEPTION_H
#define TREEDEC_EXCEPTION_H

#include <exception>
#include "trace.hpp"

class exception_invalid : public std::exception {
public:
	exception_invalid(const exception_invalid&) = default;
	exception_invalid(std::string m)
            : _m(m){ }
	const char *what() const throw() { itested();
            return _m.c_str();
	}
private:
	std::string _m;
};

#endif
