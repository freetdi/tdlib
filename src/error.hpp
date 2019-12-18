// Lukas Larisch, 2016
// Felix Salfelder, 2016
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
//
//

#ifndef TREEDEC_ERROR_HPP
#define TREEDEC_ERROR_HPP

namespace treedec{

class exception_unsuccessful : public std::runtime_error{
public:
  exception_unsuccessful() : std::runtime_error("exception_unsuccessful") { }
};

class exception_invalid_precondition : public std::runtime_error{
public:
  exception_invalid_precondition() : std::runtime_error("exception_invalid_precondition") { }
};

}
#endif
