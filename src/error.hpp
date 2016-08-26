// Lukas Larisch, 2016
// Felix Salfelder, 2016
//
// (c) 2016 Goethe-Universität Frankfurt
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
//

#ifndef TD_ERROR_HPP
#define TD_ERROR_HPP

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
