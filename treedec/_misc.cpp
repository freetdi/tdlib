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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

#include "boost_compat.h"
#include "misc.hpp"
#include "graph_py.hpp"
#include <boost/python.hpp>

namespace py = boost::python;

BOOST_PYTHON_MODULE(_misc)
{
	py::def("check_treedec", &treedec::check_treedec<_balvvu, _balvvu_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_balvvu, _balsvu_treedec>);

	py::def("check_treedec", &treedec::check_treedec<_balsvu, _balvvu_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_balsvu, _balsvu_treedec>);

#ifdef HAVE_GALA_GRAPH_H
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu16, _balvvu_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu16, _gsgvvu16_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu16, _gsgvvu32_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu16, _gsgvvu64_treedec>);

	py::def("check_treedec", &treedec::check_treedec<_gsgvvu32, _balvvu_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu32, _gsgvvu16_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu32, _gsgvvu32_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu32, _gsgvvu64_treedec>);

	py::def("check_treedec", &treedec::check_treedec<_gsgvvu64, _balvvu_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu64, _gsgvvu16_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu64, _gsgvvu32_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_gsgvvu64, _gsgvvu64_treedec>);

	py::def("check_treedec", &treedec::check_treedec<_balvvu, _gsgvvu16_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_balvvu, _gsgvvu32_treedec>);
	py::def("check_treedec", &treedec::check_treedec<_balvvu, _gsgvvu64_treedec>);
#endif
}
