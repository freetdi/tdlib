// Felix Salfelder, 2017
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
//
//
// inspired by STL_ITERATOR_EAN20051028_HPP
//             Eric Niebler 2005.

#ifndef PYITER_H
#define PYITER_H

#include <boost/python/detail/prefix.hpp>
#include <boost/python/object/stl_iterator_core.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "exception.hpp"

namespace detail {

struct dummy_checker{
	template<class T>
	bool operator()(T) const {return true;}
};

// An STL input iterator over a python sequence
template<typename ValueT, typename checker=dummy_checker>
struct pyedge_input_iterator
  : boost::iterator_facade<pyedge_input_iterator<ValueT, checker>,
                           ValueT,
									std::input_iterator_tag,
									ValueT>
{
    pyedge_input_iterator()
      : impl_(), _ok(checker())
    {
    }

    pyedge_input_iterator(boost::python::object const &ob, checker c=dummy_checker())
      : impl_(ob), _ok(c)
    {
    }

private:
    friend class boost::iterator_core_access;

    void increment() {
		 trace0("inc");
        this->impl_.increment();
    }

	 ValueT dereference() const {
		 namespace py = boost::python;

		 auto const& what=this->impl_.current().get();
		 auto e=boost::python::extract<boost::python::tuple>(what)();
		 // trace1("edg", std::string(py::extract<std::string>(py::str(e))()));
		 typedef typename ValueT::first_type ft;
		 auto a=boost::python::extract<ft>(e[0])();
		 auto b=boost::python::extract<ft>(e[1])();

		 if(!_ok(a)){
			 throw exception_invalid("source invalid: "
					 + py::extract<std::string>(py::str(e))());
		 }else if(!_ok(b)){ untested();
			 throw exception_invalid("target invalid: "
					 + py::extract<std::string>(py::str(e))());
		 }else{
			 trace0("makepair");
			 return std::make_pair(a, b);
		 }
	 }

    bool equal(pyedge_input_iterator const &that) const {
        return this->impl_.equal(that.impl_);
    }

	 boost::python::objects::stl_input_iterator_impl impl_;
	 checker _ok;
};

}

#endif
