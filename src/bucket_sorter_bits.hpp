
#ifndef BUCKET_SORTER_DETAILS
#define BUCKET_SORTER_DETAILS

#include "container_traits.hpp"
#include "bucket_sorter.hpp"

namespace treedec {

namespace detail{

    template<class C>
      struct container_iter<C,
      typename std::enable_if<std::is_same<typename C::base,
      boost::bucket_sorter<
        typename C::base::bucket_type,
      typename C::base::value_type,
      typename C::base::Bucket_type,
      typename C::base::value_index_map > >::value >::type >
      {//

        // front and back. pass elements, not references.
        // use begin/rbegin if you need references.
        static typename C::value_type front(C const& c)
        { untested();
          return *c.begin();
        }
        static typename C::value_type back(C const& c)
        { untested();
          // does not have back access.
          return *c.begin();
        }
      };

} // detail

} // treedec

#endif

// vim:ts=8:sw=2:et
