#ifndef TD_GENERIC_BASE_H
#define TD_GENERIC_BASE_H
// TODO: cleanup. sort out etc.
#include "bits/bool.hpp"
#include "algo.hpp"
#include "marker.hpp"
#include "config_traits.hpp"
#include <boost/iterator/filter_iterator.hpp>
//#include "generic_elimination_search_overlay.hpp"

namespace treedec{

namespace gen_search {

template<class A, class B>
class overlay;

// implements graph whitelist + elim + undo_elim.
// what is CFGT_t?
// rearrange...?!
template <typename G_t, class CFG_t, template<class G, class ...> class CFGT_t>
class generic_elimination_search_base
  : public treedec::algo::draft::algo1 {
	 typedef treedec::algo::draft::algo1 baseclass;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
protected:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
public:
    typedef CFGT_t<G_t> UC; // user config
    typedef overlay<G_t, G_t> default_overlay_type;
    typedef typename treedec::config::get::edge_overlay_graph<UC, default_overlay_type>::type
		               internal_graph_type;
    typedef typename internal_graph_type::adjacency_iterator overlay_adjacency_iterator;

//	typedef typename detail::eot<G_t, internal_graph_type>::type backend_type;
protected: // construct
    //TODO: better use iterators for elim_vertices
    generic_elimination_search_base(internal_graph_type&,
		                              std::vector<BOOL>& active, // here?
                                    std::vector<vd> &best_ordering_input,
                                    std::vector<vd> &current_ordering_input,
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth_input, unsigned nodes_generated_input,
                                    unsigned orderings_generated_input);

    generic_elimination_search_base(internal_graph_type&,
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth_input, unsigned nodes_generated_input,
                                    unsigned orderings_generated_input);

    generic_elimination_search_base(G_t const &g,
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth_input, unsigned nodes_generated_input,
                                    unsigned orderings_generated_input);

    ~generic_elimination_search_base(){
		 if(_need_cleanup & 1u){
			 delete &_active;
			 delete &_best_ordering;
			 delete &_current_ordering;
		 }
		 if(_need_cleanup & 2u){
			 delete &_g;
		 }
	 }
protected: // recursion.
    generic_elimination_search_base(generic_elimination_search_base& o);
public:

    virtual void do_it() = 0;

	 std::vector<vd> const& ordering() { return _best_ordering; }
    void elimination_ordering(std::vector<vd> &ordering) { ordering = _best_ordering; }

    unsigned global_lower_bound() const{ return _global_lb; }
    unsigned global_upper_bound() const{ return _global_ub; }

    unsigned get_nodes_generated() const{ return _nodes_generated; }
    unsigned get_orderings_generated() const{ return _orderings_generated; }
protected:
    struct active_filter{
        active_filter(std::vector<BOOL> const & v) : _v(v) {}
        bool operator()(vertex_descriptor v) const{
            return _v[v];
        }
        std::vector<BOOL> const& _v;
    };
    typedef boost::filter_iterator<active_filter, overlay_adjacency_iterator> adj_it;
    typedef std::pair<adj_it, adj_it> adj_range;
#if 0 //incomplete
    std::pair<boost::filter_iterator<active_filter, base_iterator>,
              boost::filter_iterator<active_filter, base_iterator> >
    vertex_range() const{ untested();
        active_filter p(Overlay._active);
        auto q=boost::vertices(Overlay);
        return std::make_pair(filter_iter_first(p, q.first, q.second),
                              filter_iter_first(p, q.second, q.second));
    }
#endif
	 adj_range adjacent_vertices(vertex_descriptor v) const;
protected:
	 // vertices_size_type?
	 unsigned eliminate(vertex_descriptor v);
	 void undo_eliminate(vertex_descriptor v);
	 const std::vector<BOOL> &active() const{
		 return _active;
	 }
	 std::vector<BOOL>& active(){
		 return _active;
	 }
protected:
    internal_graph_type& _g;
    std::vector<BOOL>& _active; // active and current_ordering -> numbering.
    std::vector<vd> &_best_ordering;
    std::vector<vd> &_current_ordering;

    unsigned _global_lb; //lb for the original graph
    unsigned _global_ub; //ub for the original graph

    unsigned _depth;

    unsigned _nodes_generated;
    unsigned _orderings_generated;

private:
    marker_type _marker;
    unsigned char _need_cleanup; // yuck. ugly
}; // generic_elimination_search_base

} // gen_search

} // treedec

namespace boost {

template<class A, class O, template<class G, class ...> class B>
std::pair<typename treedec::gen_search::generic_elimination_search_base<A, O, B>::adjacency_iterator,
          typename treedec::gen_search::generic_elimination_search_base<A, O, B>::adjacency_iterator>
adjacent_vertices(
          typename treedec::gen_search::generic_elimination_search_base<A, O, B>::vertex_descriptor v,
			 treedec::gen_search::generic_elimination_search_base<A, O, B> const& o)
{
	return o.adjacent_vertices(v);
}

} // boost
#endif
