// Felix Salfelder, 2016
//
// (c) 2016 Goethe-Universität Frankfurt
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//  boost min-degree thread
//
// BUG: use these globally
typedef boost::property<treedec::bag_t, std::vector<uint16_t> > uint16_bag_p;
typedef boost::property<treedec::bag_t, std::vector<uint32_t> > uint32_bag_p;
typedef boost::property<treedec::bag_t, std::vector<uint64_t> > uint64_bag_p;

// #define tree_directedness_ boost::bidirectionalS //  broken?
#define tree_directedness_ boost::undirectedS

// BUG: not here
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint16_bag_p> _gsgvvu16_treedec;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint32_bag_p> _gsgvvu32_treedec;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint64_bag_p> _gsgvvu64_treedec;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint64_bag_p> _gsgvvu64_treedec;

template<class G, template<class H, class ... > class cfgt=treedec::algo::default_config>
class BMD_THREAD : public TWTHREAD<G, cfgt> {
public:
	using TWTHREAD<G, cfgt>::get_result;
public:
    typedef std::vector<int> iorder_t;
	typedef TWTHREAD<G, cfgt> base;
	typedef cfgt<G> CFG;
    BMD_THREAD( G& g, const std::string& name="BMD", mag_t m=M64)
        : base(g, name, 0), _g(g), _mag(m)
    {
        base::go(); // ??
    }

    void do_print_results(std::ostream& o) {
#ifdef HAVE_GALA_GRAPH_H
        // _g.make_symmetric(true); // ???
#endif
        // base::print_results_order(o, _elimord);

        std::cerr<< "c size " << boost::num_vertices(_g) << "\n";
        // auto &g=TWTHREAD<G>::_g;
		  treedec::grtdprinter<G> P(o, _g);

#if 0
		  size_t numbags = boost::num_vertices(_t);
		  P.head(numbags, get_result());
		  boost::copy_graph(_t, P);
#else
		  // make -q work. possibly too slow
		  if(_a16){
			  _gsgvvu16_treedec t;
			  _a16->get_tree_decomposition(t);
//			  _a16->get_tree_decomposition(P); // TODO
			  size_t numbags = boost::num_vertices(t);
			  P.head(numbags, get_result());
			  boost::copy_graph(t, P);
		  }else if(_a32){ untested();
			  _gsgvvu32_treedec t;
			  _a32->get_tree_decomposition(t);
//			  _a16->get_tree_decomposition(P); // TODO
			  size_t numbags = boost::num_vertices(t);
			  P.head(numbags, get_result());
			  boost::copy_graph(t, P);
		  }else{
			  unreachable();
		  }
#endif
    }

    void run() {
#ifdef HAVE_GALA_GRAPH_H
        // TODO: faster with "remove-only" stuffed vector-graph...
        sg_dvv16* pg16;
        sg_dvv32* pg32;
        if(_mag<M16){
            pg16=new sg_dvv16(base::_g);
            treedec::check(*pg16);
            pg16->make_symmetric(true);
            assert_symmetric(*pg16);
        }else{itested();
            pg32=new sg_dvv32(base::_g);
            treedec::check(*pg32);
            pg32->make_symmetric(true);
            assert_symmetric(*pg32);
        }
        auto& g16=*pg16;
        auto& g32=*pg32;
        //auto ne=boost::num_edges(g);
        //assert(ne*2==boost::num_edges(g));
#else
        balvvd_t g16;
        boost::copy_graph(base::_g, g16);
        treedec::make_symmetric(g16);
        balvvd_t& g32(g16);
#endif
        unsigned s;
        if(_mag<M16){
            trace1("BMD", boost::num_edges(g16));
            _a16 = new treedec::impl::bmdo<sg_dvv16>(g16);
				_a16->do_it();
				s = _a16->bagsize();
        }else{itested();
            trace1("BMD", boost::num_edges(g32));
            _a32 = new treedec::impl::bmdo<sg_dvv32>(g32);
				_a32->do_it();
				s = _a32->bagsize();
        }
        trace1("BMD", s);
        base::commit_result(s);
        base::unlock_results();
    }
private:
   G& _g; // print hack.
//   iorder_t _elimord;
   decomp_t<G> _t;
   mag_t _mag;
	treedec::impl::bmdo<sg_dvv16>* _a16{nullptr};
	treedec::impl::bmdo<sg_dvv32>* _a32{nullptr};
};
