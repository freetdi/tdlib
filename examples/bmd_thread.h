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

template<class G, template<class H, class ... > class cfgt=treedec::algo::default_config>
class BMD_THREAD : public TWTHREAD<G, cfgt> {
public:
    typedef std::vector<int> iorder_t;
	typedef TWTHREAD<G, cfgt> base;
	typedef cfgt<G> CFG;
    BMD_THREAD( G& g, const std::string& name="BMD", mag_t m=M64)
        : base(g, name, 0), _g(g), _mag(m)
    {
        base::go(); // ??
    }

    void do_print_results(std::ostream& o)
    {
#ifdef HAVE_GALA_GRAPH_H
        _g.make_symmetric(true);
#endif
        base::print_results_order(o, _elimord);
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
        }else{ untested();
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
            treedec::impl::bmdo<sg_dvv16> A(g16, _elimord);
				A.do_it();
				s = A.bagsize();
        }else{ untested();
            trace1("BMD", boost::num_edges(g32));
            treedec::impl::bmdo<sg_dvv32> A(g32, _elimord);
				A.do_it();
				s = A.bagsize();
        }
        trace1("BMD", s);
        base::commit_result(s);
        base::unlock_results();
    }
private:
   G& _g; // print hack.
   iorder_t _elimord;
   mag_t _mag;
};
