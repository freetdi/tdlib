
// temporary.
#ifndef untested
#define untested()
#endif
#ifndef itested
#define itested()
#endif
#ifndef incomplete
#define incomplete()
#endif


namespace misc {

template<class G>
class DEGS{
    DEGS(const DEGS&){}

public:
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
    typedef std::set<vertex_descriptor> bag_type;
    typedef typename bag_type::iterator bag_iterator;
    typedef std::vector<std::set<vertex_descriptor> > container_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    DEGS(const G& g): _degs(boost::num_vertices(g)), _g(g)
    {
        vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
            _degs[boost::degree(*vIt, g)].insert(*vIt);
        }
    }

    size_t num_nodes() const{ untested();
        unsigned N=0;
        for(const_iterator i=_degs.begin(); i!=_degs.end(); ++i) { itested();
            N+=i->size();
        }
        return N;
    }

    void check()
    { // sometimes required when debugging fancy callbacks :/
#ifdef EXCESSIVE_DEG_DEBUG
            DEGS degs(_g);
            assert(_degs.size()==degs.size());
            assert(size()==boost::num_vertices(_g));
            assert(num_nodes()==degs.num_nodes());

            iterator j=degs._degs.begin();
            unsigned N=0;
            for(iterator i=_degs.begin(); i!=_degs.end();) {
                assert(N<boost::num_vertices(_g));
                unsigned I=i->size();
                unsigned J=j->size(); //actual _g

                if(I>J){
                    std::cerr<<"mismatch " << I << " " << J << "\n";
                    std::cerr<<"extra node " << *i->begin() << " of deg " << N << " in degs\n";
                }else if(I<J){
                    std::cerr<<"mismatch " << I << " " << J << " in " << N << "\n";
                    std::cerr<<"extra node " << *j->begin() << " of deg " << N << " in g\n";
                }
                assert(I==J);
                ++i;
                ++j;
                ++N;
            }
            assert(N==boost::num_vertices(_g));
#endif
    } //void check()

    std::set<vertex_descriptor>& operator[](size_t x){return _degs[x];}

    size_t size() const {return _degs.size();}

//private: // later.
    container_type _degs;
private:
    const G& _g;
}; // DEGS

template<class G>
struct deg_chooser{
    typedef typename misc::DEGS<G> degs_type;
    typedef typename misc::DEGS<G> type;
    static void alloc_init(size_t){
    }
};

template<class VC, class G, class CB>
void make_clique(VC V, G& g, CB* cb)
{
    //FIXME: permit any container...
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

    typename boost::graph_traits<G>::adjacency_iterator nIt1, nIt2, nEnd;
    typename misc::DEGS<G> &degs=*cb->_degs;

    std::set<vertex_descriptor> redeg;
    for(boost::tie(nIt1, nEnd) = V; nIt1 != nEnd; nIt1++){
        if(cb){
            unsigned deg = boost::degree(*nIt1,g);
            size_t n=degs[deg].erase(*nIt1);
            (void)n;
            assert(n==1);
        }
    }
    for(boost::tie(nIt1, nEnd) = V; nIt1 != nEnd; nIt1++){
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){
            bool newedge = boost::add_edge(*nIt1, *nIt2, g).second;
                assert(boost::degree(*nIt1,g));
                assert(boost::degree(*nIt2,g));

            if(cb){
                // need to redegree all vertices.
            }
            else if(newedge){ untested();
                // need to redegree *nIt1...
                if(redeg.insert(*nIt1).second){ untested();
                    unsigned deg = boost::degree(*nIt1,g);
                    // std::cerr << "queued 4 redeg " << *nIt1 << " deg" << deg-1 << "\n";
                    size_t n=degs[deg-1].erase(*nIt1);
                    (void)n;
                    assert(n==1);
                }else{ untested();
                    // already queued for redegree
                }

                // need to redegree *nIt2...
                if(redeg.insert(*nIt2).second){ untested();
                    unsigned deg = boost::degree(*nIt2,g);
                    // std::cerr << "queued 4 redeg " << *nIt2 << " deg" << deg-1 << "\n";
                    size_t n=degs[deg-1].erase(*nIt2);
                    (void)n;
                    assert(n==1);
                }else{ untested();
                    // already queued for redegree
                }
            }else{ untested();
                // nothing to do.
            }
        }
    }
    if(cb){
        for(boost::tie(nIt1, nEnd) = V; nIt1 != nEnd; nIt1++){
            (*cb)(*nIt1);
        }
    }else{ incomplete();
        for(typename std::set<vertex_descriptor>::iterator i=redeg.begin(); i!=redeg.end(); ++i){
            // std::cerr << "redegging " << *i << "\n";
            if(cb){ untested();
            }else{ incomplete();
                // just reinsert
                size_t deg=boost::degree(*i,g);
                assert(deg>0);
                bool done=degs[deg].insert(*i).second;
                (void)done;
                assert(done);
            }
        }
    }

    if(cb){
        // degs might be in an inconsistent state, depending on what cb is
    }else{ untested();
        degs.check();
    }
}

} //namespace MISC

