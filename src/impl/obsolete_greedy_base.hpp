#ifndef TREEDEC_OBSOLETE_GREEDY_BASE_HPP
#define TREEDEC_OBSOLETE_GREEDY_BASE_HPP

namespace treedec{

namespace impl{

// use greedy_base instead.
template<class G_t, template<class G, class...> class CFGT_t=algo::default_config>
class greedy_heuristic_base : public ::treedec::algo::draft::algo1{
public:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename treedec::graph_traits<G_t>::treedec_type T_t;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename std::vector<vertex_descriptor> bag_t;
    typedef typename std::vector<vertex_descriptor> O_t;

    greedy_heuristic_base(G_t &G, unsigned ub, bool ignore_isolated_vertices=false)
      : algo1("."), _g(G), _t(NULL), _own_o(),
        _ub_in(ub), _iiv(ignore_isolated_vertices), _i(0),
        _min(0), _ub(0), _current_N(&bag_i), _num_vert(boost::num_vertices(_g))
    {
        _o = new O_t;
        _o->resize(_num_vert);
        _do_tree_decomposition=true; // for now.
    }


    ~greedy_heuristic_base(){
        if(_own_o){
            delete _o;
        }
    }

    O_t& get_elimination_ordering() {

        return *_o;
    }

#if 0
    template<class O>
    void get_elimination_ordering(O& o) const { untested();
        O = *_o; // doesnt work like this.
    }
#endif

    template<class T>
    void get_tree_decomposition(T& t)const{
        // std::cerr << "hmm " << _o->size() << " " << _bags.size() << "\n";
        assert(_o->size()<=_bags.size()); // this is obsolete anyway
        assert(_o->size()==_num_vert); // this is relevant.

        // yuck... will be obsolete with FI rework
        typename std::vector<
            std::pair<vertex_descriptor, bag_t>
                > bags(_num_vert);
        typename std::vector<unsigned> io(_num_vert);

        // stuff center and friends into "skeleton"
        // _num_vert can be less than order/bags size
        for(unsigned i = 0; i < _num_vert; i++){
            bags[i].first = (*_o)[i];
            bags[i].second = _bags[i];
            // io[ (*_o)[i] ] = i;
        }

        treedec::detail::skeleton_to_treedec(_g, t, bags, *_o, _i);
    }

    vertices_size_type get_bagsize(){
        return _ub+1;
    }

    virtual void initialize() = 0;
    virtual void next(vertex_descriptor &c) = 0;
    virtual void eliminate(vertex_descriptor v) = 0;
    virtual void postprocessing() = 0;

    void disable_td(){
        // stupid hack
        _do_tree_decomposition = true;
    }

    void do_it() {
        if(_do_tree_decomposition){
            _t=new T_t;
            // bags seem to be unnecessary
            _bags.resize(_num_vert);
        }else{untested();
        }


        timer_on();

        if(!_num_vert){
            timer_off();
            return;
        }else{
        }

        assert(_o);
        O_t& elim_vertices = *_o;

#ifndef NDEBUG
        check(_g);
#endif

        initialize();

        _o->resize(_num_vert);

        assert(elim_vertices.size() == _num_vert);

        while(boost::num_edges(_g) > 0){
            vertex_descriptor c;

            next(c);

            //Abort if the width of this decomposition would be larger than 'ub'.
            if(_min >= _ub_in){ untested();
                assert(_t); // ouch?
                _t->clear(); //could be also not the case
                throw exception_unsuccessful();
            }

            elim_vertices[_i] = c;

            if(_t){
                _current_N = &_bags[_i];
            }else{
            }

            _ub = (boost::out_degree(c, _g)>_ub)?boost::out_degree(c, _g):_ub;

            // assert(bags_i);?!

            eliminate(c);

            if(!_t){
                _current_N->clear();
            }else{
            }

            ++_i;
        }

        postprocessing();

//        assert(_i == num_vert); //or _i-1==num_vert?!

        timer_off();

    } // do_it


protected:
    G_t &_g;
    T_t* _t;
    O_t* _o;
    bool _own_o;

    vertices_size_type _ub_in;
    bool _iiv;
    size_t _i;
    unsigned _min;

    std::vector<bag_t> _bags;

    vertices_size_type _ub;

    bag_t bag_i;
    bag_t* _current_N;

    unsigned _num_vert;
private:
    bool _do_tree_decomposition;
};

}

}

#endif
