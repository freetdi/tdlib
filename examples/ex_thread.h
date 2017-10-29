
template<class G>
class EX_THREAD : public TWTHREAD<ssg_16i> {
public:
    typedef TWTHREAD<ssg_16i> base;
#if 1 // slow preprocessing
    typedef ssg_16i G_work;
#else // broken preprocessing
typedef sg_dvv16 G_work;
#endif

    EX_THREAD( G g /* create local copy ... */, const std::string& name="EX")
        : base(g, name, 0), _work(g) // <- .. store it here
    {
        assert(boost::num_edges(g) == boost::num_edges(_work));
        treedec::check(g);
        treedec::check(_work);
        go();
    }

    void do_print_results(std::ostream& o)
    {
        trace2("EXR", base::_result, treedec::get_bagsize(_t));
        base::print_results_tree(o, _t, &_work);
    }

    void run()
    {
        auto &g=_work; // nonconst, local
        treedec::exact_decomposition_cutset(g, _t);
        unsigned r = treedec::get_bagsize(_t); // inefficient
        trace1("excut done", r);
        assert(boost::num_vertices(g) || r==0);
        base::commit_result(r);
        base::unlock_results();
        kill(getpid(), SIGTERM);
        base::just_wait();
    }
private:
//   result_t _result;
   G_work _work;
   decomp_t<G_work> _t;
};
