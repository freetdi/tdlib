// (c) 2016, 2017 Felix Salfelder
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
// compute gr to td using various algorithms (in parallel)
//

#include "config.h"

#ifdef DEBUG
#undef NDEBUG
#endif

#include <algorithm>
#include <assert.h>
#include <atomic>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <random>
#include <setjmp.h>
#include <signal.h>
#include <stdlib.h>


#define assert_symmetric(g) { \
    unsigned i=0; \
    auto E=boost::edges(g); \
    for(;E.first!=E.second; ++E.first){ \
        ++i;\
        assert( boost::edge(boost::source(*E.first, g), \
                    boost::target(*E.first, g),g).second); \
        assert( boost::edge(boost::target(*E.first, g), \
                    boost::source(*E.first, g),g).second); \
    }\
    trace1("symmetric", i); \
}

/// old...
enum thread_n{
    nMD = 0,
    nFI = 1,
    nMSVS = 2,
    nBMD = 3,
    nSOME = 4,
    nEX = 5,
    nEX17 = 6,
    nP17 = 7,
    nTH = 8,
    nPPFITM = 9,
    nFITM = 10,
    nPPMD = 11,
    nPPFI = 12,
    nTOTAL = 13
};

std::mutex best_mutex;
std::condition_variable cv;
static std::atomic<unsigned> global_result;
// static std::atomic<unsigned> threads_running;
bool trace=false;

#include <boost/graph/graph_traits.hpp>
#include <tdlib/graph_traits.hpp>

//
#ifdef HAVE_GALA_GRAPH_H
#include <gala/boost.h>
#include <tdlib/directed_view.hpp>
#ifdef USE_RANDOM_MD
#include <gala/examples/ssg_random.h>
#endif
// #include <boost/graph/minimum_degree_ordering.hpp>
#include <gala/examples/ssg32i.h>
#include <gala/examples/ssg16i.h>
#include <gala/examples/ssg16ia.h>
#include <gala/examples/svbs.h>
#include <gala/examples/svbs_random.h>
#include <gala/immutable.h>
#else
#warning no gala. does not fully work
#endif

#include "timer.h"

#include <tdlib/printer.hpp>
#include <tdlib/combinations.hpp>
#include <tdlib/elimination_orderings.hpp>
#include <boost/graph/copy.hpp>

// undirected simple loopless graph
template<class G>
struct uvv_config : gala::graph_cfg_default<G> {
    static constexpr bool is_directed=false;
    static constexpr bool force_simple=true;
    static constexpr bool force_loopless=true; // not used yet?
    // static constexpr bool force_symmetric=true; // meaninngless (undirected)
    // typedef tdDEGS<G> degs_type; // obsolete.
};
typedef gala::graph<std::vector, std::vector, uint16_t, uvv_config> sg_dvv16;
typedef gala::graph<std::vector, std::vector, uint32_t, uvv_config> sg_dvv32;

struct test{
    test(){
        sg_dvv16 a;
   treedec::draft::printer<sg_dvv16> x(std::cout, a);
   boost::add_vertex(x);
    }
} x;

#ifdef HAVE_GALA_GRAPH_H
#include <gala/boost_copy.h>
#include <gala/td.h>
#endif
#include "twthread.hpp"

using treedec::TWTHREAD;
using treedec::draft::TWTHREAD_BASE;

unsigned TWTHREAD_BASE::_running;

#include "grparse.h"

// bug, still tdlib

////// why is this necessary? //////
using boost::out_edges;
using boost::out_degree;
using boost::degree;
using boost::source;
using boost::target;
////////////////////////////////////

#ifdef USE_GALA
BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<ssg_16i> ));
typedef typename treedec::graph_traits<ssg_16i>::immutable_type check_type;
// boost::iterator_traits<check_type::out_edge_iterator>::value_type A;
std::iterator_traits<check_type::out_edge_iterator>::value_type A;
BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<check_type> ));
typedef typename boost::graph_traits<ssg_16i>::vertex_descriptor tttt;
#endif

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> balu_t;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> bald_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> balvvd_t;


template<class G>
using decomp_t = typename treedec::graph_traits<G>::treedec_type;

enum mag_t{
    Munknown=0,
    M6=6, // <= 64 verts
    M7=7, // <= 128 verts
    M8=8,
    M16=16,
    M15=15,
    M31=31,
    M32=32,
    M64=64
};

// fixme: move to other compilation unit
#include <stdarg.h>
#include <stdio.h>
int errorlevel=bLOG;
template<class X=balu_t, class ... rest>
struct grtd_algo_config : treedec::algo::default_config<X, rest...>{
    static void message(int badness, const char* fmt, ...) { untested();
        trace3("message", badness, fmt, errorlevel);
        if (badness >= errorlevel){ untested();
            char buffer[2048] = "c ";
            va_list arg_ptr;
            va_start(arg_ptr,fmt);
            vsprintf(buffer+2,fmt,arg_ptr);
            va_end(arg_ptr);
            std::cout << buffer;
        }else{ untested();
        }
    }
};
/*--------------------------------------------------------------------------*/


#ifdef USE_GALA
#define FIgraph_t ssg_16i
#define MDgraph_t ssg16_random
#else
#define FIgraph_t balu_t
#define MDgraph_t balu_t
#endif

#ifdef USE_GALA
#include "gala_graphs.h"

#endif

#ifdef USE_GALA
#ifdef USE_BMD
#include "bmd_thread.h"
#endif
#if defined USE_SOME
#include "some_thread.h"
#endif
#endif
#ifdef USE_THORUP
#include "th_thread.h"
#endif

struct sigaction sa;

sig_atomic_t int_received;
sig_atomic_t term_received;
bool quiet=false;
volatile unsigned finished;

static void nohandler(int, siginfo_t *, void *)
{
}
static void term_handler(int signum, siginfo_t *, void *)
{
    if(signum==SIGTERM){
        fprintf(stderr, "TERM handler\n");
        sigset_t mask;
        sigemptyset (&mask);
        sigaddset (&mask, signum);
        if (sigprocmask(SIG_BLOCK, &mask, NULL)!=0) { untested();
            printf("cannot reinstall TERM handler");
            exit(3);
        }else{ untested();
        }
        term_received = 1;
    }else{
        unreachable();
    }
}
static void int_handler(int signum, siginfo_t *, void *)
{
    if(signum==SIGUSR1 || signum==SIGINT){
        fprintf(stderr, "USR1 handler\n");
        int_received = 1;
    }else{ untested();
    }
}

#if 0 // later
static void fin(volatile unsigned & finished)
{
    std::unique_lock<std::mutex> scoped_lock(best_mutex);
    /// mutex running?!
    if(trace) std::cerr << "countdown" << finished << " " << TWTHREAD_BASE::_running << "\n";
    ++finished;
    if(finished==TWTHREAD_BASE::_running){
        if(trace) std::cerr << "signalling end\n";
//        scoped_lock.wait(cv);
        scoped_lock.unlock();
        kill(getpid(), SIGINT); // pass by "wait" in mainloop??
        cv.notify_all(); // needed after interrupt. main thread might be waiting for us.
    }
}
#endif

#ifdef USE_RANDOM_MD
#include "random_md.h"
#endif

#ifdef USE_GALA
#ifdef USE_EX17
#include "ex17_thread.h"
#endif
#endif

#ifdef USE_GALA
#ifdef USE_FI
#include "fi_thread.h"
#endif

#ifdef USE_EX
#include "ex_thread.h"
#endif
#endif

#ifdef USE_MSVS_TRIVIAL
#include "msvs_thread.h"
#endif

#ifdef USE_P17
#include "he17_thread.h"
#endif


static bool sig_received()
{
    return term_received || int_received;
}

template<class T>
void join_and_cleanup(T* t)
{
    if(t){ untested();
        t->join();
        delete t;
        t=NULL;
    }else{ untested();
    }
}
template<class T>
void shutdown_gently(T* t)
{
    if(t){ untested();
        t->interrupt();
        t->join();
    }else{ untested();
    }
}


template<class T>
void cleanup(T const& t)
{
#ifndef GENTLE_SHUTDOWN
    (void)t;
    if(trace){
        std::cerr << "forced shutdown\n";
    }
    {
        sa.sa_sigaction = nohandler;
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = SA_RESTART | SA_SIGINFO;
        if (sigaction(SIGSEGV, &sa, NULL) == -1){
            exit(3);
        }
    }
#else // requires interruption points
// does not really work, as threads do not come back...
    if(trace){
        std::cerr << "slowly shutting down.";
    }
    for(auto ti : t){
        std::cerr << ".";
        shutdown_gently(ti);
    }
    if(trace){
        std::cerr << "done..\n";
    }
#endif
}

void mainloop()
{
    sigset_t mask;
    sigemptyset (&mask);
    sigaddset (&mask, SIGTERM);

    while(true){
        if(term_received){
            if(trace) std::cerr << "term received\n";
            if (sigprocmask(SIG_BLOCK, &mask, NULL)<0) { untested();
                std::cerr << "error blocking SIGTERM\n";
                exit(99);
            }
            return;
        }else if(int_received){
            std::lock_guard<std::mutex> scoped_lock(best_mutex);
            std::cout << "c " << global_result << "\n";
            // if(! threads_running){
            //     break;
            // }else{
            // }
        }else{
        }
        int_received=0;
        if(!TWTHREAD_BASE::_running){
            // all threads masked?
            // all threads completed (that was quick!)
            break;
        }else if(!sig_received()){ untested();
            trace0("pause");
            pause();
        }
        // get here after the sighandler has been executed
    }
}

template<class V>
unsigned find_best(V const& v)
{
    unsigned best=-1;
    unsigned bestindex=nTOTAL;

    for(unsigned i=0; i<v.size(); ++i){
        if(!v[i]){
            continue;
        }
        auto result=v[i]->get_result();
        if (result < best){
            bestindex = i;
            best = result;
        }
    }
    return bestindex;
}

template<class TS, class T>
void reg_thread(TS& threads, unsigned number, T ptr){
    threads[number] = ptr;
    grtd_algo_config<balu_t>::message(bLOG, "thread %s\n", ptr->name().c_str());
//    ++threads_running;
}

// wrap graph context in edge iterator
template<class G>
class raw_edges{
public: // types
    typedef typename boost::graph_traits<G>::edge_iterator edge_iterator;
    class edge_iter{
    public:
        edge_iter(G const& g) : _g(g){
        }
    public:
        edge_iter& operator=(edge_iterator const& i){
            _i = i;
            return *this;
        }
        bool operator!=(edge_iter const& o){
            return _i!=o._i;
        }
        bool operator==(edge_iter const& o){
            return _i==o._i;
        }
        edge_iter& operator++(){
            ++_i;
            return *this;
        }
        std::pair<size_t, size_t> operator*(){
            auto s=boost::source(*_i, _g);
            auto t=boost::target(*_i, _g);
            return std::make_pair(s, t);
        }
    private:
        edge_iterator _i;
        G const& _g;
    };
public:
    raw_edges(G& g) : first(g), second(g) {
        auto p=boost::edges(g);
        first = p.first;
        second = p.second;
    }

public:
    edge_iter first;
    edge_iter second;
};


template<class P>
void twh(P& p, mag_t m, unsigned mask)
{
    std::vector<TWTHREAD_BASE*> threads(nTOTAL);
    TWTHREAD_BASE::_running=1;

#ifdef USE_GALA
//    typedef gala::graph<std::vector, std::vector, uint32_t, uvv_config> sg_dvv32;

    typedef sg_dvv16 uG16;
    typedef sg_dvv32 uG32;

    typedef sg_dpvv16 G16p;
    typedef sg_dpvv32 G32p;

    size_t n=boost::num_vertices(p);
    size_t e=boost::num_edges(p);

    uG16 g16;
    uG32 g32;
    sg_odsvv16 tg;
    auto BE=raw_edges<P>(p);
    auto B=BE.first;
    auto E=BE.second;

    if(m>M15){ untested();
        G32p pg32(B, E, n, e);
        g32 = std::move(pg32);

        assert(boost::num_edges(g32)==e); // usually.
    }else{
        G16p pg16(B, E, n, e);
//		boost::print_graph(pg16);
        assert(boost::num_edges(pg16)==e);
#ifdef USE_THORUP
        // thorup only seems to work on half graphs (oriented directed)
        if((mask & ( 1 << nTH ))) {
            boost::copy_graph(pg16, tg);
            reg_thread(threads, nTH, new TH_THREAD<sg_odsvv16, grtd_algo_config>(tg, "TH16"));
        }else{
        }
#endif
        g16 = std::move(pg16);
//        assert(boost::num_edges(g16)==e); // not if there were multiedges
    }

/*--------------------------------------------------------------------------*/
#else
    incomplete();
    typedef balu_t G;
    typedef balu_t uG16;
    typedef balu_t uG32;
    G g(p->begin(), p->end(), p->num_vertices(), p->num_edges());
    G& g16(g);
    G& g32(g);
#endif
/*--------------------------------------------------------------------------*/


    std::cout << "c n: " << n << ", e: " << e << std::endl;
#ifdef USE_GALA
    std::cout << "c gala on" << std::endl;
#endif


    if(trace){
        std::cerr << "starting threads for " << m << " bit mode\n";
    }

/*--------------------------------------------------------------------------*/
#ifdef USE_THORUP_BROKEN
    if(!(mask & ( 1 << nTH ))) {
    }else if( m < M16){ untested();
            // g16.hacksort();
        reg_thread(threads, nTH, new TH_THREAD<uG16, grtd_algo_config>(g16, "TH16"));
    }else{
    }
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_FIPPTM
    if(!(mask & ( 1 << nPPFITM ))) {
    }else if( m < M16){ untested();
        reg_thread(threads, nPPFITM, new PPFITM_THREAD<uG16, grtd_algo_config>(g16, "FIPPTM_16"));
    }else{ untested();
        reg_thread(threads, nPPFITM, new PPFITM_THREAD<uG32, grtd_algo_config>(g32, "FIPPTM_32"));
    }
/*--------------------------------------------------------------------------*/
#if 1
    if(!(mask & ( 1 << nFITM ))) {
    }else if( m < M16){ untested();
        reg_thread(threads, nFITM, new FITM_THREAD<uG16, grtd_algo_config>(g16, "FITM_16"));
    }else{ untested();
        reg_thread(threads, nFITM, new FITM_THREAD<uG32, grtd_algo_config>(g32, "FITM_32"));
    }
/*--------------------------------------------------------------------------*/
    if(!(mask & ( 1 << nPPFI ))) {
    }else if( m < M16){ untested();
        reg_thread(threads, nPPFI, new PPFI_THREAD<uG16, grtd_algo_config>(g16, "PPFI_16"));
    }else{ untested();
        reg_thread(threads, nPPFI, new PPFI_THREAD<uG32, grtd_algo_config>(g32, "PPFI_32"));
    }
/*--------------------------------------------------------------------------*/
    if(!(mask & ( 1 << nPPMD ))) {
    }else if( m < M16){ untested();
        reg_thread(threads, nPPMD, new PPMD_THREAD<uG16, grtd_algo_config>(g16, "PPMD_16"));
    }else{ untested();
        reg_thread(threads, nPPMD, new PPMD_THREAD<uG32, grtd_algo_config>(g32, "PPMD_32"));
    }
#endif
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_RANDOM_MD
    if(!(mask & ( 1 << nMD ))) {
    }else if( m <= M16){ untested();
        reg_thread(threads, nMD, new MD_THREAD<uG16>(g16, "MD16"));
    }else{ untested();
        reg_thread(threads, nMD, new MD_THREAD<uG32>(g32, "MD32"));
    }
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_MSVS_TRIVIAL
    if(!(mask & ( 1 << nMSVS ))){
    }else if( m <= M16){
        reg_thread(threads, nMSVS, new MSVS_THREAD<ssg_16i>(g16, "MSVS16"));
    }else{
        reg_thread(threads, nMSVS, new MSVS_THREAD<ssg_32i>(g32, "MSVS32"));
    }
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_FI
    if(! ( mask & ( 1 << nFI ))) { untested();
    }else if( m < M16){
        reg_thread(threads, nFI, new FI_THREAD<uG16, grtd_algo_config>(g16, "FI16"));
    }else if( m < M32){ untested();
        reg_thread(threads, nFI, new FI_THREAD<uG32, grtd_algo_config>(g32, "FI32"));
    }else{ untested();
        incomplete();
//        threads[nFI] = new FI_THREAD<uG32>(g32, "FI32");
    }
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_P17
    if(! ( mask & ( 1 << nP17 ))) {
    }else if( m < M16){ untested();
        // need "less than M16", due to bucket sorter quirk
        reg_thread(threads, nP17, new P17_THREAD<uG16, grtd_algo_config>(g16, "P17_16"));
    }else{ untested();
        reg_thread(threads, nP17, new P17_THREAD<uG32, grtd_algo_config>(g32, "P17_32"));
    }
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_SOME // 16 bit
    if(!(mask & ( 1 << nSOME ))) { untested();
    }else if( m <= M16){
        threads[nSOME] = new SOME_THREAD<uG16>(g16, "SOME");
    }else{
        incomplete();
    }
#endif
/*--------------------------------------------------------------------------*/
#ifdef USE_BMD
    if(! ( mask & ( 1 << nBMD ))) {
    }else if(m < M16){
        reg_thread(threads, nBMD, new BMD_THREAD<uG16, grtd_algo_config>(g16, "BMD16", m));
    }else if(m < M32){ untested();
        reg_thread(threads, nBMD, new BMD_THREAD<uG32, grtd_algo_config>(g32, "BMD32", m));
    }else{
        incomplete();
        // need 64 bit BMD...
        // perhaps gala<vec,vec,void*,directed>...?
        // (distribute?)
//        threads[nBMD] = new BMD_THREAD<G64>(g32, results[nBMD]);
    }
#endif
/*--------------------------------------------------------------------------*/
#if defined(USE_EX17) && defined(USE_GALA)
    if(! ( mask & ( 1 << nEX17 ))) {
    }else if(m > M15){ untested();
        // does this even make sense?
        // maybe for very sparse graphs...
        threads[nEX17] = new EX17_THREAD<uG32, grtd_algo_config>(g32, "EX17_32");
    }else{
        threads[nEX17] = new EX17_THREAD<uG16, grtd_algo_config>(g16, "EX17_16");
    }
#endif
/*--------------------------------------------------------------------------*/
#if defined(USE_EX) && defined(USE_GALA)
    if(! ( mask & ( 1 << nEX ))) {
    }else if(m > M16){ untested();
        // does this even make sense?
        // maybe for very sparse graphs...
        threads[nEX] = new EX_THREAD<uG32, grtd_algo_config>(g32, "EX32");
    }else{
        threads[nEX] = new EX_THREAD<uG16, grtd_algo_config>(g16, "EX16");
    }
#endif
/*--------------------------------------------------------------------------*/

    --TWTHREAD_BASE::_running;

//     while(running!=finished){ untested();
//         std::unique_lock<std::mutex> lk(best_mutex);
//         cv.wait(lk);
//     }


    mainloop();
    trace0("exited mainloop");

    if(trace){
        unsigned x=1;
        for(unsigned i=0; i<nTOTAL; ++i){
            if(! (x&mask) ){
            }else if(!threads[i]){
            }else
                //if(threads[i]->get_result()!=-1u)
               {
                std::cout << "c " << threads[i]->name() << ": "
                    << threads[i]->get_result() << "\n";
           // }else{
            }
            x*=2;
        }
    }

    trace0("acquiring lock");
    std::lock_guard<std::mutex> scoped_lock(best_mutex); // needed?

    unsigned best_tid=find_best(threads);

    if(best_tid<threads.size()){
    grtd_algo_config<balu_t>::message(bDEBUG, "best: %s\n",
            threads[best_tid]->name().c_str());
    // auto const& best_thread = threads::best();
    auto& best_thread = *threads[best_tid];

    if(quiet){ untested();
    }else{

        // FIXME: no switch here.
        switch(best_tid){
        case nSOME:
        case nBMD:
#ifdef USE_GALA // tmp hack
            // g.make_symmetric(true);
#endif
        default:
            best_thread.interrupt(); // uuh, wait here until copy is finished?!
//            grtd_algo_config<balu_t>::message(bDEBUG, "printing result\n");
            best_thread.print_results(std::cout);
            break;
        case nTOTAL:
            std::cout << "no result\n";
        } // switch

        std::cout << "c done\n";

    }
    }
    // THREADS..::cleanup();
    cleanup(threads);
}

static void setup_handlers()
{
    sa.sa_sigaction = term_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART | SA_SIGINFO | SA_RESETHAND;
    if (0
//     || sigaction(SIGINT,  &sa, NULL) == -1
     || sigaction(SIGTERM, &sa, NULL) == -1
     ) {
        std::cerr << "error installing sighandler\n";
        exit(1);
    }
    sa.sa_sigaction = int_handler;
    sa.sa_flags = SA_RESTART | SA_SIGINFO;
    if (sigaction(SIGUSR1, &sa, NULL) == -1) {
        std::cerr << "error installing sighandler\n";
        exit(1);
    }
    if (sigaction(SIGINT, &sa, NULL) == -1) {
        std::cerr << "error installing sighandler\n";
        exit(1);
    }
}

mag_t m=Munknown;
unsigned mask=-1;
unsigned mask_in=0;

static void parseargs(int argc, char * const * argv)
{

    int i=1;
    while(i<argc){
        if(!strncmp("-q", argv[i], 2)){ untested();
            quiet=true;
        }else if(!strncmp("-L", argv[i], 2)){ untested();
            errorlevel=bLOG;
        }else if(!strncmp("-D", argv[i], 2)){ untested();
            errorlevel=bDEBUG;
        }else if(!strncmp("-N", argv[i], 2)){ untested();
            errorlevel=bNOERROR;
        }else if(!strncmp("--dot", argv[i], 5)){ untested();
            fformat = f_DOT;
        }else if(!strncmp("--he17", argv[i], 6)){ untested();
            mask_in |= (1<<nP17);
        }else if(!strncmp("--ex17", argv[i], 6)){ untested();
            mask_in |= (1<<nEX17);
        }else if(!strncmp("--thorup", argv[i], 8)){ untested();
            mask_in |= (1<<nTH);
        }else if(!strncmp("--ppfitm", argv[i], 8)){ untested();
            mask_in |= (1<<nPPFITM);
        }else if(!strncmp("--ppfi", argv[i], 6)){ untested();
            mask_in |= (1<<nPPFI);
        }else if(!strncmp("--fi", argv[i], 4)){ untested();
            mask_in |= (1<<nFI);
        }else if(!strncmp("--fitm", argv[i], 6)){ untested();
            mask_in |= (1<<nFITM);
        }else if(!strncmp("--ppmd", argv[i], 6)){ untested();
            mask_in |= (1<<nPPMD);
        }else if(!strncmp("-T", argv[i], 2)){ untested();
            trace = true;
            errorlevel=bTRACE;
        }else if(!strncmp("-t", argv[i], 2)){ untested();
            trace = true;
            errorlevel=bTRACE;
            std::cerr << "tracing ON\n";
        }else if(!strncmp("-w", argv[i], 2)){ untested();
            m = M32;
            std::cerr << "wide ON\n";
        }else if(!strncmp("-m", argv[i], 2)){
            ++i; // incomplete: range check?
            mask_in = atoi(argv[i]);
            trace2("mask arg", argv[i], mask_in);
        }else if(!strncmp("-s", argv[i], 2)){ untested();
            // not yet
            ++i;
        }else{ untested();
            std::cerr << "error parsing args " << argv[i] << "\n";
            exit(1);
        }
        ++i;
    }
    if(mask_in){
        mask = mask_in;
        trace1("mask set", mask);
    }else{
        untested();
    }
 }

	 template<class P>
 mag_t choose_m(mag_t m, P const&){

    if(m){ untested();
    }else if(global_result >= (1l<<31)-1){ untested();
        m = M32;
    }else if(global_result >= (1<<16)-1){ untested();
        m = M31;
    }else if(global_result >= (1<<15)-1){ untested();
        m = M16;
    }else{
        m = M15;
    }
    return m;
 }

// vim:ts=8:sw=4:et:
