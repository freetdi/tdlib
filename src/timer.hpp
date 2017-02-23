/*
 * Copyright (C) 2017 Felix Salfelder
 * Author: Felix Salfelder
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * CPU time accounting, inspired by gnucap/l_timer.h
 */
#ifndef TD_TIMER_H
#define TD_TIMER_H
#include <iostream>
#include <time.h>
#include "trace.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif
/*--------------------------------------------------------------------------*/
// struct timeval
// {
//   __time_t tv_sec;    /* Seconds.  */
//   __suseconds_t tv_usec; /* Microseconds.  */
// };
/*--------------------------------------------------------------------------*/
inline double to_double(const timespec t){
  return t.tv_sec + 1e-9*t.tv_nsec;
}
/*--------------------------------------------------------------------------*/
inline double to_double(const timeval t){
  return t.tv_sec + 1e-6*t.tv_usec;
}
/*--------------------------------------------------------------------------*/
class TIMER_BASE {
  public:
    TIMER_BASE(const std::string& label):_name(label){}
    virtual ~TIMER_BASE(){}
    virtual TIMER_BASE* clone() const = 0;
    virtual void fullreset() {untested(); return;}
    virtual void reset() {untested(); return;}
    virtual void start() = 0; // {untested(); return;}
    virtual void stop(){untested(); return;}
    virtual void check() {
      if (_running) {
	stop();
	start();
      }else{
	untested();
      }
    }
    virtual double	elapsed()const
    {itested();
      return 0.;
    }
    virtual double	total()const
    { untested();
      return 0;
    }
    virtual bool	is_running()const
    {untested();
      return false;
    }

    virtual std::ostream& print(std::ostream& s) const {
      s << _name.c_str();
      for(unsigned i=_name.size(); i<13; ++i){
	s << " ";
      }
      s <<  total() << "\n";
      return s;
    }

//    std::ostream& operator<<(std::ostream& s) {print(s); return s;}
    TIMER_BASE* set_name(const std::string& n){ _name=n; return this; }

  private:
    std::string _name;
  protected:
    bool _running;
};
inline std::ostream& operator<<(std::ostream& s, const TIMER_BASE& t) {t.print(s); return s;}
/*--------------------------------------------------------------------------*/
class DOUBLE_TIMER : public TIMER_BASE { //
protected:
  DOUBLE_TIMER(const DOUBLE_TIMER& t)
    : TIMER_BASE(t),
      _ref(t._ref),
      _last(t._last),
      _total(t._total)
  { untested();
  }
public:
  explicit	DOUBLE_TIMER(const std::string& s) : TIMER_BASE(s)
  {
    fullreset();
  }
  explicit	DOUBLE_TIMER() : TIMER_BASE("noname")
  {
    fullreset();
  }
  ~DOUBLE_TIMER() {}
  TIMER_BASE* clone() const{return new DOUBLE_TIMER(*this);}
  void fullreset()
  {
    _total=.0;
    reset();
  }
  void reset() {
    _last = .0;
    _ref  = .0;
    _running = false;
  }
  void start() {
    assert(!_running);
    if (_running) { untested();
      stop();
    }
    _ref = run_time();
    _running = true;
  }
  void stop(){
    if (_running) {
      double runtime = run_time() - _ref;
      _ref = 0.0;
      _last  += runtime;
      _total += runtime;
      _running = false;
    }
  }
  void check() {
    if (_running) {
      stop();
      start();
    }else{
      untested();
    }
  }
  double elapsed()const	{return _last;}
  double total()const	{return _total;}
  bool is_running()const {return _running;}

  template<class S>
  void 	operator=(const TIMER_BASE*){ incomplete();}
//  friend DOUBLE_TIMER operator-(const DOUBLE_TIMER&,const DOUBLE_TIMER&);
private:
  double run_time() const
  { itested();
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return static_cast<double>(clock()) / static_cast<double>(CLOCKS_PER_SEC);
#endif
  }
  double _ref;		// time the clock was started
  double _last;		// time of timed operation
  double _total;	// time since program start
  bool	 _running;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
