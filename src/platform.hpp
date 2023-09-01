// platform/machine dependent definitions for treedec
#ifndef TREEDEC_PLATFORM_HPP
#define TREEDEC_PLATFORM_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// TODO: include system headers

#ifdef __APPLE__
#ifdef howmany
#undef howmany
#endif
#endif

#if __cplusplus >= 201103L
#define MOVE(x) std::move(x)
#else
#define MOVE(x) x
#endif

#ifndef INTERRUPTION_POINT
#define INTERRUPTION_POINT
#endif

#ifndef PROPAGATION_POINT
#define PROPAGATION_POINT
#endif

#endif
