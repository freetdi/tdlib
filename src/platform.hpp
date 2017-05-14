// platform/machine dependent definitions for tdlib
#ifndef TD_PLATFORM_HPP
#define TD_PLATFORM_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
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
