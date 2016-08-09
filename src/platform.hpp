// platform/machine dependent definitions for tdlib
#ifndef TD_PLATFORM_HPP
#define TD_PLATFORM_HPP

#if __cplusplus >= 201103L
#define MOVE(x) std::move(x)
#else
#define MOVE(x) x
#endif

#ifndef INTERRUPTION_POINT
#define INTERRUPTION_POINT
#endif

#endif
