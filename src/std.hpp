
#if __cplusplus >= 201103L
#define MOVE(x) std::move(x)
#else
#define MOVE(x) x
#endif
