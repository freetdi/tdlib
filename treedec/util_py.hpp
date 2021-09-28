// (c) 2021 Felix Salfelder
// GPLv3+

#ifndef TREEDEC_PYTHON_UTIL_PY_H
#define TREEDEC_PYTHON_UTIL_PY_H
// https://ctrpeach.io/posts/cpp20-string-literal-template-parameters/
//
// // is this possible with c++<20?
template<size_t N>
struct StringLiteral {
    constexpr StringLiteral(const char (&str)[N]) {
        std::copy_n(str, N, value);
    }
    char value[N];
};

#endif
