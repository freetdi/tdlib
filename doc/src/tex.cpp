#include <fstream>
#include <string>
#include <boost/algorithm/string/replace.hpp>

template <typename OUT>
void tex_tabular_begin(OUT &out, unsigned c){
    out  << "\\begin{tabular}{|";
    for(unsigned i = 0; i < c; i++){
        out << "l|";
    }
    out << "}" << std::endl;
}

template <typename OUT>
void tex_tabular_end(OUT &out){
    out << "\\end{tabular}" << std::endl;
}

template <typename OUT>
void tex_hline(OUT &out){
    out << "\\hline" << std::endl;
}

template <typename OUT>
void tex_tabular_entry(OUT &out, auto &arg1, auto &arg2, auto &arg3, auto &arg4, auto &arg5, auto &arg6, auto &arg7)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " & "
        << arg6 << " & "
        << arg7 << " \\\\" << std::endl;
}

template <typename OUT>
void tex_tabular_entry(OUT &out, auto &arg1, auto &arg2, auto &arg3, auto &arg4, auto &arg5, auto &arg6)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " & "
        << arg6 << " \\\\" << std::endl;
}

template <typename OUT>
void tex_tabular_entry(OUT &out, auto &arg1, auto &arg2, auto &arg3, auto &arg4, auto &arg5)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " \\\\" << std::endl;
}

template <typename OUT>
void tex_tabular_entry(OUT &out, auto &arg1, auto &arg2, auto &arg3, auto &arg4, auto &arg5, auto &arg6, auto &arg7, auto &arg8, auto &arg9)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " & "
        << arg6 << " & "
        << arg7 << " & "
        << arg8 << " & "
        << arg9 << " \\\\" << std::endl;
}



