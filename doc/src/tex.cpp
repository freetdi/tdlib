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

template <typename OUT,
	typename AUTO1,
	typename AUTO2,
	typename AUTO3,
	typename AUTO4,
	typename AUTO5,
	typename AUTO6,
	typename AUTO7>
void tex_tabular_entry(OUT &out,
		AUTO1 &arg1,
		AUTO2 &arg2,
		AUTO3 &arg3,
		AUTO4 &arg4,
		AUTO5 &arg5,
		AUTO6 &arg6,
		AUTO7 &arg7)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " & "
        << arg6 << " & "
        << arg7 << " \\\\" << std::endl;
}

template <typename OUT,
	typename AUTO1,
	typename AUTO2,
	typename AUTO3,
	typename AUTO4,
	typename AUTO5,
	typename AUTO6>
void tex_tabular_entry(OUT &out,
		AUTO1 &arg1,
		AUTO2 &arg2,
		AUTO3 &arg3,
		AUTO4 &arg4,
		AUTO5 &arg5,
		AUTO6 &arg6)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " & "
        << arg6 << " \\\\" << std::endl;
}

template <typename OUT,
	typename AUTO1,
	typename AUTO2,
	typename AUTO3,
	typename AUTO4,
	typename AUTO5>
void tex_tabular_entry(OUT &out,
		AUTO1 &arg1,
		AUTO2 &arg2,
		AUTO3 &arg3,
		AUTO4 &arg4,
		AUTO5 &arg5)
{
    out << arg1 << " & "
        << arg2 << " & "
        << arg3 << " & "
        << arg4 << " & "
        << arg5 << " \\\\" << std::endl;
}

template <typename OUT,
	typename AUTO1,
	typename AUTO2,
	typename AUTO3,
	typename AUTO4,
	typename AUTO5,
	typename AUTO6,
	typename AUTO7,
	typename AUTO8,
	typename AUTO9>
void tex_tabular_entry(OUT &out,
		AUTO1 &arg1,
		AUTO2 &arg2,
		AUTO3 &arg3,
		AUTO4 &arg4,
		AUTO5 &arg5,
		AUTO6 &arg6,
		AUTO7 &arg7,
		AUTO8 &arg8,
		AUTO9 &arg9)
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



