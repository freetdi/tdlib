namespace treedec{

template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t &T,
                   unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
{
    return minDegree_decomp(
        G, T,
        (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL,
        ub, ignore_isolated_vertices
       );
}

template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
   minDegree_decomp(G_t &G)
{
    return minDegree_decomp(G, (typename treedec_chooser<G_t>::type*)NULL);
}


template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  fillIn_decomp(G_t &G, T_t &T,
                unsigned ub = UINT_MAX, bool ignore_isolated_vertices=false)
{
    return fillIn_decomp
      (G, T,
       (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL,
       ub, ignore_isolated_vertices);
}

template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
   fillIn_decomp(G_t &G)
{
    return fillIn_decomp(G, (typename treedec_chooser<G_t>::type*)NULL);
}

} //namespace treedec
