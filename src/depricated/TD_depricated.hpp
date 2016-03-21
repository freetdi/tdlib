
template <typename G_t>
void all_k_colorings(unsigned int n, unsigned int k, unsigned int max,
        std::set<unsigned int> &M, std::vector<std::vector<int> > &colorings, G_t &G)
{
    std::vector<int> coloring(n, -1);
    std::set<unsigned int>::iterator iM = M.begin();
    while(iM != M.end()){
        coloring[*(iM++)]++;
    }

    iM = M.begin();
    unsigned int c = 0;

    if(treedec::app::is_valid_coloring(G, M, coloring)){
        colorings[c++] = coloring;
    }

    while(iM != M.end() && c < colorings.size()){
        if(coloring[*iM] < k-1){
            if(iM == M.end()){ break; }
            coloring[*iM]++;

            if(treedec::app::is_valid_coloring(G, M, coloring)){
                colorings[c++] = coloring;
            }
        }
        else{
            while(coloring[*iM] == k-1 && iM != M.end()){
                coloring[*iM] = 0;
                iM++;
            }
            if(iM == M.end()){ break; }

            coloring[*iM]++;

            if(treedec::app::is_valid_coloring(G, M, coloring)){
                colorings[c++] = coloring;
            }

            iM = M.begin();
        }
    }

    colorings.resize(c);
}
