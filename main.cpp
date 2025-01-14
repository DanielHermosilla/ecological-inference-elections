#include <bits/stdc++.h>

using namespace std;
using ulli = unsigned long long int;
using lli = long long int;
using ld = long double;
using vlli = vector<lli>;
using vvlli = vector<vlli>;
using vvvlli = vector<vvlli>;
using vld = vector<ld>;
using vvld = vector<vld>;
using vvvld = vector<vvld>;

//struct VectorHasher {
//    lli operator()(const vlli &V) const {
//        lli hash = V.size();
//        for(auto &i : V) {
//            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
//        }
//        return hash;
//    }
//};

void fun_sum_1(vvld x, const string& nombre, ld eps = 1.0e-7, bool exit_code = true){
    for(ulli g=x.size()-1; g<x.size(); g++){
        ld sum_of_elems = accumulate(x[g].begin(), x[g].end(), (ld) 0.0);
        if (abs(sum_of_elems - 1.0) > eps) {
            if (exit_code) {
                for (int i = 0; i < x[g].size(); i++) {
                    cout << "\n\n " << nombre
                         << "[" << to_string(g) << "]"
                         << "[" << to_string(i) << "] = "
                         << to_string(x[g][i]) << "\n\n";
                }
                cout << "EM -> FUN_SUM_1: el valor de suma de " << nombre << "[" << to_string(g) << "] = " << to_string(sum_of_elems) << " es diferente a 1.0";
                return;
            }
            else {
                cout << "EM -> FUN_SUM_1: el valor de suma de " << nombre << "[" << to_string(g) << "] = " << to_string(sum_of_elems) << " es diferente a 1.0";
            }
        }
    }
}


void fun_fix_to_sum_1(vvld x, const string& nombre, ld eps = 1.0e-3, bool exit_code = true){
    fun_sum_1(x, nombre, eps, exit_code);
    for(auto & g : x){
        ld sum_of_elems = accumulate(g.begin(), g.end(), (ld) 0.0);
        for(long double & i : g){
            i /= sum_of_elems;
        }
    }
}

vlli comb2key(vlli c, ulli I){
    ulli n = c.size();
    vlli res;
    for(int j=0; j<I; j++){
        lli sum_of_elems = 0;
        for(int i=0; i<n; i++){
            if(c[i] == j) sum_of_elems++;
        }
        res.push_back(sum_of_elems);
    }
    return res;
}

vvlli combinations_with_replacement(const vlli &nums, ulli combination_length, ulli nums_length) {
    vvlli res;
    vector<int> pos(combination_length + 1, 0);

    while (true){
        for (ulli i = combination_length; i > 0; i--) {
            if (pos[i] > nums_length - 1){ // if number spilled over: xx0(nums_length-1)xx
                pos[i - 1] += 1; // set xx1(nums_length-1)xx
                for (ulli j = i; j <= combination_length; j++)
                    pos[j] = pos[j - 1]; // set xx11..1
            }
        }
        if (pos[0] > 0) // stop condition: 1xxxx
            break;
        vlli aux;
        for(int c=1; c<pos.size(); ++c){
            aux.push_back(nums[pos[c]]);
        }
        res.push_back(aux);
        pos[combination_length] += 1; // xxxxN -> xxxxN+1
    }
    return res;
}

vvvlli create_H_g(
        vector<int> n_votantes_g,
        vector<int> limites_i,
        bool cumulative = false,
        int only_g = -1
){
    ulli I = limites_i.size();
    ulli G = n_votantes_g.size();
    vvvlli H_g;
    if(only_g < 0){
        vvlli zeros(1, vlli(I, 0));
        H_g.push_back(zeros);
        for(int g=1; g<G+1; g++){
            vvlli comb;
            vlli I_vector(I);
            iota(begin(I_vector), end(I_vector), 0);
            if(!cumulative){
                comb = combinations_with_replacement(I_vector, n_votantes_g[g-1], I_vector.size());
            }
            else {
                lli sum_of_elems = 0;
                for(int e=0; e<g; e++){
                    sum_of_elems += n_votantes_g[e];
                }
                comb = combinations_with_replacement(I_vector, sum_of_elems, I_vector.size());
            }
            vvlli aux_H_g;
            for(auto & c : comb){
                vlli key = comb2key(c, I);
                lli sum_of_keys = 0;
                for(int i=0; i<I; i++){
                    if(key[i] <= limites_i[i]) sum_of_keys++;
                    else break;
                }
                if(sum_of_keys == I){
                    aux_H_g.push_back(key);
                }
            }
            H_g.push_back(aux_H_g);
        }
    }
    return H_g;
}

//tuple<vvvlli, vector<unordered_map<vlli, lli, VectorHasher>>> create_K_g(
//        vector<int> n_votantes_g,
//        vector<int> limites_i,
//        bool cumulative = false,
//        int only_g = -1
//){
//    ulli I = limites_i.size();
//    ulli G = n_votantes_g.size();
//    vvvlli K_g;
//    vector<unordered_map<vlli, lli, VectorHasher>> indice_K_g;
//    if(only_g < 0){
//        vvlli zeros(1, vlli(I, 0));
//        K_g.push_back(zeros);
//        unordered_map<vlli, lli, VectorHasher> map_zeros;
//        map_zeros[zeros[0]] = 0;
//        indice_K_g.push_back(map_zeros);
//        for(int g=1; g<G+1; g++){
//            vvlli comb;
//            vlli I_vector(I);
//            iota(begin(I_vector), end(I_vector), 0);
//            if(!cumulative){
//                comb = combinations_with_replacement(I_vector, n_votantes_g[g-1], I_vector.size());
//            }
//            else {
//                lli sum_of_elems = 0;
//                for(int e=0; e<g; e++){
//                    sum_of_elems += n_votantes_g[e];
//                }
//                comb = combinations_with_replacement(I_vector, sum_of_elems, I_vector.size());
////                cout << "comb = " << endl;
////                for(auto c: comb){
////                    cout << "\t[ ";
////                    for(auto i: c){
////                        cout << i << " ";
////                    }
////                    cout << "]" << endl;
////                }
//            }
//            vvlli aux_K_g;
//            unordered_map<vlli, lli, VectorHasher> aux_map;
//            for(int c=0; c<comb.size(); c++){
//                vlli key = comb2key(comb[c], I);
//                lli sum_of_keys = 0;
//                for(int i=0; i<I; i++){
//                    if(key[i] <= limites_i[i]) sum_of_keys++;
//                    else break;
//                }
//                if(sum_of_keys == I){
//                    aux_K_g.push_back(key);
//                    aux_map[key] = c;
//                }
//            }
////            cout << "aux_K_g = " << endl;
////            for(auto c: aux_K_g){
////                cout << "\t[ ";
////                for(auto i: c){
////                    cout << i << " ";
////                }
////                cout << "]" << endl;
////            }
//            K_g.push_back(aux_K_g);
//            indice_K_g.push_back(aux_map);
//        }
//    }
//    return {K_g, indice_K_g};
//};

vector<int> argsort(const vector<int> v) {

    // initialize original index locations
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

vvld fun_TU_cap_v4(vvld p_gi, int mesa, const vector<int>& v_i, vector<vector<int>> b_mg, ulli G, ulli I){
    cout << "v_i = \n\t[ ";
    for(int i : v_i){
        cout << to_string(i) << " ";
    }
    cout << "]" << endl;
    cout << "b_mg[" << to_string(mesa) << "] = \n\t[ ";
    for(int i : b_mg[mesa]){
        cout << to_string(i) << " ";
    }
    cout << "]" << endl;
    vector<int> indices_ord_g = argsort(b_mg[mesa]);

    vector<int> b_g_ord(G);
    vvld p_gi_ord(G, vld(I));
    for(int g=0; g<G; g++){
        b_g_ord[g] = b_mg[mesa][indices_ord_g[g]];
        for(int i=0; i<I; i++){
            p_gi_ord[g][i] = p_gi[indices_ord_g[g]][i];
        }
    }

//    auto [K_g, indices_K_g] = create_K_g(b_g_ord, v_i, true);
    vvvlli K_g = create_H_g(b_g_ord, v_i, true);

    vvld T_gk(G+1);
    for(int f=0; f<G+1; f++){
        vld aux_T_gk(K_g[f].size());
        T_gk[f] = aux_T_gk;
    }
    T_gk[0][0] = 1.0;
    vector<vvvld> U_igfk(I-1, vvvld(G, vvld(G+1)));
    for(int i=0; i<I-1; i++){
        for(int g=0; g<G; g++){
            for(int f=0; f<G+1; f++){
                vld aux_U_k(K_g[f].size());
                for(int k=0; k<K_g[f].size(); k++){
                    if((k == 0) && (f == 0) && (g == 0)) aux_U_k[k] = 1.0;
                }
                U_igfk[i][g][f] = aux_U_k;
            }
        }
    }
    int N = *max_element(b_g_ord.begin(), b_g_ord.end()) + 1;
    vld lgac_n(N);
    for(int n=0; n<N; n++){
        ld sum = 0.0;
        for(int k=0; k<n+1; k++){
            sum += log(max(k,1));
        }
        lgac_n[n] = sum;
    }

    vvvld hlg_p_gin(G, vvld(I, vld(N)));
    for(int g=0; g<G; g++){
        for(int i=0; i<I; i++){
            for(int n=0; n<N; n++){
                hlg_p_gin[g][i][n] = ((ld) n) * log(p_gi_ord[g][i]);
            }
        }
    }

    vvld hb_gn(G, vld(N));
    for(int g=0; g<G; g++){
        for (int n=0; n<N; n++) {
            if(b_g_ord[g] > 0) hb_gn[g][n] = ((ld) n)/b_g_ord[g];
        }
    }

    vvvlli H_g = create_H_g(b_g_ord, v_i, false);

    for(int f=1; f<G+1; f++){
        //unordered_map<vlli, lli, VectorHasher> dict_indices;
        map<vlli, lli> dict_indices;
        for(int k_tuple=0; k_tuple<K_g[f-1].size(); k_tuple++){
            dict_indices[K_g[f-1][k_tuple]] = k_tuple;
        }
        for(int index_k=0; index_k<K_g[f].size(); index_k++){
            vlli k = K_g[f][index_k];
            ld suma_T = 0.0;
            vvld suma_U_ig(I, vld(f));
            for (int index_h=0; index_h<H_g[f].size(); index_h++){
                vlli h = H_g[f][index_h];
                lli sum_of_h = 0;
                for(int i=0; i<I; i++){
                    if(h[i] <= k[i]) sum_of_h++;
                    else break;
                }
                if(sum_of_h < I){
                    continue;
                }
                vlli resta_tuplas(h.size());
                for(int i=0; i<h.size(); i++){
                    resta_tuplas[i] = k[i] - h[i];
                }
                lli indice = dict_indices[resta_tuplas];
//                int indice = indices_K_g[f-1][resta_tuplas];

                ld suma_hlg_p_gin = 0.0;
                for (int i=0; i<I; i++){
                    suma_hlg_p_gin += hlg_p_gin[f-1][i][h[i]] - lgac_n[h[i]];
                }
                ld a = exp(lgac_n[b_mg[mesa][f-1]] + suma_hlg_p_gin);
                suma_T += T_gk[f - 1][indice] * a;
                for(int i=0; i<I-1; i++){
                    ld a_h_b = a * hb_gn[f - 1][h[i]];
                    for(int g=0; g<f-1; g++){
                        suma_U_ig[i][g] += U_igfk[i][g][f - 1][indice] * a;
                    }
                    if(h[i] > 0){
                        int g = f-1;
                        suma_U_ig[i][g] += U_igfk[i][g][f - 1][indice] * a_h_b;
                    }
                }
            }
            T_gk[f][index_k] = suma_T;
            for(int i=0; i<I-1; i++){
                for(int g=0; g<f; g++){
                    U_igfk[i][g][f][index_k] = suma_U_ig[i][g];
                }
            }
            if(f<G){
                for(int i=0; i<I-1; i++){
                    U_igfk[i][f][f][index_k] = suma_T;
                }
            }
        }
    }

//    int indice = indices_K_g[G][v_i];
    int indice = 0; // TODO: fix
    for (int g=0; g<G; g++) {
        if(b_g_ord[g] == 0){
            for(int i=0; i<I-1; i++){
                U_igfk[i][g][G][indice] = T_gk[G][indice] / ((ld) I);
            }
        }
    }
    vld U_I(G);
    for(int g=0; g<G; g++){
        ld expresion = 0.0;
        for(int i=0; i<I-1; i++){
            expresion += U_igfk[i][g][G][indice] / T_gk[G][indice];
        }
        U_I[g] = 1-expresion;
    }
    vector<int> nuevos = argsort(indices_ord_g);
    vvld res;
    for (auto g: nuevos){
        vld aux_res;
        for (int i=0; i<I; i++){
            if(i < I-1) aux_res.push_back(U_igfk[i][g][G][indice]/T_gk[G][indice]);
            else aux_res.push_back(U_I[g]);
        }
        res.push_back(aux_res);
    }
    return res;
}

tuple<vvld, vvvld, ld, vld, int> EM(
        vector<vector<int>> b_mg,
        vector<vector<int>> n_mi,
        bool printear = false,
        ld epsilon = 1.0e-4,
        bool printear_resultados = true
){
    cout << "\n############ EM ALGORITHM ############" << endl;

    // Inicializacion
    ulli G = b_mg[0].size();
    ulli M = b_mg.size();
    ulli I = n_mi[0].size();
    int iteracion = 0;
    ld ll = -1.0e10;
    vld ll_m(M, 0.0);
    ld delta = 1.0e+10;
    vvld p_gi_old(G, vld(I));
    vvld p_gi(G, vld(I, 1.0/((ld)I)));
    if(printear){
        cout << "iteracion = " << iteracion << endl;
        cout << "p_gi = " << endl;
        for(const auto& i: p_gi){
            cout << "\t";
            for(auto j: i){
                cout << to_string(j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    vvvld q_mgi(M);
    while (delta > epsilon){
        time_t start, end;

        time(&start);

        ios_base::sync_with_stdio(false);

        for(int m=0; m<M; m++){
            q_mgi[m] = fun_TU_cap_v4(p_gi, m, n_mi[m], b_mg, G, I);
        }

        time(&end);

        cout << "\nTiempo ejecucion E-Step = " << to_string(ld(end-start)) << "\n" << endl; // Calcular tiempo
        for(int m=0; m<M; m++){
            // fun_sum_1(q_mgi[m], "q_mgi["+str(m)+"]", eps = 1e-7, exit_code = True)
            fun_fix_to_sum_1(q_mgi[m], "q_mgi["+to_string(m)+"]", 1.0e-2, true);
        }
        p_gi_old = p_gi;
        for(int g=0; g<G; g++){
            for (int i=0; i<I; i++) {
                ld sum_b = 0;
                ld sum_bq = 0;
                for (int m=0; m<M; m++) {
                    sum_b += b_mg[m][g];
                    sum_bq += (b_mg[m][g] * q_mgi[m][g][i]);
                }
                p_gi[g][i] = sum_bq/sum_b;
            }
        }
        ld sum_ll_m = 0;
        for (int m=0; m<M; m++) {
            ld sum_b_mg = 0;
            for (int g=0; g<G; g++) {
                ld aux_sum_b_mg = 0;
                for (int i=0; i<I; i++){
                    aux_sum_b_mg += (q_mgi[m][g][i] - p_gi[g][i]) * log(p_gi[g][i]);
                }
                sum_b_mg += b_mg[m][g] * aux_sum_b_mg;
            }
            ll_m[m] = sum_b_mg;
            sum_ll_m += sum_b_mg;
        }
        ll = sum_ll_m;
        iteracion++;

        ld max_diff = -1.0;
        for (int i=0; i<p_gi.size(); i++) {
            for (int j=0; j<p_gi[i].size(); j++) {
                max_diff = max(abs(p_gi[i][j] - p_gi_old[i][j]), max_diff);
            }
        }
        delta = max_diff;
        if(printear){
            cout << "iteracion = " << to_string(iteracion) << endl;
            cout << "ll = " << to_string(ll) << endl;
            cout << "delta = " << to_string(delta) << endl;
            cout << "p_gi = " << endl;
            vld suma_p;
            for(const auto& i: p_gi){
                cout << "\t";
                ld suma_p_gi = 0.0;
                for(auto j: i){
                    cout << to_string(j) << " ";
                    suma_p_gi += j;
                }
                suma_p.push_back(suma_p_gi);
                cout << endl;
            }
            cout << endl;
            for (int i=0; i<suma_p.size(); i++){
                cout << "p_gi[" << to_string(i) << "] = " << to_string(suma_p[i]) << endl;
            }
            cout << endl;
        }
        fun_sum_1(p_gi, "p_gi", 1.0e-7, true);
    }
    if(printear_resultados){
        cout << endl;
        cout << "p_gi = " << endl;
        for(const auto& i: p_gi){
            cout << "\t";
            for(auto j: i){
                cout << to_string(j) << " ";
            }
            cout << endl;
        }
        cout << endl;
        cout << "ll_m = \n\t[ ";
        for(auto i: ll_m){
            cout << to_string(i) << " ";
        }
        cout << "]" << endl;
    }
    return {p_gi, q_mgi, ll, ll_m, iteracion};
}

int main() {
    int M = 50;
    int J = 30;
    ld lambda = 0.8;
//    vvld p_gi{
//        {0.2, 0.3, 0.1, 0.4},
//        {0.3, 0.2, 0.1, 0.4},
//        {0.1, 0.1, 0.2, 0.6},
//        {0.2, 0.2, 0.3, 0.3},
//        {0.8, 0.1, 0.05, 0.05}
//    };

    vector<vector<int>> n_mi{
            {16, 10, 7, 17},
            {14, 12, 7, 17},
            {12, 10, 8, 20},
            {16, 10, 8, 16},
            {14, 6, 12, 18}
    };

    vector<vector<int>> b_mg{
            {17, 12, 7, 10, 4},
            {3, 23, 5, 13, 6},
            {7, 7, 16, 18, 2},
            {3, 18, 4, 20, 5},
            {7, 9, 10, 22, 2}
    };

    auto [p_gi, q_mgi, ll, ll_m, iteraciones] = EM(b_mg, n_mi, true, 1.0e-3, true);

    return 0;
}
