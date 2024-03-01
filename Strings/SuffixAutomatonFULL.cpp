const int N = 3e5 + 9;

// len -> largest string length of the corresponding endpos-equivalent class
// link -> longest suffix that is another endpos-equivalent class.
// firstpos -> 1 indexed end position of the first occurrence of the largest string of that node
// minlen(v) -> smallest string of node v = len(link(v)) + 1
// terminal nodes -> store the suffixes
struct SA {
    int sz, l; 
    // C: indicador si es clone
    // Lk: link
    // Ilk: inverse links
    // L: maxlen of state
    // S: Size of state
    // F: first ocurrence of the state
    vl L, Lk, S, C, T, F;
    vector<vl> N, Ilk;

    SA(string s, int n) : L(2 * n), Lk(2 * n), C(2 * n), F(2 * n), N(2 * n, vl(26, -1)) {
        l = L[0] = 0, Lk[0] = -1, sz = 1; 
        int p;
        for(char c : s) p = extend(c - 'a');
        Ilk.resize(sz); S.assign(sz, -1);
        rep(i, sz) if(i) Ilk[Lk[i]].pb(i);

        T.assign(sz, 0); 
        while(p != -1) T[p] = 1, p = Lk[p];
    }
    int extend(int c) {
        int cur = sz++, p = l; 
        C[cur] = 0, L[cur] = L[l] + 1, F[cur] = L[cur] - 1;
        while(p != -1 && N[p][c] == -1) N[p][c] = cur, p = Lk[p];
        if(p == -1) { 
            Lk[cur] = 0, l = cur;
            return cur;
        }
        int q = N[p][c];
        if(L[p] + 1 == L[q]) { Lk[cur] = q, l = cur; return cur; }
        int w = sz++; 
        C[w] = 1, L[w] = L[p] + 1, Lk[w] = Lk[q], N[w] = N[q];
        F[w] = F[q];
        while (p != -1 && N[p][c] == q) N[p][c] = w, p = Lk[p];
        Lk[q] = Lk[cur] = w, l = cur; return cur;
    }
    int size(int p) {
        if(S[p] != -1) return S[p];
        for(int i : Ilk[p]) S[p] += size(i);
        return S[p] += (1 - C[p]) + 1;
    }

    int find(string& s) {
        int p = 0;
        for(auto u: s) {
            int x = u - 'a';
            if(N[p][x] == -1) return -1;
            p = N[p][x];
        }
        return p;
    }

    vl all_occur(string& s) {
        vl aux;
        int p = find(s);
        if(p == -1) return aux;
        vl Q = {p};
        int large = s.size();
        while(!Q.empty()) {
            int u = Q.back(); Q.pop_back();
            if(!C[u]) aux.pb(F[u] - large + 1);
            for(int v: Ilk[u]) Q.pb(v);
        }
        return aux;
    }
};

int32_t main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    int t; cin >> t;
    while (t--) {
        string s; cin >> s;
        int n = s.size();
        SA sa(s, n);
        long long ans = 0; //number of unique substrings
        for(int i = 1; i < sa.sz; i++) ans += sa.L[i] - sa.L[sa.Lk[i]];
        cout << ans << '\n';
    }
    return 0;
}
