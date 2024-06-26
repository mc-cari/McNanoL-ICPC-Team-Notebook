#include "../Header.cpp"

struct RH {
    int binpow(int a, int b, int MOD) {
        int res = 1;
        while(b) {
            if(b & 1) res = (1LL * res * a) % MOD;
            a = (1LL * a * a) % MOD;
            b >>= 1;
        }
        return res;
    }
    // choose base B random to avoid hacks 33 37 41
    // randomize V(s[i])
    const int B = 119, M[2] = {999727999, 1070777777}, P[2] = {325255434, 10018302};
    vector<int> H[2], I[2];
    RH(string &s) {
        int N = s.size(); 
        rep(k, 2) {
            H[k].resize(N + 1), I[k].resize(N + 1);
            H[k][0] = 0; ll b = 1;
            rep(i, N) {
                H[k][i + 1] = (H[k][i] + b * (s[i] - 'a' + 1)) % M[k];
                b = (b * B) % M[k];
            }
            I[k][N] = binpow(b, M[k] - 2, M[k]);
            for(int i = N - 1; i >= 0; i--) I[k][i] = (1LL * I[k][i + 1] * B) % M[k]; 
        }
    }
    ll get(int l, int r)  // inclusive - exclusive
    {
        ll h0 = (H[0][r] - H[0][l] + M[0]) % M[0];
        h0 = (1LL * h0 * I[0][l]) % M[0];
        ll h1 = (H[1][r] - H[1][l] + M[1]) % M[1];
        h1 = (1LL * h1 * I[1][l]) % M[1];
        return (h0 << 32) | h1;
    }
};
bool compare(int a, int b){
    ll l = 0, r = n-1, p, res = -1;
    while(l <= r){
        p = (l + r) / 2;
        if(rhs[a].get(0, p) == rhs[b].get(0, p))l = p+1;
        else {
            res = p;
            r = p-1;
        }
    }
    if(res == -1)return a < b;
    return s[a][res] < s[b][res];
}

//Suffix Array O(N log^2 N)
rep(n) sa[i] = i;
sort(ALL(sa), compare)

template<const int LIM, const int MOD, const int BASE>
struct Hash {
    vector<int> H, P;
 
    int binpow(int a, int b) {
        int res = 1;
        while(b) {
            if(b & 1) res = (1LL * res * a) % MOD;
            a = (1LL * a * a) % MOD;
            b >>= 1;
        }
        return res;
    }
 
    void init(const string& s) {
        H.resize(LIM); P.resize(LIM);
        int n = s.length(), b = 1;
        H[0] = 0;
        rep(i, n) {
            H[i + 1] = (H[i] + 1LL * b * (s[i] - 'a' + 1)) % MOD;
            b = (1LL * b * BASE) % MOD;
        }
        P[n] = binpow(b, MOD - 2);
        for(int i = n - 1; i >= 0; i--) P[i] = (1LL * P[i + 1] * BASE) % MOD; 
    }
    int query(int l, int r) { return (1LL * (H[r] - H[l] + MOD) * P[l]) % MOD; }
};

struct MultiHash {    
    Hash<200200, 1000000007, 123> A; 
    Hash<200200, 1000000009, 1234> B;
    Hash<200200, 1000000021, 12345> C;

    void init(const string &s) {
        A.init(s); B.init(s); C.init(s);
    }

    vector<int> query(int l, int r) {
        vector<int> q;
        q.pb(A.query(l, r)); q.pb(B.query(l, r)); q.pb(C.query(l, r));
        return q; 
    }    
};
