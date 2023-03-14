#include "../Header.cpp"

// Modular inverse of an array
// x_{i} ^ -1 = prefix_{i-1} * suffix{i+1} * (x1 * x2 * ... * x_n) ^ -1
#include "Euclid.cpp"
vl invs(vl &a, int m) {
    int n = a.size();
    if (n == 0) return {};
    vl b(n);
    ll v = 1;
    for (int i = 0; i != n; ++i) {
        b[i] = v;
        v = (v * a[i]) % m;
    }
    pll p = euclid(v, m); 
    ll x = p.ff, y = p.ss;
    x = (x % m + m) % m;
    for (int i = n - 1; i >= 0; --i) {
        b[i] = (x * b[i]) % m;
        x = (x * a[i]) % m;
    }
    return b;
}

// LINEAR CONGUENCE EQUATION
// find x such that a * x % m = b
ll linear_congruence(ll a, ll b, ll m) {
    a %= m, b %= m;
    int g = gcd(a, n);
    if(g % b) return -1;
    return (minv(a, m) * b) % m;
    // more solutions: x_i = (x + i * m) % m for all i in [0, g - 1]
}

// CRT for coprime modules
// Use garner for non coprimes modules
#include "Modular.cpp"
struct Congruence { ll a, m; };

ll chinese_remainder_theorem(vector<Congruence> const& congruences) {
    ll M = 1, solution = 0;
    for (auto const& congruence : congruences) M *= congruence.m;
    for (auto const& congruence : congruences) {
        ll a_i = congruence.a;
        ll M_i = M / congruence.m;
        ll N_i = inv(M_i, congruence.m);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}

// Benja CRT
#include "Euclid.cpp"
pll CRT(pll a, pll b)
{
    if (a.ss < b.ss) swap(a, b);
    ll x, y; tie(x, y) = euclid(a.ss, b.ss); 
    ll g = a.ss * x + b.ss * y, l = a.ss / g * b.ss;
    if ((b.ff - a.ff) % g) return {-1, -1}; // no solution
    x = (b.ff - a.ff) % b.ss * x % b.ss / g * a.ss + a.ff; 
    return {x + (x < 0) * l, l};
}

pll CRT(vector<pll> &v)
{
    int N = v.size(); pll ans = v[0];
    rep(i, N) if (i) ans = CRT(ans, v[i]);
}


// DISCRETE LOGARITHM
// Returns minimum x for which a ^ x % m = b % m in O(sqrt(m))
// Includes case when a and m are not coprime dividing by the gcd
int solve(int a, int b, int m) {
    a %= m, b %= m;
    int k = 1, add = 0, g;
    while ((g = __gcd(a, m)) > 1) {
        if (b == k)
            return add;
        if (b % g)
            return -1;
        b /= g, m /= g, ++add;
        k = (k * 1ll * a / g) % m;
    }

    int n = sqrt(m) + 1;
    int an = 1;
    rep(i, n) an = (an * 1ll * a) % m;

    unordered_map<int, int> vals;
    for (int q = 0, cur = b; q <= n; ++q) {
        vals[cur] = q;
        cur = (cur * 1ll * a) % m;
    }

    for (int p = 1, cur = k; p <= n; ++p) {
        cur = (cur * 1ll * an) % m;
        if (vals.count(cur)) {
            int ans = n * p - vals[cur] + add;
            return ans;
        }
    }
    return -1;
}

