#include "../Header.cpp"

// stores smallest prime factor for every number
int spf[MAXN];
 
// Calculating SPF (Smallest Prime Factor) for every
// number till MAXN.
// Time Complexity : O(nloglogn)
void sieve()
{
    spf[1] = 1;
    for (int i=2; i<MAXN; i++)
 
        // marking smallest prime factor for every
        // number to be itself.
        spf[i] = i;
 
    // separately marking spf for every even
    // number as 2
    for (int i=4; i<MAXN; i+=2)
        spf[i] = 2;
 
    for (int i=3; i*i<MAXN; i++)
    {
        // checking if i is prime
        if (spf[i] == i)
        {
            // marking SPF for all numbers divisible by i
            for (int j=i*i; j<MAXN; j+=i)
 
                // marking spf[j] if it is not
                // previously marked
                if (spf[j]==j)
                    spf[j] = i;
        }
    }
}
 
// A O(log n) function returning primefactorization
// by dividing by smallest prime factor at every step
vector<int> getFactorization(int x)
{
    vector<int> ret;
    while (x != 1)
    {
        ret.push_back(spf[x]);
        x = x / spf[x];
    }
    return ret;
}

void primeFactors(ll n) {    
    while (n % 2 == 0) {  
        cout << 2 << " ";  
        n = n/2;  
    }  
  
    for (int i = 3; i <= sqrt(n); i = i + 2) {   
        while (n % i == 0) {  
            cout << i << " ";  
            n = n/i;  
        }  
    }   
    if (n > 2)  
        cout << n << " ";  
}


// FACTORIZATION
// 33.3% of factors O(sqrt(n) / 3)
vl trial_division3(ll n) {
    vl factorization;
    for (int d : {2, 3, 5}) {
        while (n % d == 0) {
            factorization.pb(d);
            n /= d;
        }
    }
    int increments[] = {4, 2, 4, 2, 4, 6, 2, 6};
    int i = 0;
    for (ll d = 7; d * d <= n; d += increments[i++]) {
        while (n % d == 0) {
            factorization.pb(d);
            n /= d;
        }
        if (i == 8) i = 0;
    }
    if (n > 1) factorization.pb(n);
    return factorization;
}

// Pollard p - 1
// Probabilistic to find divisors
// O(B log(B) log^2(n))
ll pollards_p_minus_1(ll n) {
    int B = 10;
    ll g = 1;
    while (B <= 1000000 && g < n) {
        ll a = 2 + rand() %  (n - 3);
        g = gcd(a, n);
        if (g > 1)
            return g;

        // compute a^M
        for (int p : primes) {
            if (p >= B) continue;
            ll p_power = 1;
            while (p_power * p <= B) p_power *= p;
            a = power(a, p_power, n);
            g = gcd(a - 1, n);
            if (g > 1 && g < n) return g;
        }
        B *= 2;
    }
    return 1;
}
// PRIMALITY DETERMINISTIC UNTIL 10^18 
// 2, 3, 5, 13, 19, 73, 193, 407521, 299210837 


// POLLARD RHO
// Pollard Rho
ll mult(ll a, ll b, ll mod) { return (__int128)a * b % mod; }

// For a, b <= 10 ^ 18
ll mult2(ll a, ll b, ll mod) {
    ll result = 0;
    while (b) {
        if (b & 1) result = (result + a) % mod;
        a = (a + a) % mod;
        b >>= 1;
    }
    return result;
}


ll f(ll x, ll c, ll mod) {
    return (mult(x, x, mod) + c) % mod;
}

ll rho(ll n, ll x0=2, ll c=1) {
    ll x = x0, y = x0, g = 1;
    while (g == 1) {
        x = f(x, c, n);
        y = f(y, c, n);
        y = f(y, c, n);
        g = gcd(abs(x - y), n);
    }
    return g;
}

// Brent Factorization
ll brent(ll n, ll x0=2, ll c=1) {
    ll x = x0, g = 1, q = 1, xs, y;
    int m = 128;
    int l = 1;
    while (g == 1) {
        y = x;
        repx(i, 1, l)  x = f(x, c, n);
        int k = 0;
        while (k < l && g == 1) {
            xs = x;
            for (int i = 0; i < m && i < l - k; i++) {
                x = f(x, c, n);
                q = mult(q, abs(y - x), n);
            }
            g = gcd(q, n);
            k += m;
        }
        l *= 2;
    }
    if (g == n) {
        do {
            xs = f(xs, c, n);
            g = gcd(abs(xs - y), n);
        } while (g == 1);
    }
    return g;
}

// Eulers Totient Function (PHI function)

// phi(n) := number of integers in [1, n] coprime to n
// phi(p) = p - 1,
// phi(p^k) = p^k - p^{k-1}
// phi(a * b) = phi(a) * phi(b) * g / phi(g) with g = gcd(a, b)
// for n O(sqrt(n)) using factorization

// Eulers theorem: a and m relatively prime -> a^phi(m) % m = 1 
// m is prime -> Fermat little theorem: a^{m-1} % m = 1
// a^n % m = a^{n % phi(m)}
// Generalization: x and m not coprime and n >= log2(n)  -> x^n % m = x^{phi(m) + [n % phi(m)]}

int phi(int n) {
    int result = n;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if (n > 1) result -= result / n;
    return result;
}

// Phi function for [1, n] in O(n log(log(n)))
void phi_1_to_n(int n) {
    vector<int> phi(n + 1);
    rep(i, n + 1) phi[i] = i;
    repx(i, 2, n + 1) {
        if (phi[i] == i) {
            for (int j = i; j <= n; j += i)
                phi[j] -= phi[j] / i;
        }
    }
}


// Divisor sum property: sum_{d|n} phi(d) = n
// phi [1, n] with the divisor sum property in O(n log(n))
void phi_1_to_n(int n) {
    vector<int> phi(n + 1);
    for (int i = 0; i <= n; i++)
        phi[i] = i;

    for (int i = 2; i <= n; i++) {
        if (phi[i] == i) {
            for (int j = i; j <= n; j += i)
                phi[j] -= phi[j] / i;
        }
    }
}

// NUMBER OF DIVISORS AND SUM OF DIVISORS
// n = p1^e1 * p2^e2 * ... * pk^ek
// d(n) := number of divisors
// d(n) = (e1 + 1) * (e2 + 1) * ... * (ek + 1)
// sigma(n) := sum of divisors
// 1 + p1 + p1^2 + ... + p1^e1 = (p1^{e1+1} - 1) / (p1 - 1)
// sigma(n) = (p1^{e1+1} - 1) / (p1 - 1) * (p2^{e2+1} - 1) / (p2 - 1) * ... * (pk^{ek+1} - 1) / (pk - 1)
// Multiplicative functions: f(a * b) = f(a) * f(b) if a and b are coprimes