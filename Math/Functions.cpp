#include "../Header.cpp"

using namespace chrono;
auto start1 = high_resolution_clock::now();
auto stop = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(stop - start1);
//cerr << duration.count()/1000 << "ms" << endl;

default_random_engine generator;
uniform_real_distribution<double> distribution(0,LLONG_MAX);
ll num = distribution(generator);


// files
ifstream input;
input.open("divisibility-tree.in");
input >> n;
ofstream output;
output.open("divisibility-tree.out");
output<<" ";  
output.close();

// suma subconjuntos
for i = 0 to n-1
for mask = 0 to (1«n) - 1
if (mask & (1«i))
dp(mask) += dp(mask - (1«i))

// suma divisores
for p ∈ P
for x ∈ S (DE MENOR A MAYOR)
if (p divide a x)
dp(x) += dp(x / p)


// hash pairs unrderedmap<pll,ll,hash_pair>
struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
}; 
int maxlog2(int x) //potecncia de 2 mayor que es menor o igual a x
{
    // ll maxlog = 63 - __builtin_clzll(x);
    int maxlog = 31 - __builtin_clz(x);
    return maxlog;
}
int A[10000];//Set con reset O(1),Tambien con Map
int t=1;
bool fin(int x)
{
    return A[x]==t;
}
void borrar()
{
    t++;
}
void inse(int x)
{
    A[x]=t;
}
int res(int a,int b)
{
    int c=1,d;
    for(int i=0;i<b;i++)
    {
        c*=10;
    }
    d=c*10;
    a=a-(a-a%d);
    a=a-(a%c);
    a=a/c;
    return a;
}


//propagate val in mask to all its submask
for (int i = 0; i < p; i++)
{
    for(int mask = 0; mask < (1 << p); mask++)
    {
        if((mask & (1 << i)) == 0)
            f[mask] += f[mask | (1 << i)];
        if(mask & (1 << i)) // to propagate from submasks to mask
            dp[mask] += dp[mask - (1 << i)];
    }
}

// O(3^n)
rep(m, (1 << n)){
    // 2^k k: number of bits in m
    // iterates over al submasks of m in descending order of value
    for(int s = m; ; s = (s-1) & m){
        cout << s << endl;
        if(s == 0) break;
    }
}

int bit_opst(ll N,ll S)
{
    //Obtain the remainder (modulo) of S when it is divided by N (N is a power of 2)
    return S &(N -1);
    //Determine if S is a power of 2.
    return (S &( S -1)) == 0;
    //Turn oﬀ the last bit in S, e.g.S = (40)10 = (101000)2 →S = (32)10 = (100000)2.
    return S &( S -1);
    // Turn on the last zero in S, e.g.S = (41)10 = (101001)2 →S = (43)10 = (101011)2.
    return S|| (S + 1);
    // Turn oﬀ the last consecutive run of ones in S
    return S &( S + 1);
    // Turn on the last consecutive run of zeroes in S
    return S|| (S -1);
    // Turn on all bits
    return S = (1 << 62)-1;
    // multiply by 2^N
    return S<<=N;
    // Divide by 2^N
    return S>>=N;
    // Turn on the N-th bit
    return S = S||(1<<N);
    // Check if N-th bit is on
    return (S & (1 << N));
    // Turn off the N-th bit
    return S &= ~(1 << N);
    // Alternate satatus of N-th bir
    return S ^= (1 << N);
    //Value of the least significant bit on from the right
    return N = (S&(-S));
}
//count numbers with i bit set [1, n-1]
ll kol(ll n, ll i)
{
    return (n / (1ll << (i + 1))) * (1ll << i) + max((n % (1ll << (i + 1))) - (1ll << i), 0ll);
}
kol(r+1, i) - kol(l, i);
