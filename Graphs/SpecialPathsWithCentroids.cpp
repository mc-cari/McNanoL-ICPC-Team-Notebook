int n, k, sz[mxN], mxtree[mxN], ans; 
vector<pair<int, int>> adj[mxN];
vector<int> new_vis;
bool used[mxN];

int dfs(int u, int p = -1){
    sz[u] = 1;
    for(pair<int, int> v : adj[u]){
        if(v.second != p && !used[v.second]){
            sz[u] += dfs(v.second, u); 
        }
    }
    return sz[u];
}

int get_centroid(int u, int ms, int p = -1){ // u = node, ms = size of the current tree delimited by previous centroids, p = parent
    for(pair<int, int> v : adj[u]){
        if(v.second != p && !used[v.second]){
            if(sz[v.second]*2 > ms){
                return get_centroid(v.second, ms, u);
            }
        }
    }
    return u;
}

int flag[mxN], stamp[mxN];
int timestamp;

void solve(int u){
    ++timestamp;
    int sz_curr_tree = dfs(u);
    int centroid = get_centroid(u, sz_curr_tree);
    used[centroid] = 1;
    flag[0] = 0;
    stamp[0] = timestamp;
    // do something with the centroid 
    for(pair<int, int> v : adj[centroid]){
        
    }
    // end that something with the centroid xdd

    // solve for child trees centroids
    for(pair<int, int> v : adj[centroid]){
        if(used[v.second]) continue;
        solve(v.second);
    }
}
