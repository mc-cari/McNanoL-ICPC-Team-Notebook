short mu[MAXN] = {0,1};
void mobius(){
	fore(i,1,MAXN)if(mu[i])for(int j=i+i;j<MAXN;j+=i)mu[j]-=mu[i];
}

ll fun(ll N) {
	ll ans=0;
	for(ll i=1;i*i<=N;i++) ans+=mu[i]*(N/(i*i));
	return ans;
}
