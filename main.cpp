#include <stdio.h>
#include <vector>

#define gc getchar
#define MOD 1000000007
#define MAX_SIZE 100000003
 
int fact[MAX_SIZE];
int invFact[MAX_SIZE];
int f[MAX_SIZE];
bool isprime[MAX_SIZE];
 
int qpow(int a, long long b) {
	int ans = 1;
	while (b) {
		if (b & 1) ans = 1ll * ans * a % MOD;
		a = 1ll * a * a % MOD;
		b >>= 1;
	}
	return ans;
}
 
int add(int a, int b) 
{  
	int w = a + b;
	return (w < MOD) ? w : w - MOD; 
}
 
long long read()
{
	long long value = 0; bool ne = 0;
	char c = gc();
	while(c==' ' or c=='\n') c = gc();
	if(c=='-'){ne = 1;c = gc();}
	while(c>='0' and c<='9')
	{
		value = (value<<3)+(value<<1)+c-'0'; c = gc();
	}
	if(ne) value*=-1;
	return value;
}
 
 
void precomputeFactorial()
{
	fact[0] = 1;
	for(int i = 1;i <= MAX_SIZE-1; i++){
		fact[i] = 1LL * fact[i-1] * i % MOD;
		f[i] = i;
	}
    invFact[MAX_SIZE-1] = qpow(fact[MAX_SIZE-1], MOD-2);
	for(int i = MAX_SIZE-1; i; i--)invFact[i-1] = 1LL * invFact[i] * i % MOD;
}
 
 
int C(int n, int k) {
    if (k > n || k < 0) return 0;
    return 1ll * fact[n] * invFact[k] % MOD * invFact[n - k] % MOD;
}
 
int lagrange(int n, long long z) {
	z = z % MOD;
	int prefix = z, suffix = 1;
	for (int i = 0; i <= n; i++) {
	    int j = n - i;
	    if (i > 0) {
	        prefix = 1LL * prefix * (z - i + MOD) % MOD; 
	        f[i+1] = 1ll * f[i+1] * prefix % MOD;
	    }
	    else{
	    	f[1] = 1ll * f[1] * prefix % MOD;
		}
	    suffix = 1LL * suffix * (z - j + MOD) % MOD; 
	    if (j-1 >= 0){
	    	f[j-1] = 1ll * f[j-1]* suffix % MOD;
		}
	}
 
	int ans = 0;
	for(int i = 0; i <= n; i++) {
		int mul = ((n-i) & 1) ? MOD-1 : 1;
		int b = 1ll * (i == 0 ? 1 : invFact[i]) * invFact[n-i] % MOD;
		ans= add(ans, 1ll * mul * f[i] % MOD * b % MOD);
	}
		 
	for(int i = 0; i <= n+1; i++){
		isprime[i] = false;
		f[i] = i;
	}
 
	return ans;
}
 
 
int solve_sum(long long n, long long a, long long r) {
	

	a = a % MOD;
	std::vector<int> primes;
	f[0] = 0, f[1] = 1;
	
	if(r > 100000000){
		int ans = 0;
		for(int i = 1, now = a; i <= n; i++, now = 1ll * now * a % MOD)
		{
			ans = add(ans, 1ll * now * qpow(i, r) % MOD);
		}
		return ans;
	}
	
	for(int i = 2; i <= (r+1); i++){
	    if(!isprime[i])
	    {
	        primes.push_back(i);
	        int w = qpow(i, r);
	        for(long long j = i; j <= (r+1); j *= i)
	        {
	            isprime[j] = true;
	            f[j] = 1ll * f[j/i] * w % MOD;
	        }
	    }
	    for(int p:primes)
	        if(i * p <= (r+1) && i % p)
	        {
	            for(long long t = p; t * i <= (r+1); t *= p)
	            {
	                isprime[t*i] = true;
	                f[t*i] = 1ll * f[t] * f[i] % MOD;
	            }
	        }
	    else break;
	}
 
	int inva = qpow(a,MOD-2);
 
 
	for (int i = 0, now = a, pww = 1; i <= r + 1; i++, now = 1ll * now * a % MOD, pww = 1ll * pww * (MOD - inva) % MOD) {

	    int x = 1ll * f[i + 1] * now % MOD;
	    f[i + 1] = add(f[i], x) % MOD;
	
	    int y = 1ll * C(r + 1, i) * pww % MOD * f[i] % MOD;
	    f[0] = add(f[0], y);
	}
 
	if(n < r){
		int ans = f[n];
		for(int i = 0; i <= r+1; i++){
			isprime[i] = false;
			f[i] = i;
		}
		return ans;
	}
 
 
	f[0] = 1ll * (MOD - f[0]) % MOD * qpow(qpow((1ll - inva + MOD) % MOD, r+1), MOD-2) % MOD;
	for(int i = 1, pww = inva; i <= r; i++, pww = 1ll * pww * inva % MOD){
		int temp = add(f[0], f[i]);
		f[i] = 1ll * temp * pww % MOD;
	} 
	int ans, g_n, g_0;
	if(a == 1){
		ans = lagrange(r+1,n);
	}
	else{
		g_0 = f[0];
		g_n = lagrange(r,n);
		ans = (1ll * qpow(a,n) * g_n % MOD - g_0 + MOD) % MOD;
	}
	
	return ans;
}
 
int main()
{	
	precomputeFactorial(); 
    int t = read();
    while(t)
    {
    	long long n = read();
		long long a = read();
		long long r = read();
    	printf("%d\n", solve_sum(n, a, r));
    	t--;
	}
    return 0;
}    
