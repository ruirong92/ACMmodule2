
typedef long long LL;

//////// 素数表生成 ///////

const int maxn = 10000000 + 10;
const int maxp = 700000;

int vis[maxn];//vis[i]=1,则i是合数；vis[i]=0,则i是1或者素数
int prime[maxp];

//筛素数
void sieve(int n)
{
 int m = (int)sqrt(n+0.5);//避免浮点误差
 memset(vis,0,sizeof(vis));
 for (int i = 2; i <= m; ++i)
 {
  if (!vis[i])
  {
   for (int j = i*i; j <= n; j += i)
   {
    vis[j] = 1;
   }
  }
 }
}

//生成素数表，放在prime数组中，返回素数个数
int gen_primes(int n)
{
 sieve(n);
 int c = 0;
 for (int i = 2; i <= n; ++i)
 {
  if (!vis[i])
  {
   prime[c++] = i;
  }
 }
 return c;
}



//////欧几里得算法与扩展欧几里得算法/////

//返回gcd(a,b)
LL gcd(LL a,LL b)
{
 return b == 0 ? a : gcd(b,a%b);
}

//求整数x和y,使得ax+by=d,且|x|+|y|最小。其中d=gcd(a,b)
//即使a,b在int范围内，x和y有可能超出int范围
void gcd(LL a,LL b,LL& d,LL& x, LL& y)
{
 if (!b)
 {
  d = a;
  x = 1;
  y = 0;
 }
 else
 {
  gcd(b,a%b,d,y,x);
  y -= x*(a/b);
 }
}


//二进制的gcd算法
/*
二进制最大公约数基于下述事实：
若a、b都是偶数，则gcd(a,b)=2*gcd(a/2,b/2)
若a是奇数、b是偶数，则gcd(a,b)=gcd(a/2,b/2)
若a、b都是奇数，则gcd(a,b)=gcd((a-b)/2,b)
*/

LL binary_gcd(LL x,LL y)
{
  if(!x)
    return y;
  if(!y)
    return x;
  LL k = 0;
  while(!(x&1) && !(y&1))
  {
    x >>= 1;
    y >>= 1;
    k ++;
  }
  LL t = -v;
  if(!(x&1))
    t = u;
  while(t)
  {
    while(!(t&1))
      t >>= 1;
    if(t > 0)
      u = t;
    else
      v = -t;
    t = u - v;
  }
  return u * (1 << k);
}


/////////模运算////////////

//返回 ab mod n 。 要求0<=a,b<n
LL mul_mod(LL a,LL b,int n)
{
 return a * b % n;
}

//返回 a^p mod n ,要求 0<=a<n
LL pow_mod(LL a,LL p,LL n)
{
 if (!p)
 {
  return 1;
 }
 LL ans = pow_mod(a,p/2,n);
 ans = ans * ans % n;
 if(p%2)
 {
  ans = ans * a % n;
 }
 return ans;
}

//计算模n下a的逆。如果不存在逆，返回-1
LL inv(LL a,LL n)
{
 LL d,x,y;
 gcd(a,n,d,x,y);
 return d == 1 ? (x+n)%n : -1;
}


/*公式：A^x % m = A^(x%phi(m)+phi(m)) % m (x >= phi(m)) */

//////////欧拉函数//////////////
/*
    φ(n)为不超过n且与n互质的正整数个数
    欧拉函数的公式为
    φ(n)=n(1-1/p1)(1-1/p2)...(1-1/p3)
*/

//计算欧拉phi函数。phi(n)为不超过n且与n互质的正整数个数
int euler_phi(int n)
{
 int m = (int)sqrt(n+0.5);
 int ans = n;
 for (int i = 2; i <= m; ++i)
 {
  if (!(n%i))
  {
   ans = ans / i * (i-1);
   while(!(n%i))
   {
    n /= i;
   }
  }
 }
 if (n > 1)
 {
  ans = ans / n * (n-1);
 }
}

//用类似筛法的方法计算phi[1],phi[2],...,phi[n]
int phi[maxn];
void phi_table(int n)
{
 memset(phi,0,sizeof(phi));
 phi[1] = 1;
 for (int i = 2; i <= n; ++i)
 {
  if (!phi[i])
  {
   for (int j = i; j <= n; j += i)
   {
    if (!phi[j])
    {
     phi[j] = j;
    }
    phi[j] = phi[j] / i * (i-1);
   }
  }
 }
}


////////////模方程//////////////

////////中国余弦定理///////////
//n个方程：x=a[i](mod m[i])(0<=i<n)
LL china(int n,int *a,int *m)
{
 LL M = 1,d,y,x=0;
 for (int i = 0; i < n; ++i)
 {
  M *= m[i];
 }
 for (int i = 0; i < n; ++i)
 {
  LL w = M / m[i];
  gcd(m[i],w,d,d,y);
  x = (x + y * w * a[i]) % M;
 }
 return (x+M) % M;
}

//求解模方程a^x = b(mod n) . n为素数，无解返回 -1(大步小步算法)
int log_mod(int a,int b,int n)
{
 int m,v,e = 1,i;
 m = (int)sqrt(n+0.5);
 v = inv(pow_mod(a,m,n),n);
 map<int,int> x;
 x[1] = 0;
 for (int i = 1; i < m; ++i)
 {
  if (x.count(b))
  {
   return i * m + x[b];
  }
  b = mul_mod(b,v,n);
 }
 return -1;
}

