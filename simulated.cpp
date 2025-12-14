#pragma GCC optimize("Ofast,no-stack-protector")
#include <bits/stdc++.h>
using namespace std;
#define endl "\n"
#define int long long
using ull=unsigned long long;
using pii=pair<int,int>;
const int dx[4] = {1,0,-1,0}, dy[4] = {0,1,0,-1};
#define OVL(x,s) for(auto y:x) cout<<y<<s; cout<<"\n";
template <typename T> istream& operator>>(istream& is, vector<T> &a) {
    copy_n(istream_iterator<T>(is), a.size(), a.begin()); return is;}
#ifdef IOI
template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream &os, const T_container &v) { os << '{'; string sep; for (const T &x : v) os << sep << x, sep = ", "; return os << '}'; }
 
void dbg_out() { cout << endl; }
template<typename Head, typename... Tail> void dbg_out(Head H, Tail... T) { cout << ' ' << H; dbg_out(T...); }
#define dbg(...) cout << "(" << #__VA_ARGS__ << "):", dbg_out(__VA_ARGS__);
#else
#define dbg(...) 1337;
#endif
#define pb push_back
#define F first
#define S second
#define all(v) v.begin(),v.end()

const int maxn = 2e5+100;
const int mod = 1e9+7;


double tle = 4000.0;  // Time limit of execution in ms. Intially set at 4seconds, but change it. The higher the better


int n, m;
vector<pair<int,int>> edges;

bool P(double E,double E_next,double T,mt19937& rng){
    if(E_next<=E) return true;
    if(T<=0) return false;
    double prob =  exp(-(E_next-E)/T);
    if(prob > 1.0) return true;
    else{
        bernoulli_distribution d(prob);
        return d(rng);
    }
}


class state {
    public:
    vector<int> perm;
    vector<int> pos;
    mt19937 mt{ static_cast<mt19937::result_type>(
        chrono::steady_clock::now().time_since_epoch().count()
        )};

    state() {
        // Generate the initial state
        perm.resize(n);
        pos.resize(n);
        for(int i = 0;i!=n;i++){
            perm[i]=i;
        }
        shuffle(perm.begin(),perm.end(),mt);
        for (int i = 0; i < n; i++) pos[perm[i]]=i;
    }

    state next(double T,double beg) {
        state s_next;
        // Modify s_next to a random neighboring state
        mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
        uniform_int_distribution<int> rando(0, n-1);
        int i=rando(rng);
        int j=rando(rng);
        int radius=1+(int)((n-1)*T/(beg+1e-9));
        if (radius < 1) radius = 1;
        int l=max(0LL,i-radius);
        int r=min(n-1,i+radius);
        uniform_int_distribution<int> radint(l, r);
        j = radint(rng);
        if (j == i) j = (i + 1 < n ? i + 1 : i - 1);
        int a=i;
        int b=j;        
        s_next.perm=perm;
        s_next.pos=pos;
        int va=s_next.perm[a];
        int vb=s_next.perm[b];
        swap(s_next.perm[a],s_next.perm[b]);
        s_next.pos[va]=b;
        s_next.pos[vb]=a;
        return s_next;
    }

    double E() {
        // Implement the energy function here
        int bw = 0;
        for (auto [u,v] : edges) {
            int d=abs(pos[u]-pos[v]);
            bw = max(bw,d);
        }
        return (double)bw;
    }
};
pair<double, state> simAnneal() {
    state s = state();
    state best = s;
    double E_cur=s.E();
    double E_best=E_cur;
    auto timebeg=chrono::steady_clock::now();
    double initT=n;
    double finaltemp=1e-4;
    mt19937 rng((uint64_t)chrono::steady_clock::now().time_since_epoch().count());

    while (true) {
        auto now = chrono::steady_clock::now();
        double elapsed_ms =
        chrono::duration<double, milli>(now - timebeg).count();
        if (elapsed_ms >= tle) break;
        double progress = elapsed_ms / tle;
        double T=initT*pow(finaltemp / initT, progress);
        state nxt=s.next(T, initT);
        double E_next=nxt.E();
        if (P(E_cur,E_next,T,rng)) {
            s = nxt;
            E_cur =E_next;
            if (E_next< E_best) {
                E_best= E_next;
                best =s;
            }
        }
    }

    return {E_best, best};
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cin>>n>>m;
    for (int i = 0; i < m; i++) {
        int u,v;
        cin>>u>>v;
        --u;--v;
        edges.pb({u,v});
    }

    pair<double, state> res=simAnneal();
    double E_best = res.first;
    state best = res.second;

    cout<<(int)E_best<<endl;
    for (int i = 0; i < n; i++) {
        cout<<best.perm[i]+1<<" ";
    }
    return 0;
}
