# 멀티패스: g_ch(t) = sum_l c_l g(t - tau_l), C[k0+d] = A e^{jphi} sum_p g_ch(d+8p) X'_p
import math, cmath, random
from multiprocessing import Pool
from rifsim import G, cfo_cov
K0=127; WS=list(range(119,176)); WN=list(range(16,97)); HW=[0.5,1.0,0.5]
TL=10**0.96; TH=10**1.44
def channel(rng, trms, kfac_db):
    # LOS(지연0) + 지수 감쇠 확산 경로(1~4탭 간격 랜덤), 확산 총전력 = LOS / K
    taps={0:1.0+0j}
    if trms>0:
        pdiff=10**(-kfac_db/10); t=0; paths=[]
        while t<6*trms:
            t+=rng.randint(1,4); paths.append(t)
        w=[math.exp(-t/trms) for t in paths]; s=sum(w)
        for t,wi in zip(paths,w):
            amp=math.sqrt(pdiff*wi/s/2)
            taps[t]=taps.get(t,0)+complex(rng.gauss(0,amp),rng.gauss(0,amp))
    gch={}
    for tau,c in taps.items():
        for e,ge in G.items(): gch[tau+e]=gch.get(tau+e,0)+c*ge
    return gch
def trial(rng,N,A,h1,chol,coh,trms,kdb):
    M=(N-1)*128; a,b,d=chol; rot=cmath.exp(1j*rng.uniform(0,2*math.pi))
    gch=channel(rng,trms,kdb); tmin=min(gch); tmax=max(gch)
    X={}
    def getX(p):
        if p not in X:
            if h1 and p==0: X[p]=coh
            elif h1 and -p in X: X[p]=X[-p]
            else:
                u1=rng.gauss(0,1); u2=rng.gauss(0,1); X[p]=complex(a*u1,b*u1+d*u2)
        return X[p]
    sc=math.sqrt(M/2/sum(x*x for x in HW))
    taps=WN+WS; P={}
    u={}
    def U(k):
        if k not in u: u[k]=complex(rng.gauss(0,1),rng.gauss(0,1))
        return u[k]
    for k in taps:
        dd=k-K0; s=0j
        for p in range(math.ceil((tmin-dd)/8), (tmax-dd)//8+1):
            gv=gch.get(dd+8*p)
            if gv: s+=gv*getX(p)
        w=(HW[0]*U(k-1)+HW[1]*U(k)+HW[2]*U(k+1))*sc
        c=A*rot*s+w; P[k]=c.real*c.real+c.imag*c.imag
    F=[0.0]*8; n=[0]*8
    for k in WN: F[k%8]+=P[k]; n[k%8]+=1
    F=[F[i]/n[i] for i in range(8)]
    Pk=max(P[k] for k in WS); argk=max(WS,key=lambda k:P[k])
    return Pk/max(F), max(F)/min(F), argk
def job(args):
    N,AdB,ppm,h1,trms,kdb,ntr,seed=args; rng=random.Random(seed); chol,coh=cfo_cov(N,ppm)
    A=math.sqrt(10**(AdB/10))
    return [trial(rng,N,A,h1,chol,coh,trms,kdb) for _ in range(ntr)]
def run(N,AdB,ppm=0.25,h1=False,trms=0,kdb=0,ntr=48000,seed=1,procs=6):
    with Pool(procs) as p:
        res=p.map(job,[(N,AdB,ppm,h1,trms,kdb,ntr//procs,seed*1000+i) for i in range(procs)])
    return [r for part in res for r in part]
def q(xs,p): s=sorted(xs); return s[min(len(s)-1,int(p*len(s)))]
def dB(x): return 10*math.log10(x)
