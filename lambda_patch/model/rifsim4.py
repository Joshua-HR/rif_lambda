# 제안 알고리즘 검증용 통계 모델 (멀티패스/CFO 포함)
#  - Q_e  = Σ|r|²  (클래스 e 원시 샘플 에너지)  = b가 랜덤일 때 Var(C[k]) 의 정확값
#  - Π_e  = Σ r²   (의사분산)                  = E[C[k]^2]
#  - Λ[k] = 2(Q|C|² - Re(Π* C²)) / (Q² - |Π|²)  (비원형 백색화, H0에서 χ²₂)
import math, cmath, random
from multiprocessing import Pool
from rifsim import G, cfo_cov
from rifsim3 import channel
K0=127; WS=list(range(119,176)); WN=list(range(16,97)); HW=[0.5,1.0,0.5]
TL=10**0.96; TH=10**1.44
RHO_MAX=0.98
def lam(c,Q,Pi):
    r=abs(Pi)/Q
    if r>RHO_MAX: Pi=Pi*(RHO_MAX/r)
    return 2*(Q*(c.real*c.real+c.imag*c.imag)-(Pi.conjugate()*c*c).real)/(Q*Q-abs(Pi)**2)
def trial(rng,N,A,h1,chol,coh,pcoh,trms,kdb,need_wn=True):
    M=(N-1)*128; a,b,d=chol; ph=rng.uniform(0,2*math.pi); rot=cmath.exp(1j*ph)
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
    u={}
    def U(k):
        if k not in u: u[k]=complex(rng.gauss(0,1),rng.gauss(0,1))
        return u[k]
    taps=(WN if need_wn else [])+WS; C={}
    for k in taps:
        dd=k-K0; s=0j
        for p in range(math.ceil((tmin-dd)/8),(tmax-dd)//8+1):
            gv=gch.get(dd+8*p)
            if gv: s+=gv*getX(p)
        C[k]=A*rot*s+(HW[0]*U(k-1)+HW[1]*U(k)+HW[2]*U(k+1))*sc
    # 원시 샘플 모멘트 (클래스별): 신호 모집단값 + 열잡음 유한표본 요동
    Q={}; Pi={}
    for e in range(8):
        sQ=sum(abs(v)**2 for t,v in gch.items() if (t-e)%8==0)
        sP=sum(v*v for t,v in gch.items() if (t-e)%8==0)
        Q[e]=M*(A*A*sQ+1)+math.sqrt(M)*rng.gauss(0,1)
        Pi[e]=A*A*rot*rot*sP*pcoh+complex(rng.gauss(0,1),rng.gauss(0,1))*math.sqrt(M)
    out={}
    P={k:abs(C[k])**2 for k in taps}
    if need_wn:
        F=[0.0]*8; n=[0]*8
        for k in WN: F[k%8]+=P[k]; n[k%8]+=1
        F=[F[i]/n[i] for i in range(8)]
        out['md']=(max(P[k] for k in WS)/max(F), max(F)/min(F))
    out['zq']=max(P[k]/Q[(k-K0)%8] for k in WS)
    out['lam']=max(lam(C[k],Q[(k-K0)%8],Pi[(k-K0)%8]) for k in WS)
    out['lam_vec']=[lam(C[k],Q[(k-K0)%8],Pi[(k-K0)%8]) for k in WS]
    return out
def job(args):
    N,AdB,ppm,h1,trms,kdb,ntr,seed,keepvec=args; rng=random.Random(seed); chol,coh=cfo_cov(N,ppm)
    th=[2*math.pi*ppm*1e-6*7987.2e6*(512/499.2e6)*s for s in range(N-1)]
    pcoh=128*complex(sum(math.cos(2*t) for t in th),sum(math.sin(2*t) for t in th))
    A=math.sqrt(10**(AdB/10)) if AdB is not None else 0.0
    res=[]
    for _ in range(ntr):
        o=trial(rng,N,A,h1,chol,coh,pcoh,trms,kdb)
        if not keepvec: o.pop('lam_vec')
        res.append(o)
    return res
def run(N,AdB,ppm=0.25,h1=False,trms=0,kdb=0,ntr=60000,seed=1,procs=6,keepvec=False):
    with Pool(procs) as p:
        res=p.map(job,[(N,AdB,ppm,h1,trms,kdb,ntr//procs,seed*1000+i,keepvec) for i in range(procs)])
    return [r for part in res for r in part]
