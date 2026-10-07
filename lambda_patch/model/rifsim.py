# 통계 모델: 2x CIR, HPRF 펄스 간격 8탭, k0=127
# 탭 k=k0+8m+e (e in -3..4): C = A e^{jphi} g(e) X_m + thermal
#  X_m = sum_n b_n a_{n+m} e^{j theta_n}  (실수 ±1 시퀀스 곱의 합 -> 실수 가우시안을 CFO로 회전)
import math, cmath, random, sys
from multiprocessing import Pool

G_DB = {0:6.1, 1:1.8, 2:-10.3, 3:-7.3, 4:-18.5, -3:-14.0, -2:-9.7, -1:2.8}   # MD 5.1 표
G = {e: math.sqrt(10**(v/10)) for e,v in G_DB.items()}
K0=127
WS=list(range(119,176)); WN=list(range(16,97))
TAPS=list(range(14,178))
def me(k):
    d=k-K0; m=(d+3)//8; return m, d-8*m
TAPINFO={k:me(k) for k in TAPS}
MS=sorted({TAPINFO[k][0] for k in TAPS})
TL=10**0.96; TH=10**1.44; DB=10.0
FC=7987.2e6; TS=512/499.2e6

def cfo_cov(N, ppm):
    th=[2*math.pi*ppm*1e-6*FC*TS*s for s in range(N-1)]
    c=sum(math.cos(t)**2 for t in th); s2=sum(math.sin(t)**2 for t in th); cs=sum(math.cos(t)*math.sin(t) for t in th)
    # 실수 Y_s~N(0,128) 를 e^{j th_s} 로 회전해 합 -> [Re,Im] 공분산 128*[[c,cs],[cs,s2]]
    a=math.sqrt(128*c); b=128*cs/a; d=math.sqrt(max(128*s2-b*b,0.0))
    coh=128*complex(sum(math.cos(t) for t in th), sum(math.sin(t) for t in th))
    return (a,b,d), coh

def trial(rng, N, A, ppm, h1, hw, chol, coh, circ_self=False):
    M=(N-1)*128
    phi=rng.uniform(0,2*math.pi); rot=cmath.exp(1j*phi)
    a,b,d=chol
    X={}
    for m in MS:
        if circ_self:   # 비교용: 원형 복소 자기잡음
            X[m]=complex(rng.gauss(0,1),rng.gauss(0,1))*math.sqrt(M/2)
        else:
            u1=rng.gauss(0,1); u2=rng.gauss(0,1)
            X[m]=complex(a*u1, b*u1+d*u2)
    if h1:
        for m in MS:
            if m<0 and -m in X: X[m]=X[-m]      # 맞는 키: 자기상관 대칭
        X[0]=coh
    # 열잡음 (필터 hw, 탭당 분산 M)
    nh=len(hw); sc=math.sqrt(M/2/sum(x*x for x in hw))
    u=[complex(rng.gauss(0,1),rng.gauss(0,1)) for _ in range(len(TAPS)+nh)]
    P={}
    for i,k in enumerate(TAPS):
        if k<16 or k>175: continue
        w=sum(hw[j]*u[i+j] for j in range(nh))*sc
        m,e=TAPINFO[k]
        c=A*rot*G[e]*X[m]+w
        P[k]=c.real*c.real+c.imag*c.imag
    F=[0.0]*8; n=[0]*8
    for k in WN: F[k%8]+=P[k]; n[k%8]+=1
    F=[F[i]/n[i] for i in range(8)]
    Fmax=max(F); Fmin=min(F)
    Pk=max(P[k] for k in WS)
    Fmean=sum(P[k] for k in WN)/len(WN)
    Zpc=max(P[k]/F[k%8] for k in WS)
    return Pk/Fmax, Fmax/Fmin, Pk/Fmean, Zpc

def job(args):
    N, AdB, ppm, h1, hwname, ntr, seed, circ = args
    hw={'white':[1.0],'rx':[0.5,1.0,0.5]}[hwname]
    rng=random.Random(seed)
    chol,coh=cfo_cov(N,ppm)
    A=math.sqrt(10**(AdB/10)) if AdB is not None else 0.0
    return [trial(rng,N,A,ppm,h1,hw,chol,coh,circ) for _ in range(ntr)]

def run(N, AdB, ppm=0.25, h1=False, hw='rx', ntr=60000, seed=1, circ=False, procs=6):
    per=ntr//procs
    with Pool(procs) as p:
        res=p.map(job,[(N,AdB,ppm,h1,hw,per,seed*1000+i,circ) for i in range(procs)])
    return [r for part in res for r in part]

def q(xs, p):
    s=sorted(xs); return s[min(len(s)-1,int(p*len(s)))]
def dB(x): return 10*math.log10(x) if x>0 else -99
def md_valid(Z,D): return Z >= (TL if D<DB else TH)
