# 파형 수준 검증: 실제 ±1 STS(틀린 키) + 펄스 + CFO + 색잡음 -> 2x 상관 -> 원시 샘플 Q/Π -> Λ
import math, cmath, random, operator, sys
from multiprocessing import Pool
from rifsim import G
K0=127; WS=list(range(119,176)); WN=list(range(16,97)); HW=[0.5,1.0,0.5]
FS=998.4e6; FC=7987.2e6; RHO_MAX=0.98
GPH={e:G[e]*cmath.exp(1j*0.3*e) for e in G}
TL=10**0.96; TH=10**1.44
def cls(k): return ((k-K0+3)%8)-3
def frag(rng,N,A,ppm,rotT):
    M=(N-1)*128; Ltx=N*128; T=8*(M-1)+256+8
    a=[1 if rng.random()<0.5 else -1 for _ in range(Ltx)]
    b=[1 if rng.random()<0.5 else -1 for _ in range(M)]      # 틀린 키(독립)
    r=[0j]*T
    if A>0:
        for l in range(Ltx):
            pos=8*(l-128)+K0
            if pos+4<0 or pos-3>=T: continue
            al=a[l]
            for e,ge in GPH.items():
                t=pos+e
                if 0<=t<T: r[t]+=al*ge
        c0=A*cmath.exp(1j*rng.uniform(0,2*math.pi))
        r=[v*c0*rotT[t] for t,v in enumerate(r)]
    s=1/math.sqrt(2*sum(h*h for h in HW))
    u=[complex(rng.gauss(0,1),rng.gauss(0,1)) for _ in range(T+2)]
    r=[r[t]+(0.5*u[t]+u[t+1]+0.5*u[t+2])*s for t in range(T)]
    C={k:sum(map(operator.mul,b,r[k:k+8*M:8])) for k in WS+WN}
    Q={}; Pi={}
    for e in range(-3,5):
        seg=r[K0+e:K0+e+8*M:8]
        Q[e]=sum(v.real*v.real+v.imag*v.imag for v in seg); Pi[e]=sum(v*v for v in seg)
    lam=[]
    for k in WS:
        q=Q[cls(k)]; p=Pi[cls(k)]; c=C[k]; rho=abs(p)/q
        if rho>RHO_MAX: p=p*(RHO_MAX/rho)
        lam.append(2*(q*abs(c)**2-(p.conjugate()*c*c).real)/(q*q-abs(p)**2))
    P={k:abs(C[k])**2 for k in C}
    F=[0.0]*8; n=[0]*8
    for k in WN: F[k%8]+=P[k]; n[k%8]+=1
    F=[F[i]/n[i] for i in range(8)]
    zmd=max(P[k] for k in WS)/max(F); D=max(F)/min(F)
    return lam,zmd,D
def job(args):
    N,AdB,ppm,ntr,seed=args; rng=random.Random(seed)
    M=(N-1)*128; T=8*(M-1)+256+8; w=2*math.pi*ppm*1e-6*FC/FS
    rotT=[cmath.exp(1j*w*t) for t in range(T)]
    A=math.sqrt(10**(AdB/10)) if AdB is not None else 0.0
    ths=[10,14,18,2*math.log(57/1e-3),2*math.log(57/1e-3)+2*math.log(1.5)]
    cnt=[0]*len(ths); ntap=0; fa=0; fa_md=0; Tdes=ths[-1]
    for _ in range(ntr):
        lam,zmd,D=frag(rng,N,A,ppm,rotT)
        ntap+=len(lam)
        for i,t in enumerate(ths): cnt[i]+=sum(1 for x in lam if x>t)
        fa+= max(lam)>=Tdes; fa_md+= zmd>=(TL if D<10 else TH)
    return cnt,ntap,fa,fa_md,ntr
if __name__=='__main__':
    ntr=int(sys.argv[1])
    cases=[('열잡음',None,0.25),('틀린키 -82.5dBm (중간, CFO.25)',0,0.25),('틀린키 -42.5dBm (강, CFO 0)',40,0.0)]
    ths=[10,14,18,2*math.log(57/1e-3),2*math.log(57/1e-3)+2*math.log(1.5)]
    for i,(lab,AdB,ppm) in enumerate(cases):
        with Pool(6) as p: res=p.map(job,[(32,AdB,ppm,ntr//6,9000+100*i+j) for j in range(6)])
        cnt=[sum(r[0][t] for r in res) for t in range(len(ths))]; ntap=sum(r[1] for r in res)
        fa=sum(r[2] for r in res); fmd=sum(r[3] for r in res); n=sum(r[4] for r in res)
        tail=" ".join(f"t={t:.1f}: {c/ntap:.2e}(χ²₂ {math.exp(-t/2):.2e})" for t,c in zip(ths,cnt))
        print(f"[{lab}] n={n} | 탭별 꼬리 P(Λ>t): {tail}", flush=True)
        print(f"   프래그먼트 Pfa: 제안 Λ≥{ths[-1]:.2f} -> {fa/n:.1e} ({fa}건) | MD 규칙 -> {fmd/n:.1e} ({fmd}건)", flush=True)
