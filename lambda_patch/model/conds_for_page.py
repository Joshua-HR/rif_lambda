from rifsim4 import *
import math, json
T=2*math.log(57/1e-3)+2*math.log(1.5); Tq=T/2
cases=[('열잡음 (신호 없음)',None,0.25,0,0),('틀린 키 −88.5 dBm',-6,0.25,0,0),('틀린 키 −82.5 dBm',0,0.25,0,0),
       ('틀린 키 −76.5 dBm',6,0.25,0,0),('틀린 키 −42.5 dBm',40,0.25,0,0),('틀린 키 −42.5 dBm · 잔여 CFO 0',40,0.0,0,0),
       ('틀린 키 −42.5 dBm · 멀티패스 τ10',40,0.25,10,0),('틀린 키 −42.5 dBm · 멀티패스 τ20 K−6',40,0.25,20,-6),
       ('틀린 키 −82.5 dBm · 멀티패스 τ20 K−6',0,0.25,20,-6)]
out=[]
for i,(lab,AdB,ppm,trms,kdb) in enumerate(cases):
    R=run(32,AdB,ppm=ppm,trms=trms,kdb=kdb,ntr=60000,seed=9100+i)
    n=len(R)
    out.append({'label':lab,'n':n,
      'kMD':sum(1 for o in R if o['md'][0]>=(TL if o['md'][1]<10 else TH)),
      'kZq':sum(1 for o in R if o['zq']>=Tq),
      'kL':sum(1 for o in R if o['lam']>=T)})
    print(out[-1],flush=True)
json.dump(out,open('conds_for_page.json','w'),ensure_ascii=False)
