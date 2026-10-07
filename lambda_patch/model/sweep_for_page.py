from rifsim4 import *
import math, json
T=2*math.log(57/1e-3)+2*math.log(1.5)
out=[]
for AdB in [None,-20,-15,-12,-10,-8,-6,-4,-2,0,2,4,6,10,20,40]:
    R=run(32,AdB,ppm=0.25,ntr=60000,seed=8800+(AdB if AdB is not None else -99))
    n=len(R)
    md=sum(1 for o in R if o['md'][0]>=(TL if o['md'][1]<10 else TH))
    la=sum(1 for o in R if o['lam']>=T)
    dmed=sorted(o['md'][1] for o in R)[n//2]
    out.append({'AdB':AdB,'dBm':(None if AdB is None else AdB-82.5),'n':n,'kMD':md,'kL':la,'Dmed_dB':10*math.log10(dmed)})
    print(out[-1],flush=True)
json.dump(out,open('sweep_for_page.json','w'))
