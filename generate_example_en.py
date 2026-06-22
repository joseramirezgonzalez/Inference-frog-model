import random, math, json
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

OUT=Path('/mnt/data/frog_beta045_example_en')
K=5; T=4; beta=0.45; seed=9
random.seed(seed)

def sample_beta1b(b):
    u=random.random()
    return 1-(1-u)**(1/b)

def sample_xi(p):
    u=random.random()
    if p <= 0: return 0
    if p >= 1: return 10**18
    return int(math.floor(math.log(u)/math.log(p)))

def sample_dir():
    return 1 if random.random()<0.5 else -1

latent={}
def ensure(r,x):
    if (r,x) not in latent:
        p=sample_beta1b(beta)
        xi=sample_xi(p)
        dirs=[sample_dir() for _ in range(20)]
        latent[(r,x)]={'Pi':p,'Xi':xi,'dirs':dirs}
    return latent[(r,x)]

active={}; visited=[set([0]) for _ in range(K+1)]
order=[]; history=[]; events=[]
for r in range(1,K+1):
    ensure(r,0)
    active[(r,0)]={'A':0,'pos':0,'age':0,'alive':True,'death_time':None}
    order.append((0,r,0))
history.append({k:(v['pos'],v['alive'],v['age']) for k,v in active.items()})

for t in range(T):
    arrivals=defaultdict(list)
    for (r,x),st in list(active.items()):
        if st['A']>t or not st['alive']:
            continue
        lat=ensure(r,x); age=st['age']; old=st['pos']
        if lat['Xi']==age:
            st['alive']=False; st['death_time']=t+1
            events.append({'r':r,'x':x,'t0':t,'t1':t+1,'age':age,'type':'death','from':old,'to':old})
        else:
            d=lat['dirs'][age]
            st['pos']+=d; st['age']+=1
            arrivals[r].append(st['pos'])
            events.append({'r':r,'x':x,'t0':t,'t1':t+1,'age':age,'type':'jump','from':old,'to':st['pos'],'dir':d})
    for r,arrs in arrivals.items():
        for site in sorted(set(arrs)):
            if site not in visited[r]:
                visited[r].add(site)
                ensure(r,site)
                active[(r,site)]={'A':t+1,'pos':site,'age':0,'alive':True,'death_time':None}
                order.append((t+1,r,site))
    history.append({k:(v['pos'],v['alive'],v['age']) for k,v in active.items()})

particles=[]
for n,(a,r,x) in enumerate(order,1):
    st=active[(r,x)]; lat=latent[(r,x)]
    C=max(0,T-a)
    if a<T:
        Z=min(lat['Xi'],C)
        Delta=int(lat['Xi']<C)
    else:
        Z=0; Delta=0
    obs_dirs=lat['dirs'][:Z]
    X=[x]
    for d in obs_dirs: X.append(X[-1]+d)
    # full state-position record at t=0..T
    statepos=[]
    for t in range(T+1):
        if t<a:
            state='sl'; pos=x
        elif st['death_time'] is not None and t>=st['death_time']:
            state='de'; pos=X[-1]
        else:
            state='al'
            m=min(t-a,Z)
            pos=X[m]
        statepos.append({'t':t,'state':state,'pos':pos})
    Ys=[]; Ds=[]; Js=[]
    for j in range(T):
        Y=int(a+j<T and lat['Xi']>=j)
        D=int(a+j<T and lat['Xi']==j)
        J=Y-D
        Ys.append(Y); Ds.append(D); Js.append(J)
    particles.append(dict(n=n,r=r,x=x,A=a,Pi=lat['Pi'],Xi=lat['Xi'],C=C,Z=Z,Delta=Delta,
                          dirs=obs_dirs,X=X,statepos=statepos,Y=Ys,D=Ds,J=Js,
                          alive_T=statepos[-1]['state']=='al',death_time=st['death_time']))

R=Counter(); dct=Counter(); G=0
for p in particles:
    G += p['Z']
    for j in range(T):
        R[j]+=p['Y'][j]; dct[j]+=p['D'][j]
B=sum(int(latent[(r,0)]['Xi']>=1) for r in range(1,K+1))
N=sum(int(p['A']<T) for p in particles)
M=sum(int(p['A']<T and p['Xi']>=1) for p in particles)

# likelihood and MLE
def loglik(b):
    s=0.0
    for j in range(T):
        h=b/(b+j+1)
        s += dct[j]*math.log(h)+(R[j]-dct[j])*math.log(1-h)
    return s
bs=np.linspace(0.01,1.5,100000)
vals=np.array([loglik(float(b)) for b in bs])
idx=int(vals.argmax()); mle=float(bs[idx])

# Figure 1: trajectories, one page with 5 panels
fig, axes=plt.subplots(K,1,figsize=(9,10),sharex=True)
for r,ax in enumerate(axes,1):
    ps=[p for p in particles if p['r']==r]
    for p in ps:
        times=list(range(p['A'], T+1))
        pos=[]
        for t in times:
            rec=p['statepos'][t]
            pos.append(rec['pos'])
        line=ax.plot(times,pos,marker='o',label=f"x={p['x']}")[0]
        ax.scatter([p['A']],[p['x']],marker='o',s=70,facecolors='none',edgecolors=line.get_color(),linewidths=1.5)
        if p['Delta']==1:
            ax.scatter([p['death_time']],[p['X'][-1]],marker='x',s=70,color=line.get_color(),linewidths=2)
        else:
            ax.scatter([T],[p['statepos'][T]['pos']],marker='s',s=45,facecolors='none',edgecolors=line.get_color(),linewidths=1.5)
        ax.annotate(f"({r},{p['x']})",(p['A'],p['x']),xytext=(3,5),textcoords='offset points',fontsize=7)
    ax.set_ylabel(f"r={r}\nposition")
    ax.grid(True,alpha=.25)
    ax.set_xticks(range(T+1))
    ax.legend(loc='upper left',ncol=max(1,min(5,len(ps))),fontsize=7,frameon=False)
axes[-1].set_xlabel('calendar time $t$')
fig.suptitle(r'Labeled trajectories through $T=4$: circle = activation, x = observed death, square = right censoring',fontsize=12)
fig.tight_layout(rect=[0,0,1,.97])
fig.savefig(OUT/'trayectorias.pdf',bbox_inches='tight')
plt.close(fig)

# Figure 2 activation order
fig,ax=plt.subplots(figsize=(9,5.5))
for p in particles:
    ax.scatter(p['A'],p['n'],s=45)
    ax.annotate(f"({p['r']},{p['x']})",(p['A'],p['n']),xytext=(5,0),textcoords='offset points',va='center',fontsize=8)
ax.set_xlabel('activation time')
ax.set_ylabel(r'index $n$ in $Z_n^{(5)}$')
ax.set_xticks(range(T+1))
ax.set_yticks(range(1,len(particles)+1))
ax.grid(True,alpha=.25)
ax.set_title('Global activation order across the five realizations')
fig.tight_layout()
fig.savefig(OUT/'orden_activacion.pdf',bbox_inches='tight')
plt.close(fig)

# Figure 3 censoring bars
fig,ax=plt.subplots(figsize=(9,7))
labels=[]
for y,p in enumerate(particles):
    labels.append(f"({p['r']},{p['x']})")
    start=p['A']; end=(p['death_time'] if p['Delta'] else T)
    ax.hlines(y,start,end,linewidth=2)
    ax.scatter(start,y,marker='o',s=35)
    if p['Delta']:
        ax.scatter(end,y,marker='x',s=55,linewidths=2)
    else:
        ax.scatter(end,y,marker='s',s=35,facecolors='none')
ax.set_yticks(range(len(particles))); ax.set_yticklabels(labels,fontsize=8)
ax.invert_yaxis(); ax.set_xticks(range(T+1)); ax.grid(True,axis='x',alpha=.25)
ax.set_xlabel('calendar time')
ax.set_title('Individual observation windows: x = observed death; square = right censoring at $T$')
fig.tight_layout()
fig.savefig(OUT/'censura.pdf',bbox_inches='tight')
plt.close(fig)

# Figure 4 age counts
ages=np.arange(T)
fig,ax=plt.subplots(figsize=(8,4.8))
width=.36
ax.bar(ages-width/2,[dct[j] for j in ages],width,label='deaths $d(j)$')
ax.bar(ages+width/2,[R[j]-dct[j] for j in ages],width,label='jumps $R(j)-d(j)$')
ax.set_xticks(ages); ax.set_xlabel('individual age $j$'); ax.set_ylabel('number of visible rows')
ax.set_title('Descomposition de las filas observadas por edad')
ax.legend(); ax.grid(True,axis='y',alpha=.25)
fig.tight_layout(); fig.savefig(OUT/'conteos_edad.pdf',bbox_inches='tight'); plt.close(fig)

# Figure 5 likelihood curve normalized
fig,ax=plt.subplots(figsize=(8,4.8))
bs2=np.linspace(.03,1.2,1000); lv=np.array([loglik(float(b)) for b in bs2]); rel=np.exp(lv-lv.max())
ax.plot(bs2,rel)
ax.axvline(beta,linestyle='--',label=r'generating value $\beta_0=0.45$')
ax.axvline(mle,linestyle=':',label=fr'MLE $\widehat{{\beta}}={mle:.4f}$')
ax.set_xlabel(r'$\beta$'); ax.set_ylabel('relative likelihood')
ax.set_title('Likelihood function for the simulated record')
ax.legend(); ax.grid(True,alpha=.25)
fig.tight_layout(); fig.savefig(OUT/'likelihood.pdf',bbox_inches='tight'); plt.close(fig)

res=dict(K=K,T=T,beta=beta,seed=seed,particles=particles,R={str(j):R[j] for j in range(T)},d={str(j):dct[j] for j in range(T)},G=G,B=B,N=N,M=M,mle=mle)
(OUT/'example.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
print(json.dumps({k:res[k] for k in ['K','T','beta','seed','R','d','G','B','N','M','mle']},indent=2))
