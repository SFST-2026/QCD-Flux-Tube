import matplotlib as mpl
mpl.rcParams["pdf.fonttype"]=42
mpl.rcParams["ps.fonttype"]=42
mpl.rcParams["font.family"]="sans-serif"
import numpy as np, matplotlib.pyplot as plt
import matplotlib.patches as mp
from scipy.special import jn_zeros
plt.rcParams.update({'font.size':10,'axes.linewidth':0.8,'mathtext.fontset':'cm'})
hbar_c=0.1973269804; R0=0.81; invR0=hbar_c/R0
x0=jn_zeros(0,1)[0]; x1=jn_zeros(1,1)[0]; x2=jn_zeros(2,1)[0]
B=2.58; k_perp=x0*invR0
def dbreak(mpi,sig):
    mq=(mpi/1000)**2/(2*B); meff=np.sqrt(mq**2+k_perp**2); return 2*meff/sig*hbar_c

# ============ FIG 1: d_break(m_pi) with TWO data points ============
fig,ax=plt.subplots(figsize=(7.0,4.4))
mpi=np.linspace(0,750,300)
d1=np.array([dbreak(m,np.pi*invR0**2) for m in mpi])
d2=np.array([dbreak(m,0.19) for m in mpi])
ax.fill_between(mpi,1.00,1.50,color='0.94',zorder=0)
ax.fill_between(mpi,d2,d1,color='#cfe0f5',alpha=0.7,zorder=1,label=r'model band, $\sigma\in[\pi/R_0^2,\,0.19]$')
ax.plot(mpi,d1,'-',color='#1f4e9c',lw=2.3,zorder=3,label=r'model, $\sigma=\pi/R_0^2$')
ax.axhline(1.53*R0,color='#1f4e9c',ls=':',lw=1.1,alpha=.7,zorder=2)
ax.annotate(r'chiral plateau $\dfrac{2x_{0,1}}{\pi}R_0=1.24$ fm',(30,1.252),
            fontsize=9.5,color='#1f4e9c',va='bottom')
ax.errorbar([640],[1.25],yerr=[0.01],fmt='o',ms=8,color='#c0392b',capsize=4,
            mec='k',mew=0.5,zorder=5,label=r'Bali/SESAM $N_f{=}2$, $m_\pi\!\sim\!640$')
ax.errorbar([200],[1.13],yerr=[0.14],fmt='D',ms=7,color='#e67e22',capsize=4,
            mec='k',mew=0.5,zorder=5,label=r'Bulava et al.\ phys.\ point')
ax.set_xlabel(r'$m_\pi$  [MeV]'); ax.set_ylabel(r'$d_{\mathrm{break}}$  [fm]')
ax.set_xlim(-15,760); ax.set_ylim(1.00,1.42)
ax.legend(fontsize=8.2,loc='lower right',framealpha=.96,edgecolor='0.8')
ax.grid(alpha=.22,lw=0.6)
fig.tight_layout(); fig.savefig('Fig2.pdf', bbox_inches='tight'); plt.close()

# ============ FIG 2: cavity mode tower (the "pretty" schematic) ============
fig,(axL,axR)=plt.subplots(1,2,figsize=(7.2,3.6),gridspec_kw={'width_ratios':[1,1.25]})
# left: cylinder cross-section with first few transverse modes (J0,J1,J2 radial)
theta=np.linspace(0,2*np.pi,200)
axL.plot(np.cos(theta),np.sin(theta),'k-',lw=1.5)
axL.fill(np.cos(theta),np.sin(theta),color='#eaf1fb',zorder=0)
# illustrate field intensity ~ J0(x0 r)^2 as shading rings
from scipy.special import j0
rr=np.linspace(0,1,60)
for i in range(len(rr)-1):
    val=j0(x0*rr[i])**2
    circ=plt.Circle((0,0),rr[i],color=plt.cm.Blues(0.25+0.6*val),fill=True,zorder=0)
    axL.add_patch(circ)
axL.plot(np.cos(theta),np.sin(theta),'k-',lw=1.6,zorder=4)
axL.annotate(r'$R_0=0.81$ fm',(0,-1.32),ha='center',fontsize=10)
axL.annotate('flux-tube\ncross-section',(0,0),ha='center',va='center',fontsize=8.5,color='#16407a')
axL.set_xlim(-1.5,1.5); axL.set_ylim(-1.55,1.4); axL.set_aspect('equal'); axL.axis('off')
# right: energy-level diagram of x_{L,1}/R0
levels=[('$\\Sigma$ ($x_{0,1}$)',x0,'#1f4e9c'),
        ('$\\Pi_u$ ($x_{1,1}$)',x1,'#27ae60'),
        ('$\\Delta_g$ ($x_{2,1}$)',x2,'#c0392b')]
for name,xv,col in levels:
    E=xv*invR0
    axR.hlines(E,0.15,0.85,color=col,lw=2.4)
    axR.text(0.9,E,f'{name}: {E:.2f} GeV',va='center',fontsize=9,color=col)
axR.set_ylim(0.4,1.45); axR.set_xlim(0,2.4)
axR.set_ylabel('gluonic gap  [GeV]'); axR.set_xticks([])
axR.spines['top'].set_visible(False); axR.spines['right'].set_visible(False)
axR.spines['bottom'].set_visible(False)
fig.tight_layout(); fig.savefig('Fig1.pdf', bbox_inches='tight'); plt.close()
print("both figures saved")
