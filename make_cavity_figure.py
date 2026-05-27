import matplotlib as mpl
mpl.rcParams["pdf.fonttype"]=42
mpl.rcParams["ps.fonttype"]=42
mpl.rcParams["font.family"]="sans-serif"
import numpy as np, matplotlib.pyplot as plt
from scipy.special import jn_zeros, j0
plt.rcParams.update({'font.size':10,'mathtext.fontset':'cm'})
hbar_c=0.1973269804; R0=0.81; invR0=hbar_c/R0
x0=jn_zeros(0,1)[0]; x1=jn_zeros(1,1)[0]; x2=jn_zeros(2,1)[0]
fig,(axL,axR)=plt.subplots(1,2,figsize=(7.2,3.5),gridspec_kw={'width_ratios':[1,1.3]})
# left: filled-contour of J0(x0 r)^2 intensity on the disk
N=300; xx=np.linspace(-1,1,N); X,Y=np.meshgrid(xx,xx); Rr=np.sqrt(X**2+Y**2)
Z=np.where(Rr<=1, j0(x0*Rr)**2, np.nan)
im=axL.imshow(Z,extent=[-1,1,-1,1],origin='lower',cmap='Blues',vmin=0,vmax=1)
theta=np.linspace(0,2*np.pi,200)
axL.plot(np.cos(theta),np.sin(theta),'k-',lw=1.8)
axL.annotate(r'$R_0=0.81$ fm',(0,-1.34),ha='center',fontsize=10)
axL.annotate(r'$|{\rm field}|^2\!\propto\!J_0^2(x_{0,1}r/R_0)$',(0,1.18),ha='center',fontsize=8.5,color='#16407a')
axL.set_xlim(-1.45,1.45); axL.set_ylim(-1.5,1.4); axL.set_aspect('equal'); axL.axis('off')
# right: levels
for name,xv,col in [('$\\Sigma$ ($x_{0,1}$)',x0,'#1f4e9c'),
                    ('$\\Pi_u$ ($x_{1,1}$)',x1,'#27ae60'),
                    ('$\\Delta_g$ ($x_{2,1}$)',x2,'#c0392b')]:
    E=xv*invR0; axR.hlines(E,0.12,0.78,color=col,lw=2.6)
    axR.text(0.84,E,f'{name}: {E:.2f} GeV',va='center',fontsize=9.2,color=col)
axR.set_ylim(0.4,1.45); axR.set_xlim(0,2.5); axR.set_ylabel('gluonic gap  [GeV]'); axR.set_xticks([])
for s in ['top','right','bottom']: axR.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig('Fig1.pdf', bbox_inches='tight'); plt.close()
print("ok")
