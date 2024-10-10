import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

base = importr('base')
utils = importr('utils')
evmix = importr('evmix')
stats = importr('stats')
numpy2ri.activate()

def ppoints(n, a=0.375):
    if isinstance(n, int):
        m = np.arange(1, n + 1)
        if n > 10:
            a = 0.5
    else:
        m = np.arange(1, len(n) + 1)
        if len(n) > 10:
            a = 0.5
    pp = (m - a)/(m[-1] + (1 - a) - a)
    return pp

df = pd.read_csv(r"X:\georisk\HaRIA_B_Wind\data\raw\from_bom\2019\Daily\DC02D_Data_003003_999999999632559.txt", skipinitialspace=True)
gust = df[df['Quality of maximum gust speed']=='Y']['Speed of maximum wind gust in m/s'].dropna().values

xx = np.arange(0, 60, 0.01)

fit = evmix.fgammagpd(gust)
dfit = evmix.dgammagpd(xx, fit.rx2('gshape')[0], fit.rx2('gscale')[0], fit.rx2('u')[0], fit.rx2('sigmau')[0], fit.rx2('xi')[0])
idxfit = np.argmin(np.abs(xx - fit.rx2('u')[0]))
fit2 = evmix.fgammagpd(gust, phiu=False)
dfit2 = evmix.dgammagpd(xx, fit2.rx2('gshape')[0], fit2.rx2('gscale')[0], fit2.rx2('u')[0], fit2.rx2('sigmau')[0], fit2.rx2('xi')[0], fit2.rx2('phiu'))
idxfit2 = np.argmin(np.abs(xx - fit2.rx2('u')[0]))

fitu = evmix.fgammagpd(gust, useq=np.linspace(np.median(gust), gust.max(), 20))
dfitu = evmix.dgammagpd(xx, fitu.rx2('gshape')[0], fitu.rx2('gscale')[0], fitu.rx2('u')[0], fitu.rx2('sigmau')[0], fitu.rx2('xi')[0])
idxfitu = np.argmin(np.abs(xx - fitu.rx2('u')[0]))

fitfix = evmix.fgammagpd(gust, useq=np.linspace(np.median(gust), gust.max(), 20), fixedu=True)
dfitfix = evmix.dgammagpd(xx, fitfix.rx2('gshape')[0], fitfix.rx2('gscale')[0], fitfix.rx2('u')[0], fitfix.rx2('sigmau')[0], fitfix.rx2('xi')[0])
idxfitfix = np.argmin(np.abs(xx - fitfix.rx2('u')[0]))

print("Bulk tail fraction")
print(f"u: {fit.rx2('u')[0]:.3f}, sigma: {fit.rx2('sigmau')[0]:.4f}, xi: {fit.rx2('xi')[0]:.4f}")
print("Parameterised tail fraction")
print(f"u: {fit2.rx2('u')[0]:.3f}, sigma: {fit2.rx2('sigmau')[0]:.4f}, xi: {fit2.rx2('xi')[0]:.4f}")

print("Prof. lik. for initial value")
print(f"u: {fitu.rx2('u')[0]:.3f}, sigma: {fitu.rx2('sigmau')[0]:.4f}, xi: {fitu.rx2('xi')[0]:.4f}")
print("Prof. lik. for fixed threshold")
print(f"u: {fitfix.rx2('u')[0]:.3f}, sigma: {fitfix.rx2('sigmau')[0]:.4f}, xi: {fitfix.rx2('xi')[0]:.4f}")

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].hist(gust, bins=100, density=True, fc='white', ec='k')

ax[0].plot(xx, dfit, color='red', label=f"Bulk tail fraction: $u={{{fit.rx2('u')[0]:.1f}}}$")
ax[0].vlines(fit.rx2('u')[0], 0, dfit[idxfit], color='red')
ax[0].plot(xx, dfit2, color='blue', label=f"Parameterised tail fraction: $u={{{fit2.rx2('u')[0]:.1f}}}$")
ax[0].vlines(fit2.rx2('u')[0], 0, dfit2[idxfit2], color='blue')
ax[0].legend()

ax[1].hist(gust, bins=100, density=True, fc='white', ec='k')

ax[1].plot(xx, dfit, color='red', label="Default initial value (90% quantile)")
ax[1].vlines(fit.rx2('u')[0], 0, dfit[idxfit], color='red')

ax[1].plot(xx, dfitu, color='purple', label="Prof. lik. for initial value")
ax[1].vlines(fitu.rx2('u')[0], 0, dfitu[idxfitu], color='purple')

ax[1].plot(xx, dfitfix, color='green', label="Prof. lik. for fixed threshold")
ax[1].vlines(fitfix.rx2('u')[0], 0, dfitfix[idxfitfix], color='green')

ax[1].legend()

plt.show()

n = len(gust)
empprob = ppoints(n)
transempprob = -1/np.log(empprob)
#transempprob = 1/(1 - empprob)
minemppower = -np.log10(1 - 1/n/100)
maxemppower = np.ceil(np.log10(np.max(transempprob))) + 1
theprob = 1 - np.power(10, -np.arange(minemppower, maxemppower+0.05, 0.05))
theprob = np.sort(np.concatenate((theprob, 1 - theprob)))
transtheprob = -1/np.log(theprob)
#transtheprob = 1/(1-theprob)

gshape = fit.rx2('gshape')[0]
gscale = fit.rx2('gscale')[0]
u = fit.rx2('u')[0]
sigmau = fit.rx2('sigmau')[0]
xi = fit.rx2('xi')[0]
phiu = fit.rx2('phiu')[0]
thequant = evmix.qgammagpd(theprob, gshape, gscale, u, sigmau, xi, phiu)

thequant2 = evmix.qgammagpd(theprob, fit2.rx2('gshape')[0], fit2.rx2('gscale')[0], fit2.rx2('u')[0], fit2.rx2('sigmau')[0], fit2.rx2('xi')[0], fit.rx2('phiu')[0])
thequantu = evmix.qgammagpd(theprob, fitu.rx2('gshape')[0], fitu.rx2('gscale')[0], fitu.rx2('u')[0], fitu.rx2('sigmau')[0], fitu.rx2('xi')[0], fitu.rx2('phiu')[0])
thequantfix = evmix.qgammagpd(theprob, fitfix.rx2('gshape')[0], fitfix.rx2('gscale')[0], fitfix.rx2('u')[0], fitfix.rx2('sigmau')[0], fitfix.rx2('xi')[0], fitfix.rx2('phiu')[0])

if phiu < 0.1:
    xlims = (1/0.1/365.25, np.power(10, maxemppower)/365.25)
else:
    xlims = (1/phiu/365.25, np.power(10, maxemppower)/365.25)

ylims = (min(np.quantile(gust, 0.8), u), max(np.concatenate((gust, thequant))))
N=1000
alpha=0.05
simdata = np.zeros((N, n))
for i in range(N):
    simdata[i, :] = np.sort(evmix.rgammagpd(n, gshape, gscale, u, sigmau, xi, phiu))

ci = np.quantile(simdata, q=(alpha/2, 1-alpha/2), axis=0)

fig, ax = plt.subplots(1, 1)
ax.scatter(transempprob/365.25, np.sort(gust), color='k', marker='x', label='Data')
ax.plot(transtheprob/365.25, thequant, color='r', label='Fitted gamma-gpd')
ax.plot(transtheprob/365.25, thequant2, color='b', ls='--', label='Parameterised tail fraction')
ax.plot(transtheprob/365.25, thequantu, color='purple', ls='--', label='Prof. lik. for initial value')
ax.plot(transtheprob/365.25, thequantfix, color='green', ls='--', label='Prof. lik. for fixed threshold')
ax.plot(transempprob/365.25, ci.T, color='k', linewidth=1, linestyle='--')
ax.axhline(u, ls='--', lw=1, color='0.5')
ax.axvline(1/(365.25*phiu), ls='--', lw=1, color='0.5')
ax.grid(which='major', linestyle='-')
ax.grid(which='minor', linestyle='--', linewidth=1)
ax.set_xscale('log')
ax.set_xlabel("ARI")
ax.set_ylabel("Return level")
ax.set_xlim(xlims)
ax.set_ylim(ylims)
plt.show()