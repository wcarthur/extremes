"""
Evaluating non-stationary hazard function values

"""


import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import genpareto
import scipy.integrate as integrate



def hft(t, x, shp, loc, scale, M=2, d=100):
    """
    Nonstationary hazard function
    """
    k = scale * np.exp(np.log(M) * t / d)
    hf = 1 - genpareto.cdf(x, shp, loc=loc, scale=k)
    return hf

def chf(tmin, tmax, x, shp, loc, scale, M=2, d=100):
    """
    Cumulative hazard function
    """
    return integrate.quad(hft, tmin, tmax, args=(x, shp, loc, scale, M, d))

def ttf(tmin, tmax, x, shp, loc, scale, M, d):
    """
    Distribution of time to failure
    """
    Ck = 1 / np.sqrt(2 * shp + 1)
    p = hft(tmax, x, shp, loc, scale, M, d)
    ex = 2 * Ck**2 / (1 - Ck**2)
    def integrand(t, p, Ck, M, d):
        ex = 2 * Ck**2 / (1 - Ck**2)
        return np.power(1 - (1 - np.power(p, 1/ex))/np.power(M, t / d), 2 * Ck**2 / (1 - Ck**2))
    
    val = integrate.quad(integrand, tmin, tmax, args=(p, Ck, M, d))
    res = integrand(tmax, p, Ck, M, d) * np.exp(-val[0])
    return res

def sf(tmin, tmax, x, shp, loc, scale, M=2, d=100):
    """
    Survival function
    """
    Ck = 1 / np.sqrt(2 * shp + 1)
    p = hft(tmax, x, shp, loc, scale, M, d)
    ex = 2 * Ck**2 / (1 - Ck**2)
    def integrand(t, p, Ck, M, d):
        ex = 2 * Ck**2 / (1 - Ck**2)
        return np.power(1 - (1 - np.power(p, 1/ex))/np.power(M, t / d), 2 * Ck**2 / (1 - Ck**2))
    
    val = integrate.quad(integrand, tmin, tmax, args=(p, Ck, M, d))
    res = np.exp(-val[0])
    return res

def quantile(p, shp, scale):
    return (scale/shp) * (1 - np.power(p, shp))


scale = 10
beta = 0.001
loc = 0
shp = 0.0388888

#x = genpareto.ppf(0.999, shp, loc=loc, scale=scale)

#p0 = 0.002
#x = (scale/shp) * (1 - np.power(p0, shp))
#print(x)

tmin=0
tmax=1000

x = 50
p0 = hft(0, x, shp, loc, scale)
p0 = 0.002
x = quantile(p0, shp, scale)
Cx = 1 / np.sqrt(2 * shp + 1)
print(f"p0={p0:.6f}; Cx = {Cx:.4f}; x = {x:.2f}")




ep = [quantile(p, shp, scale) for p in np.logspace(-4, -1, endpoint=True, num=31)]
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.semilogy(ep, np.logspace(-4, -1, num=31, endpoint=True))
ax.set_ylabel("Exceedance probability")
ax.set_xlabel("Magnitude")
ax.vlines(x=x, ymin=0, ymax=p0, ls='--')
ax.grid(which='both')
plt.show()

# Hazard function for a fixed time:

xx = np.arange(0, 101)
p = genpareto.pdf(xx, shp, loc=loc, scale=scale)
c = genpareto.cdf(xx, shp, loc=loc, scale=scale)
fig, axes = plt.subplots(2, 1, figsize=(12, 8))
ax = axes.flatten()
ax[0].plot(xx, p)
ax[1].plot(xx, c)
ax[0].set_title("PDF")
ax[1].set_title("CDF")
ax[0].grid(which='both')
ax[1].grid(which='both')

plt.show()


tt = np.arange(tmax+1)

hh1 = [hft(t, x, shp, loc, scale, M=1.0, d=100) for t in tt]
hh2 = [hft(t, x, shp, loc, scale, M=1.1, d=100) for t in tt]
hh3 = [hft(t, x, shp, loc, scale, M=1.5, d=100) for t in tt]
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.semilogx(tt, hh1, label="M=1.0")
ax.semilogx(tt, hh2, label="M=1.1")
ax.semilogx(tt, hh3, label="M=1.5")
ax.legend()
ax.grid(which='both')
ax.set_title("Hazard function")
plt.show()

# Cumulative hazard function
H1 = [chf(tmin, t, x, shp, loc, scale, M=1, d=100)[0] for t in tt]
H2 = [chf(tmin, t, x, shp, loc, scale, M=1.1, d=100)[0] for t in tt]
H3 = [chf(tmin, t, x, shp, loc, scale, M=1.5, d=100)[0] for t in tt]

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.semilogx(tt, H1, label="M=1")
ax.semilogx(tt, H2, color='r', label="M=1.1")
ax.semilogx(tt, H3, color='k', label="M=1.5")

#ax.set_ylim(0, 20)
ax.legend()
ax.grid(which='both')
ax.set_title("Cumulative hazard function")
plt.show()

# PDF of time to exceedance (time-to-failure)
TTF1 = [ttf(0, t, x, shp, loc, scale, M=1.0, d=100) for t in tt]
TTF2 = [ttf(0, t, x, shp, loc, scale, M=1.02, d=100) for t in tt]
TTF3 = [ttf(0, t, x, shp, loc, scale, M=1.1, d=100) for t in tt]
TTF4 = [ttf(0, t, x, shp, loc, scale, M=1.25, d=100) for t in tt]

# Find median time to exceedance
TTF1med = np.argmin(np.abs(np.cumsum(TTF1)-0.5))
TTF2med = np.argmin(np.abs(np.cumsum(TTF2)-0.5))
TTF3med = np.argmin(np.abs(np.cumsum(TTF3)-0.5))
TTF4med = np.argmin(np.abs(np.cumsum(TTF4)-0.5))

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.plot(tt, TTF1, label="1.0 (stat)", color='r')
ax.vlines(TTF1med, ymin=0, ymax=TTF1[TTF1med], ls='--', color='r')
ax.plot(tt, TTF2, label="1.02", color='g')
ax.vlines(TTF2med, ymin=0, ymax=TTF2[TTF2med], ls='--', color='g')

ax.plot(tt, TTF3, label="1.10", color='b')
ax.vlines(TTF3med, ymin=0, ymax=TTF3[TTF3med], ls='--', color='b')

ax.plot(tt, TTF4, label="1.25", color='k')
ax.vlines(TTF4med, ymin=0, ymax=TTF4[TTF4med], ls='--', color='k')

ax.legend(title=f"Cx={np.round(Cx, 2)}, M:")
ax.set_xlim(0, 500)
ax.grid(which='both')
ax.set_title(f"PDF of time to failure (p0={np.round(p0, 4)})")
plt.show()