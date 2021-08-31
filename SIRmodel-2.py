import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 136000000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 561591, 8973373
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.065, 1/14 
# A grid of time points (in days)
t = np.linspace(0,150, 150)


##elections
beta2, gamma = 0.095, 1/14
t2 = np.linspace(90, 120, 120)


##kumbh mela
beta3,gamma=0.125,1/14
t3=np.linspace(120,150,150)



# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
num=10000000
#ax.plot(t, S/(num), 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, (I+R)/(num), 'r', alpha=0.5, lw=2, label='Confirmed')
###ax.plot(t, (I+R)/(num), 'g', alpha=0.5, lw=2, label='Confirmed')
### ax.plot(t, R/(num), 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days from dec 2020')
ax.set_ylabel('Number (10000000s)')
ax.set_ylim(0,2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)

# Initial number of infected and recovered individuals, I0 and R0.
#I0, R0 = 561591, 10973373
I0, R0 = 580001,10984373

y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t2, args=(N, beta2, gamma))
S, I, R = ret.T

ax1=ax
ax1.plot(t2, (I+R)/(num), 'b', alpha=0.5, lw=2)
legend = ax1.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax1.spines[spine].set_visible(False)


I0, R0 = (961591, 12525039)

y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t3, args=(N, beta3, gamma))
S, I, R = ret.T

ax2=ax
ax2.plot(t3, (I+R)/(num), 'g', alpha=0.5, lw=2)
legend = ax2.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax1.spines[spine].set_visible(False)
plt.title("If elections and kumbh mela(superspreader events) hadn't happened")
plt.show()








