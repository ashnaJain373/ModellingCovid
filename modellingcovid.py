import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys

# Total population, N.
N = 136000000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 561591, 8973373
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.065, 1/14 
# A grid of time points (in days)
#t = np.linspace(0,90, 90)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.065, 1/14 

# Initial conditions vector
y0 = S0, I0, R0
Svaccine=S0
Ivaccine=I0
Rvaccine=R0
yvaccine = Svaccine, Ivaccine, Rvaccine

SvaccineDouble=S0
IvaccineDouble=I0
RvaccineDouble=R0
yvaccineDouble = SvaccineDouble, IvaccineDouble, RvaccineDouble

# Initial conditions vector
y0 = S0, I0, R0
SFutureSS=S0
IFutureSS=I0
RFutureSS=R0
yFutureSS = SFutureSS, IFutureSS, RFutureSS

SFutureLD=S0
IFutureLD=I0
RFutureLD=R0
yFutureLD = SFutureLD, IFutureLD, RFutureLD

SElection=S0
IElection=I0
RElection=R0
yElection = SElection, IElection, RElection

SKumbh=S0
IKumbh=I0
RKumbh=R0
yKumbh = SKumbh, IKumbh, RKumbh


Splot=[]
Iplot=[]
Rplot=[]
Totalplot=[]
Totalplotvaccine=[]
TotalplotvaccineDouble=[]
TotalplotElection=[]
TotalplotKumbh=[]
TotalplotFutureSS=[]
TotalplotFutureLD=[]

options = str(sys.argv)
option = options[len(options)-3]

FutureLockDown = False
FutureVaccine = False
numMonths = 7;
Election=True
Kumbh=True
startKumbh=False
startElection=False

# 5 months current data till May 1
# Elections + Kumbh
if option == 'A':
        Election=True
        Kumbh=True
        numMonths = 5
# 5 months current data till May 1
# NoElections + Kumbh
elif option == 'B':
        numMonths = 5
        Election=False
        Kumbh=True
# 5 months current data till May 1
# NoElections + Kumbh
elif option == 'C':
        numMonths = 5
        Kumbh=False
# 5 months current data till May 1
# NoElections + NoKumbh
elif option == 'D':
        numMonths = 5
        Kumbh=False
        Election=False
elif option == 'E':
        FutureLockDown = True
elif option == 'F':
        FutureVaccine = True
    
daysMonths = 30
days = numMonths * daysMonths
StartVaccination = False
StartFuture = False

for time in range(days):
    # A grid of time points (in days)
    t = np.linspace(0,2,2)

    # election at 90 days
    # We need to model the case if elections did not happen
    
    if time == 3 * daysMonths:
        startElection=True
        beta = 0.095
        
        if not Election:
            yElection=y0
            
            # Need 2 copies
            for elem in Totalplot:
                TotalplotElection.append(elem)

    # kumbh at 120 days        
    elif time == 4 * daysMonths:
        startKumbh=True
        beta = 0.135
        if not Kumbh:
            yKumbh=y0

            if Election:
                # Need 2 copies
                for elem in Totalplot:
                    TotalplotKumbh.append(elem)
            else:
                # Need 2 copies
                for elem in TotalplotElection:
                    TotalplotKumbh.append(elem)
                

    if FutureLockDown:
        if time ==  5 * daysMonths:
            yFutureLD=y0
            yFutureSS=y0

            # assume normal - no lockdown; no super-spreader
            beta = 0.065
            StartFuture = True

            # Need 2 copies
            for elem in Totalplot:
                TotalplotFutureSS.append(elem)
                TotalplotFutureLD.append(elem)              
            
    elif FutureVaccine:
        if time ==  5 * daysMonths:
            yvaccine=y0
            yvaccineDouble=y0

            # assume normal - no lockdown; no super-spreader
            beta = 0.065
            StartVaccination = True

            # Need 2 copies
            for elem in Totalplot:
              Totalplotvaccine.append(elem)
              TotalplotvaccineDouble.append(elem)              
        
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma))
    S, I, R = ret.T

    S0=S[len(S)-1]
    I0=I[len(I)-1]
    R0=R[len(R)-1]

    y0 = S0, I0, R0
  
    Splot.append(S0)
    Iplot.append(I0)
    Rplot.append(R0)
    Totalplot.append(I0+R0)

    if not Election and startElection:
        if not startKumbh:
            beta = 0.065
        elif not Kumbh:
            beta = 0.065

        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, yElection, t, args=(N, beta, gamma))
        S, I, R = ret.T

        SElection=S[len(S)-1]
        IElection=I[len(I)-1]
        RElection=R[len(R)-1]
        
        yElection = SElection, IElection, RElection
        TotalplotElection.append(IElection+RElection)

        #restore beta
        if not startKumbh:
            beta = 0.095
        elif not Kumbh:
            beta = 0.135

        if not Kumbh:
            continue

    if not Kumbh and startKumbh:
        beta = 0.065

        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, yKumbh, t, args=(N, beta, gamma))
        S, I, R = ret.T

        SKumbh=S[len(S)-1]
        IKumbh=I[len(I)-1]
        RKumbh=R[len(R)-1]
        
        yKumbh = SKumbh, IKumbh, RKumbh
        TotalplotKumbh.append(IKumbh+RKumbh)

        #restore beta
        beta = 0.135

    if StartFuture:
        #Lockdown
        beta = 0.035

        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, yFutureLD, t, args=(N, beta, gamma))
        S, I, R = ret.T

        SFutureLD=S[len(S)-1]
        IFutureLD=I[len(I)-1]
        RFutureLD=R[len(R)-1]
        
        yFutureLD = SFutureLD, IFutureLD, RFutureLD
        TotalplotFutureLD.append(IFutureLD+RFutureLD)

        #Another Super Spreader
        beta = 0.135

        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, yFutureSS, t, args=(N, beta, gamma))
        S, I, R = ret.T

        SFutureSS=S[len(S)-1]
        IFutureSS=I[len(I)-1]
        RFutureSS=R[len(R)-1]
        
        yFutureSS = SFutureSS, IFutureSS, RFutureSS
        TotalplotFutureSS.append(IFutureSS+RFutureSS)
        
        #restore beta
        beta = 0.065
        
    if StartVaccination:
        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, yvaccine, t, args=(N, beta, gamma))
        S, I, R = ret.T

        Svaccine=S[len(S)-1]
        Ivaccine=I[len(I)-1]
        Rvaccine=R[len(R)-1]
        Svaccine=Svaccine-260000
        
        yvaccine = Svaccine, Ivaccine, Rvaccine
        Totalplotvaccine.append(Ivaccine+Rvaccine)

        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, yvaccineDouble, t, args=(N, beta, gamma))
        S, I, R = ret.T
        
        SvaccineDouble=S[len(S)-1]
        IvaccineDouble=I[len(I)-1]
        RvaccineDouble=R[len(R)-1]
        SvaccineDouble=SvaccineDouble-(260000*3)
        
        yvaccineDouble = SvaccineDouble, IvaccineDouble, RvaccineDouble
        TotalplotvaccineDouble.append(IvaccineDouble+RvaccineDouble)
        
    
#end of for loop

            
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#num=10000000

color = 'r'
if numMonths == 5:
    t = np.linspace(0, 90, 90)
    ax.plot(t, Totalplot[0:90], color, alpha=0.5, lw=2, label='Confirmed Cases')

    if Election and Kumbh:
        color = 'b'
        t = np.linspace(90, 120, 30)
        ax.plot(t, Totalplot[90:120], color, alpha=0.5, lw=2, label='')

        color = 'g'        
        t = np.linspace(120, 150, 30)
        ax.plot(t, Totalplot[120:150], color, alpha=0.5, lw=2, label='')
    else:
        color = 'r'
        t = np.linspace(90, 120, 30)
        ax.plot(t, Totalplot[90:120], color, alpha=0.5, lw=2, label='')

        color = 'r'        
        t = np.linspace(120, 150, 30)
        ax.plot(t, Totalplot[120:150], color, alpha=0.5, lw=2, label='')
    
    if not Election and Kumbh:
        t = np.linspace(90, 120, 30)
        ax.plot(t, TotalplotElection[90:120], 'b', alpha=0.5, lw=2, label='No Election')
        
        t = np.linspace(120, 150, 30)
        ax.plot(t, TotalplotElection[120:150], 'g', alpha=0.5, lw=2, label='')
        
    if Election and not Kumbh:
        t = np.linspace(90, 120, 30)
        ax.plot(t, TotalplotKumbh[90:120], 'b', alpha=0.5, lw=2, label='')
        
        t = np.linspace(120, 150, 30)
        ax.plot(t, TotalplotKumbh[120:150], 'g', alpha=0.5, lw=2, label='No Kumbh')

    if not Election and not Kumbh:
        t = np.linspace(90, 120, 30)
        ax.plot(t, TotalplotElection[90:120], 'g', alpha=0.5, lw=2, label='No Election and No Kumbh')

        t = np.linspace(120, 150, 30)
        ax.plot(t, TotalplotElection[120:150], 'g', alpha=0.5, lw=2, label='')
  
elif numMonths == 7:
    t = np.linspace(0, 150, 150)
    ax.plot(t, Totalplot[0:150], 'r', alpha=0.5, lw=2, label='Confirmed Cases Current')

    if FutureVaccine:
        t = np.linspace(150, 250, 60)
        ax.plot(t, Totalplot[150:210], 'y', alpha=0.5, lw=2, label='Confirmed Cases No Vaccination')
        ax.plot(t, Totalplotvaccine[150:210], 'g', alpha=0.5, lw=2, label='Confirmed Cases Current Vaccination')
        ax.plot(t, TotalplotvaccineDouble[150:210], 'b', alpha=0.5, lw=2, label='Confirmed Triple Vaccination')
        
    if FutureLockDown:
        t = np.linspace(150, 250, 60)
        ax.plot(t, Totalplot[150:210], 'r', alpha=0.5, lw=2, label='')
        ax.plot(t, TotalplotFutureLD[150:210], 'g', alpha=0.5, lw=2, label='Confirmed Cases if Lockdown')
        ax.plot(t, TotalplotFutureSS[150:210], 'b', alpha=0.5, lw=2, label='Confirmed Cases if Super Spreader')        

  
ax.set_xlabel('Time /days from dec 2020')
ax.set_ylabel('Number (10000000s)')
if FutureLockDown:
    ax.set_ylim(0,50000000)
else:
    ax.set_ylim(0,30000000)    
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)

plt.show()
