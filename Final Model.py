import numpy as np
from scipy.integrate import odeint
import pandas as pd

# Gather inputs from the user
N = int(input("Please enter the total population: "))
I0 = int(input("Please enter the initial number of infected: "))
R0 = int(input("Please enter the initial number of recovered: "))
E0 = int(input("Please enter the initial number of exposed: "))

S0 = N - E0 - I0 - R0

# Contact rate, recovery rate, exposed rate, and direct recovery factor
beta = float(input("Please enter the contact rate: "))
gamma = float(input("Please enter the recovery rate: ")) 
sigma = float(input("Please enter the exposed rate: ")) 
meu = float(input("Please enter the direct recovery factor: ")) 

n = int(input("Please enter the number of days: "))
t = np.linspace(0, n, num=n + 1)  # +1 to include the endpoint

# The SEIR model differential equations.
def SEIR_model(y, t, N, beta, gamma, sigma, meu):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma * E - meu * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I + meu * E

    return dSdt, dEdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, E0, I0, R0

# Integrate the SEIR equations over the time grid
sol = odeint(SEIR_model, y0, t, args=(N, beta, gamma, sigma, meu))
S, E, I, R = sol.T

# Convert to integer
S = S.astype(int)
E = E.astype(int)
I = I.astype(int)
R = R.astype(int)

# Create a DataFrame and save to CSV
df1 = pd.DataFrame(data={'Suspected': S, 'Exposed': E, 'Infected': I, 'Recovered': R})
df1.index = df1.index + 1  # Start index at 1
df1.to_csv('Covid.csv', index=True)
