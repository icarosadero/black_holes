import plotly.offline as py
import plotly.graph_objs as go
import numpy as np
#import pandas as pd
from scipy.integrate import odeint

a,b,T,R = 1.0,5.0,100.0,23.0
tau = np.linspace(0,10.0,100000.0)

def Delta(r): return a**2+r**2-r*b

def pp(r): return (1-b/r)*R-a*b*T/r
def tt(r): return -(r**2+a**2+a*a*b/r)*T-a*b*R/r

def rDotSq(r): return -(Delta(r)+T*tt(r)+R*pp(r))

def rDdot(r): return -(1/(2*r))*(rDotSq(r)+(1/r)*(2*r*(1-T**2)-b+(b*(a*T+R)**2)/(r**2)))

def computeDerivatives(y,tau):
    r = y[2]
    D = Delta(r)
    return np.array([tt(r)/D,pp(r)/D,y[3],rDdot(r)])

r0 = 80.0
D0 = Delta(r0)
y0 = np.array([tt(r0)/D0,pp(r0)/D0,r0,0.0])

sol = odeint(computeDerivatives, y0, tau)

s = np.sqrt(sol[:,2]**2+a**2)
x = s*np.sin(sol[:,1])
y = s*np.cos(sol[:,1])

data = [
    go.Scatter(
        x=x,
        y=y,
        mode="lines",
        line = dict(
            width = 2,
            color = 'rgb(0, 0, 0)'
        )
    )
]

layout = go.Layout(
    title='T={},R={},a={},b={},r={},u={}'.format(T,R,a,b,r0,y0[-1]),
    showlegend = False,
    yaxis=dict(scaleanchor="x", scaleratio=1)
)

fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename = 'plot')#, image="svg")
