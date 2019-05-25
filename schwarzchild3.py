import numpy as np
import plotly.offline as py
from plotly import tools
import plotly.graph_objs as go
from scipy.integrate import odeint

"""
Paramemters
m2, r0 = -2.0, 1.0
r1,r2 = -1.0, 1.0
M = 2*np.sqrt(-r0*m2*(r1+r2)/2)
"""
r0,M,m2 = 1.0, 5.0, 1.0
a,b,c = -r0*m2/2.0, (M*M)/4.0, -(3.0/8)*M*M*r0
#Try circular orbit
#a,b,c = -1, 0, 100 works with M=5**4

#t = np.linspace(0,100000,10000)
t = np.linspace(0,1e+5,1000)
y0 = np.array([243,0.0,0.0])

#def poly(r): return (r0*m2/2)*(r-r1)*(r-r2)
def poly(r): return a*r**2.0 + b*r + c
roots = np.roots([a,b,c])
print('Roots: ',roots)

def R(r):
    return poly(r)/r**4

def P(r):
    return M/(2*r*r) #for some reason works for M->M**4

def f(y,t):
    """
    Order: r,u,phi
    """
    r,u,phi = y
    return np.array([u,R(r),P(r)])

sol = odeint(f, y0, t)

np.savetxt('orbitValues3.txt', sol)

#PLOT-----------------------------------------------------

rex = np.linspace(min(roots)-2.0, max(roots)+2.0)
data = [
    go.Scatter(
        x = rex,
        y = poly(rex),
        mode = 'lines',
        line=dict(
            color='purple'
        ) 
),
    go.Scatterpolar(
        r = sol[:,0],
        theta = sol[:,2]*180./np.pi,
        mode = 'lines',
        line =  dict(
            color = 'violet'
        ),
        subplot='polar'
    )
]

layout = go.Layout(
    title = "r0:{:2.3f}, m2:{:2.3f}, M:{:2.3f}. Roots: {:2.3f}, {:2.3f}".format(
        y0[0], m2, M, roots[0], roots[1]),
    polar = dict(
      radialaxis = dict(
        angle = 45
      ),
      angularaxis = dict(
        direction = "clockwise",
        period = 6
      )
    )
)

py.plot(go.Figure(data = [data[0]]), filename = 'quadratic4.html')
py.plot(go.Figure(data = [data[1]], layout=layout), filename = 'orbit4.html')