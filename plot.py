import plotly.offline as py
import plotly.graph_objs as go
import numpy as np
import pandas as pd

data = pd.read_csv('data.csv')
data = data.iloc[1:]
data = [
    go.Scatter(
        x=data['x'],
        y=data['y'],
        mode="lines",
        line = dict(
            width = 2,
            color = 'rgb(0, 0, 0)'
        )
    )
]

layout = go.Layout(
    title='T=0.0,R=-1.0,a=2.0,b=0.1,t=0.0,phi=0.0,r=1.0,u=0.0',
    showlegend = False,
    yaxis=dict(scaleanchor="x", scaleratio=1)
)

fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename = 'plot.html')