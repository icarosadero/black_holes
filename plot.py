import plotly.offline as py
import plotly.graph_objs as go
import numpy as np
import pandas as pd

p = pd.read_csv('constants.csv').columns
p = list(map(float,p))
string = 'T={:.2f},R={:.2f},a={:.2f},b={:.2f},r={:.2f},u={}'.format(p[0],p[1],p[2],p[3],p[4],0.0)
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
    title=string,
    showlegend = False,
    yaxis=dict(scaleanchor="x", scaleratio=1)
)

fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename = string+'.html')#, image="svg")
