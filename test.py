import plotly.express as px
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd

# Load example data
df = px.data.tips()

# Create box plot with title
fig = px.box(df, x="day", y="total_bill", title="Box Plot of Total Bill by Day")
fig.update_layout(
    margin=dict(t=200) # adjust top margin to increase distance between title and plot
)
# Create Dash app
app = dash.Dash(__name__)

# Define app layout
app.layout = html.Div(children=[
    html.H1(children='Box Plot Example'),
    dcc.Graph(
        id='example-graph',
        figure=fig
    )
])

# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
