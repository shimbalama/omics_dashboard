import pandas as pd
import plotly.express as px
import dash
import dash_core_components as dcc
import dash_html_components as html

# Load data
df = pd.read_csv('data.csv')

# Calculate log2 fold change and average expression
df['log2FC'] = df['sample1'] - df['sample2']
df['avgExpr'] = (df['sample1'] + df['sample2']) / 2

# Calculate p-values
df['pvalue'] = [0.01, 0.1, 0.001, 0.05, 0.2] # example p-values

# Define significance threshold
alpha = 0.05

# Color points based on significance
df['color'] = ['red' if pval <= alpha else 'blue' for pval in df['pvalue']]

# Create MA plot
fig = px.scatter(df, x='avgExpr', y='log2FC', color='color',
                 labels={'avgExpr': 'Average Expression', 'log2FC': 'Log2 Fold Change'},
                 title='MA Plot')

# Set layout
fig.update_layout(
    xaxis_title='Average Expression',
    yaxis_title='Log2 Fold Change',
    plot_bgcolor='white'
)

# Create app
app = dash.Dash(__name__)

# Create layout
app.layout = html.Div([
    dcc.Graph(id='ma-plot', figure=fig)
])

# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
