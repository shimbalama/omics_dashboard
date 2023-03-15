import plotly.express as px
df = px.data.tips()
print(df)
fig = px.box(df, x="time", y="total_bill")
fig.show()