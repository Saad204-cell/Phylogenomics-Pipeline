import pandas as pd
import plotly.graph_objects as go
from Bio import Phylo
import io

df = pd.read_csv("output/metadata.csv")
tree = Phylo.read("output/tree.nwk", "newick")

# Render ASCII tree for embedding
buf = io.StringIO()
Phylo.draw_ascii(tree, file=buf)
tree_ascii = buf.getvalue()

fig = go.Figure(data=[go.Bar(
    x=df["Accession"],
    y=df["Length"],
    text=df["Organism"],
    hovertext=df["Organism"]
)])
fig.update_layout(
    title="Sequence Length by Accession (hover for organism)",
    template="plotly_dark"
)

html_parts = [
    "<html><body style='background:#111;color:#eee;font-family:monospace;padding:20px;'>",
    "<h1>Phylogenomic Report</h1>",
    "<h2>Phylogenetic Tree (ASCII View)</h2>",
    f"<pre style='background:#222;padding:15px;border-radius:5px;'>{tree_ascii}</pre>",
    "<h2>Sequence Metadata</h2>",
    fig.to_html(full_html=False, include_plotlyjs="cdn"),
    "</body></html>"
]

with open("output/phylogenetic_report.html", "w") as f:
    f.write("\n".join(html_parts))

print("Report saved to output/phylogenetic_report.html")
