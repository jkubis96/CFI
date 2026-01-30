import pandas as pd
from collections import Counter

import sys
import os
import warnings


_old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from gedspy import enrichment_heatmap

sys.stdout.close()
sys.stdout = _old_stdout

import numpy as np
import networkx as nx

def gene_interaction_network(idata:pd.DataFrame, min_con: int = 2):
        
        """
        Creates a gene or protein interaction network graph.
    
        The network is built from gene/protein interaction data. Nodes represent genes or proteins, 
        edges represent interactions, and edge colors indicate the type of interaction.
    
        Args:
            idata (pd.DataFrame): A DataFrame containing the interaction data with columns:
                - "A" (str): first gene/protein in the interaction
                - "B" (str): second gene/protein in the interaction
                - "connection_type" (str): interaction type, e.g., "gene -> protein"
            min_con (int, optional): Minimum number of connections (node degree) required 
                                     for a gene/protein to be included in the network. Default is 2.
    
        Returns:
            nx.Graph: A NetworkX graph representing the interaction network.
                - Nodes have attributes:
                    - "size": node size based on connection count (log-scaled)
                    - "color": node color (default is 'khaki')
                - Edges have attributes:
                    - "color": edge color based on interaction type
    
        Example:
            >>> G = gene_interaction_network(interactions_df, min_con=3)
            >>> nx.draw(G, with_labels=True, node_size=[G.nodes[n]['size'] for n in G.nodes()])
        """

        inter = idata
        inter = inter[["A", "B", "connection_type"]]

        dict_meta = pd.DataFrame(
            {
                "interactions": [
                    ["gene -> gene"],
                    ["protein -> protein"],
                    ["gene -> protein"],
                    ["protein -> gene"],
                    ["gene -> gene", "protein -> protein"],
                    ["gene -> gene", "gene -> protein"],
                    ["gene -> gene", "protein -> gene"],
                    ["protein -> protein", "gene -> protein"],
                    ["protein -> protein", "protein -> gene"],
                    ["gene -> protein", "protein -> gene"],
                    ["gene -> gene", "protein -> protein", "gene -> protein"],
                    ["gene -> gene", "protein -> protein", "protein -> gene"],
                    ["gene -> gene", "gene -> protein", "protein -> gene"],
                    ["protein -> protein", "gene -> protein", "protein -> gene"],
                    [
                        "gene -> gene",
                        "protein -> protein",
                        "gene -> protein",
                        "protein -> gene",
                    ],
                ],
                "color": [
                    "#f67089",
                    "#f47832",
                    "#ca9213",
                    "#ad9d31",
                    "#8eb041",
                    "#4fb14f",
                    "#33b07a",
                    "#35ae99",
                    "#36acae",
                    "#38a9c5",
                    "#3aa3ec",
                    "#957cf4",
                    "#cd79f4",
                    "#f35fb5",
                    "#f669b7",
                ],
            }
        )

        genes_list = list(inter["A"]) + list(inter["B"])

        genes_list = Counter(genes_list)

        genes_list = pd.DataFrame(genes_list.items(), columns=["features", "n"])

        genes_list = genes_list.sort_values("n", ascending=False)

        genes_list = genes_list[genes_list["n"] >= min_con]

        inter = inter[inter["A"].isin(list(genes_list["features"]))]
        inter = inter[inter["B"].isin(list(genes_list["features"]))]

        inter = inter.groupby(["A", "B"]).agg({"connection_type": list}).reset_index()

        inter["color"] = "black"

        for inx in inter.index:
            for inx2 in dict_meta.index:
                if set(inter["connection_type"][inx]) == set(
                    dict_meta["interactions"][inx2]
                ):
                    inter["color"][inx] = dict_meta["color"][inx2]
                    break

        G = nx.Graph()

        for _, row in genes_list.iterrows():
            node = row["features"]
            color = "khaki"
            weight = np.log2(row["n"] * 500)
            G.add_node(node, size=weight, color=color)

        for _, row in inter.iterrows():
            source = row["A"]
            target = row["B"]
            color = row["color"]
            G.add_edge(source, target, color=color)

        return G






def encrichment_cell_heatmap(data:pd.DataFrame,
                             fig_size = (35,25), 
                             sets = None,
                             top_n = 2,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering: str | None = 'ward',
                             scale:bool = True):
    
    """
    Creates a functional enrichment heatmap for cell types.

    This function visualizes the most significant functional terms (GO, KEGG, REACTOME, specificity)
    across different cell types.

    Args:
        data (pd.DataFrame): Input data containing columns dependent on the source (GO-TERM, KEGG, REACTOME, specificity).
        fig_size (tuple, optional): Figure size (width, height). Default is (35, 25).
        sets (list, optional): List of specific cell sets to include. Default is None (all sets).
        top_n (int, optional): Number of top terms to include per cell. Default is 2.
        test (str, optional): Name of the statistical test column. Default is 'FISH'.
        adj (str, optional): P-value adjustment method. Default is 'BH'.
        parent_inc (bool, optional): Whether to include parent terms in labels. Default is False.
        font_size (int, optional): Font size for the heatmap. Default is 16.
        clustering (str | None, optional): Clustering method for rows/columns ('ward', 'single', None). Default is 'ward'.
        scale (bool, optional): Whether to scale values before plotting. Default is True.

    Returns:
        matplotlib.figure.Figure: A heatmap figure of functional enrichment per cell type.

    Raises:
        ValueError: If the 'source' column in the data is not one of ['GO-TERM', 'KEGG', 'REACTOME', 'specificity'].

    Example:
        >>> fig = encrichment_cell_heatmap(data_df, top_n=3, parent_inc=True)
        >>> fig.savefig('cell_heatmap.svg', bbox_inches='tight')
    """
    
    if not any(x in data['source'].iloc[0] for x in ('GO-TERM', 'KEGG', 'REACTOME', 'specificity')):
        raise ValueError("Invalid value for 'source' in data. Expected: 'GO-TERM', 'KEGG', 'REACTOME' or 'specificity'.")

    set_col = 'cell'
    if  data['source'].iloc[0] == 'GO-TERM':
        term_col = 'child_name'
        parent_col = 'parent_name'

    elif data['source'].iloc[0] == 'KEGG':
        term_col = '3rd'
        parent_col = '2nd'
        
    elif data['source'].iloc[0] == 'REACTOME':
        term_col = 'pathway'
        parent_col = 'top_level_pathway'
        
    elif data['source'].iloc[0] == 'specificity':
        term_col = 'specificity'
        parent_col = 'None'

        
    title= f"Cells - {data['source'].iloc[0]}"
    
    if isinstance(sets, list):
       data[data['cell'].isin(sets)]
    
    stat_col = [x for x in data.columns if test in x and adj in x and parent_col.upper() not in x.upper()][0]
    
    if parent_inc and data['source'].iloc[0] != 'specificity':
        data[term_col] = data.apply(lambda row: f"{row[parent_col]} -> {row[term_col]}", axis=1)
        
    data = data.loc[
        data.groupby([set_col, term_col])[stat_col].idxmin()
    ].reset_index(drop = True)
    
    data = (
        data
        .sort_values(stat_col, ascending=True)
        .groupby(set_col)
        .head(top_n)
    ).reset_index(drop = True)
    

        
    if sets is None and len(list(set(data['cell']))) < 2:
        if clustering is not None:
            clustering = None
            print(
                'Clustering could not be conducted, because only one group is available in this analysis data.'
            )
        
    figure = enrichment_heatmap(data = data, 
                           stat_col = stat_col, 
                           term_col = term_col,
                           set_col = set_col,
                           sets = sets,
                           title = title,
                           fig_size = fig_size,
                           font_size = font_size,
                           scale = scale,
                           clustering = clustering) 

    return figure
      


    
    
def draw_cell_conections(data, top_n = 15):
    
    """
    Creates a cell-cell interaction network graph based on co-occurrence frequency.

    The function generates a NetworkX graph where nodes represent cell types,
    and edges represent the frequency of interactions between cells.

    Args:
        data (pd.DataFrame): A DataFrame containing columns:
            - "cell1" (str): source cell type
            - "cell2" (str): target cell type
        top_n (int, optional): Number of top interactions to include per source cell. Default is 15.

    Returns:
        nx.Graph: A NetworkX graph with attributes:
            - Nodes:
                - "size": node size (default 10)
                - "color": node color (default "#FFA07A")
            - Edges:
                - "weight": edge weight (log-transformed from frequency)
                - "color": edge color (default '#DCDCDC')
                - "alpha": edge transparency (default 0.05)

    Example:
        >>> G = draw_cell_conections(cell_interactions_df, top_n=10)
        >>> nx.draw(G, with_labels=True, node_size=[G.nodes[n]['size'] for n in G.nodes()])
    """
    
    cell_cell_df = (
        data
        .groupby(["cell1", "cell2"])
        .size()
        .reset_index(name="weight")
        .sort_values("weight", ascending=False)
    )
        
    cell_cell_df['weight'] = np.log(cell_cell_df['weight'])
        
    cell_list = list(set(list(cell_cell_df['cell1']) + list(cell_cell_df['cell2'])))
    
    df_top = (
        cell_cell_df.sort_values('weight', ascending=False)
          .groupby('cell1')
          .head(top_n)
    )
    
    
    
    G = nx.Graph()
    
    for c in cell_list:
        node = c
        color = "#FFA07A"
        weight = 10
        G.add_node(node, size=weight, color=color)
    
    for _, row in df_top.iterrows():
        source = row["cell1"]
        target = row["cell2"]
        color = '#DCDCDC'
        weight = row["weight"]/10 
    
        G.add_edge(
            source,
            target,
            weight=weight,
            color=color,
            alpha=0.05
        )
        
    nx.spring_layout(
        G,
        weight="weight",
        k=0.1,         
        iterations=500
    )
    
    return G
            
        

from matplotlib.lines import Line2D
from adjustText import adjust_text
import matplotlib.pyplot as plt


def volcano_plot_conections(
    deg_data: pd.DataFrame,
    p_adj: bool = True,
    top: int = 25,
    top_rank: str = "p_value",
    p_val: float | int = 0.05,
    lfc: float | int = 0.25,
    rescale_adj: bool = True,
    image_width: int = 12,
    image_high: int = 12,
):
    """
    Generate a volcano plot from differential expression results.

    A volcano plot visualizes the relationship between statistical significance
    (p-values or standarized p-value) and log(fold change) for each gene, highlighting
    genes that pass significance thresholds.

    Parameters
    ----------
    deg_data : pandas.DataFrame
        DataFrame containing differential expression results from calc_DEG() function.

    p_adj : bool, default=True
        If True, use adjusted p-values. If False, use raw p-values.

    top : int, default=25
        Number of top significant genes to highlight on the plot.

    top_rank : str, default='p_value'
        Statistic used primarily to determine the top significant genes to highlight on the plot. ['p_value' or 'FC']

    p_val : float | int, default=0.05
        Significance threshold for p-values (or adjusted p-values).

    lfc : float | int, default=0.25
        Threshold for absolute log fold change.

    rescale_adj : bool, default=True
        If True, rescale p-values to avoid long breaks caused by outlier values.

    image_width : int, default=12
        Width of the generated plot in inches.

    image_high : int, default=12
        Height of the generated plot in inches.

    Returns
    -------
    matplotlib.figure.Figure
        The generated volcano plot figure.

    """

    if top_rank.upper() not in ["FC", "P_VALUE"]:
        raise ValueError("top_rank must be either 'FC' or 'p_value'")

    if p_adj:
        pv = "adj_pval"
    else:
        pv = "p_val"

    deg_df = deg_data.copy()

    shift = 0.25

    p_val_scale = "-log(p_val)"

    min_minus = min(deg_df[pv][(deg_df[pv] != 0) & (deg_df["log(FC)"] < 0)])
    min_plus = min(deg_df[pv][(deg_df[pv] != 0) & (deg_df["log(FC)"] > 0)])

    zero_p_plus = deg_df[(deg_df[pv] == 0) & (deg_df["log(FC)"] > 0)]
    zero_p_plus = zero_p_plus.sort_values(by="log(FC)", ascending=False).reset_index(
        drop=True
    )
    zero_p_plus[pv] = [
        (shift * x) * min_plus for x in range(1, len(zero_p_plus.index) + 1)
    ]

    zero_p_minus = deg_df[(deg_df[pv] == 0) & (deg_df["log(FC)"] < 0)]
    zero_p_minus = zero_p_minus.sort_values(by="log(FC)", ascending=True).reset_index(
        drop=True
    )
    zero_p_minus[pv] = [
        (shift * x) * min_minus for x in range(1, len(zero_p_minus.index) + 1)
    ]

    tmp_p = deg_df[
        ((deg_df[pv] != 0) & (deg_df["log(FC)"] < 0))
        | ((deg_df[pv] != 0) & (deg_df["log(FC)"] > 0))
    ]

    del deg_df

    deg_df = pd.concat([zero_p_plus, tmp_p, zero_p_minus], ignore_index=True)

    deg_df[p_val_scale] = -np.log10(deg_df[pv])

    deg_df["top100"] = None

    if rescale_adj:

        deg_df = deg_df.sort_values(by=p_val_scale, ascending=False)

        deg_df = deg_df.reset_index(drop=True)

        eps = 1e-300
        doubled = []
        ratio = []
        for n, i in enumerate(deg_df.index):
            for j in range(1, 6):
                if (
                    n + j < len(deg_df.index)
                    and (deg_df[p_val_scale][n] + eps)
                    / (deg_df[p_val_scale][n + j] + eps)
                    >= 2
                ):
                    doubled.append(n)
                    ratio.append(
                        (deg_df[p_val_scale][n + j] + eps)
                        / (deg_df[p_val_scale][n] + eps)
                    )

        df = pd.DataFrame({"doubled": doubled, "ratio": ratio})
        df = df[df["doubled"] < 100]

        df["ratio"] = (1 - df["ratio"]) / 5
        df = df.reset_index(drop=True)

        df = df.sort_values("doubled")

        if len(df["doubled"]) == 1 and 0 in df["doubled"]:
            df = df
        else:
            doubled2 = []

            for l in df["doubled"]:
                if l + 1 != len(doubled) and l + 1 - l == 1:
                    doubled2.append(l)
                    doubled2.append(l + 1)
                else:
                    break

            doubled2 = sorted(set(doubled2), reverse=True)

        if len(doubled2) > 1:
            df = df[df["doubled"].isin(doubled2)]
            df = df.sort_values("doubled", ascending=False)
            df = df.reset_index(drop=True)
            for c in df.index:
                deg_df.loc[df["doubled"][c], p_val_scale] = deg_df.loc[
                    df["doubled"][c] + 1, p_val_scale
                ] * (1 + df["ratio"][c])

    deg_df.loc[(deg_df["log(FC)"] <= 0) & (deg_df[pv] <= p_val), "top100"] = "bisque"
    deg_df.loc[(deg_df["log(FC)"] > 0) & (deg_df[pv] <= p_val), "top100"] = "skyblue"
    deg_df.loc[deg_df[pv] > p_val, "top100"] = "lightgray"

    if lfc > 0:
        deg_df.loc[
            (deg_df["log(FC)"] <= lfc) & (deg_df["log(FC)"] >= -lfc), "top100"
        ] = "lightgray"

    down_int = len(
        deg_df["top100"][(deg_df["log(FC)"] <= lfc * -1) & (deg_df[pv] <= p_val)]
    )
    up_int = len(deg_df["top100"][(deg_df["log(FC)"] > lfc) & (deg_df[pv] <= p_val)])

    deg_df_up = deg_df[deg_df["log(FC)"] > 0]

    if top_rank.upper() == "P_VALUE":
        deg_df_up = deg_df_up.sort_values([pv, "log(FC)"], ascending=[True, False])
    elif top_rank.upper() == "FC":
        deg_df_up = deg_df_up.sort_values(["log(FC)", pv], ascending=[False, True])

    deg_df_up = deg_df_up.reset_index(drop=True)

    n = -1
    l = 0
    while True:
        n += 1
        if deg_df_up["log(FC)"][n] > lfc and deg_df_up[pv][n] <= p_val:
            deg_df_up.loc[n, "top100"] = "dodgerblue"
            l += 1
        if l == top or deg_df_up[pv][n] > p_val:
            break

    deg_df_down = deg_df[deg_df["log(FC)"] <= 0]

    if top_rank.upper() == "P_VALUE":
        deg_df_down = deg_df_down.sort_values([pv, "log(FC)"], ascending=[True, True])
    elif top_rank.upper() == "FC":
        deg_df_down = deg_df_down.sort_values(["log(FC)", pv], ascending=[True, True])

    deg_df_down = deg_df_down.reset_index(drop=True)

    n = -1
    l = 0
    while True:
        n += 1
        if deg_df_down["log(FC)"][n] < lfc * -1 and deg_df_down[pv][n] <= p_val:
            deg_df_down.loc[n, "top100"] = "tomato"

            l += 1
        if l == top or deg_df_down[pv][n] > p_val:
            break

    deg_df = pd.concat([deg_df_up, deg_df_down])
    
        


    que = ["lightgray", "bisque", "skyblue", "tomato", "dodgerblue"]

    deg_df = deg_df.sort_values(
        by="top100", key=lambda x: x.map({v: i for i, v in enumerate(que)})
    )

    

    
    deg_df = deg_df.reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(image_width, image_high))

    plt.scatter(
        x=deg_df["log(FC)"], y=deg_df[p_val_scale], color=deg_df["top100"], zorder=2
    )
    
    
    tl = deg_df[p_val_scale][deg_df[pv] >= p_val]
    
    if len(tl) > 0:
        
        line_p = np.max(tl)
    
    else:
        line_p = np.min(deg_df[p_val_scale])
    
    plt.plot(
        [max(deg_df["log(FC)"]) * -1.1, max(deg_df["log(FC)"]) * 1.1],
        [line_p, line_p],
        linestyle="--",
        linewidth=3,
        color="lightgray",
        zorder=1,
    )

    if lfc > 0:
        plt.plot(
            [lfc * -1, lfc * -1],
            [-3, max(deg_df[p_val_scale]) * 1.1],
            linestyle="--",
            linewidth=3,
            color="lightgray",
            zorder=1,
        )
        plt.plot(
            [lfc, lfc],
            [-3, max(deg_df[p_val_scale]) * 1.1],
            linestyle="--",
            linewidth=3,
            color="lightgray",
            zorder=1,
        )

    plt.xlabel("log(FC)")
    plt.ylabel(p_val_scale)
    plt.title("Volcano plot")

    plt.ylim(min(deg_df[p_val_scale]) - 5, max(deg_df[p_val_scale]) * 1.3)

    texts = [
        ax.text(deg_df["log(FC)"][i], deg_df[p_val_scale][i], deg_df["feature"][i])
        for i in deg_df.index
        if deg_df["top100"][i] in ["dodgerblue", "tomato"]
    ]

    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", alpha=0.5))

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="top-upregulated",
            markerfacecolor="dodgerblue",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="top-downregulated",
            markerfacecolor="tomato",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="upregulated",
            markerfacecolor="skyblue",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="downregulated",
            markerfacecolor="bisque",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="non-significant",
            markerfacecolor="lightgray",
            markersize=10,
        ),
    ]

    ax.legend(handles=legend_elements, loc="upper right")
    ax.grid(visible=False)

    ax.annotate(
        f"\nmin {pv} = " + str(p_val),
        xy=(0.025, 0.975),
        xycoords="axes fraction",
        fontsize=12,
    )

    if lfc > 0:
        ax.annotate(
            "\nmin log(FC) = " + str(lfc),
            xy=(0.025, 0.95),
            xycoords="axes fraction",
            fontsize=12,
        )

    ax.annotate(
        "\nDownregulated: " + str(down_int),
        xy=(0.025, 0.925),
        xycoords="axes fraction",
        fontsize=12,
        color="black",
    )

    ax.annotate(
        "\nUpregulated: " + str(up_int),
        xy=(0.025, 0.9),
        xycoords="axes fraction",
        fontsize=12,
        color="black",
    )

    plt.show()

    return fig


