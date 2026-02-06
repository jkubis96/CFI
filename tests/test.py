import os

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

import networkx as nx
from jdti import COMPsc
from matplotlib.figure import Figure

from cfi_toolkit import (
    CellFunCon,
    compare_connections,
    draw_cell_conections,
    encrichment_cell_heatmap,
    gene_interaction_network,
)


@pytest.fixture(scope="session")
def data_dir():
    return os.path.join(os.getcwd(), "data")


@pytest.fixture(scope="session")
def jseq_s1(data_dir):
    jseq = COMPsc.project_dir(data_dir, ["s1"])
    jseq.load_sparse_from_projects(normalized_data=True)
    return jseq


@pytest.fixture(scope="session")
def jseq_s2(data_dir):
    jseq = COMPsc.project_dir(data_dir, ["s2"])
    jseq.load_sparse_from_projects(normalized_data=True)
    return jseq


@pytest.fixture(scope="session")
def instance1(jseq_s1, tmp_path_factory):
    tmp = tmp_path_factory.mktemp("cellfuncon_projects")
    os.chdir(tmp)

    inst = CellFunCon(jseq_s1)
    inst.calculate_cell_connections()
    inst.save_project("test")

    inst = CellFunCon.load_project("test.psc")

    inst.calculate_cells_markers(min_exp=0, min_pct=0.05, n_proc=2)

    inst.enrich_cells_fucntionality(p_value=0.05, log_fc=0.25, top_max=200)

    return inst


@pytest.fixture(scope="session")
def instance2(jseq_s2):
    inst = CellFunCon(jseq_s2)
    inst.calculate_cell_connections()
    return inst


# =========================
# TESTY PODSTAWOWE
# =========================


def test_cell_connections(instance1):
    df = instance1.get_cell_connections()
    assert isinstance(df, pd.DataFrame)
    assert len(df.index) > 1


def test_included_cells(instance1):
    cells = instance1.get_included_cells()
    assert isinstance(cells, list)
    assert len(cells) > 0


# =========================
# TESTY ENRICHMENT
# =========================


@pytest.mark.parametrize("data_type", ["GO-TERM", "KEGG", "REACTOME", "specificity"])
def test_enrichment_data(instance1, data_type):
    enr = instance1.get_enrichment_data(
        data_type=data_type,
        p_value=1,
        test="FISH",
        adj="BH",
        parent_inc=False,
        top_n=20,
    )

    assert isinstance(enr, pd.DataFrame)
    assert len(enr.index) > 1


def test_heatmap_plot(instance1):
    enr = instance1.get_enrichment_data(
        data_type="GO-TERM", test="FISH", adj="BH", parent_inc=False, top_n=10
    )

    fig = encrichment_cell_heatmap(data=enr, fig_size=(3, 3), top_n=3, scale=True)

    assert isinstance(fig, Figure)


def test_gene_interaction_network(instance1):
    cell = instance1.get_included_cells()[0]
    gi = instance1.get_gene_interactions(cell)

    fig = gene_interaction_network(idata=gi, min_con=1)

    assert isinstance(fig, nx.Graph)


def test_draw_cell_connections(instance1):
    cell_con = instance1.get_cell_connections()
    fig = draw_cell_conections(cell_con)

    assert isinstance(fig, nx.Graph)


def test_compare_connections(instance1, instance2):
    comparison = compare_connections(
        instances_dict={"healthy": instance2, "disease": instance1},
        connection_type=[
            "Adhesion-Adhesion",
            "Gap-Gap",
            "Ligand-Ligand",
            "Ligand-Receptor",
            "Receptor-Receptor",
            "Undefined",
        ],
    )

    assert isinstance(comparison, pd.DataFrame)
    assert len(comparison.index) > 1
