from cfi import CellFunCon

from jdti import COMPsc

import os

jseq_object = COMPsc.project_dir(os.path.join(os.getcwd(), 'data'), ['s1'])
jseq_object.load_sparse_from_projects(normalized_data=True)

instance1 = CellFunCon(jseq_object)

instance1.calculate_cell_connections()


instance1.save_project('test')

instance1 = CellFunCon.load_project('test.psc')

instance1.calculate_cells_markers(min_exp = 0, 
                        min_pct = 0.05, 
                        n_proc=10)



instance1.enrich_cells_fucntionality(p_value = 0.05, 
                                     log_fc = 0.25, 
                                     top_max = 500)



enr = instance1.get_enrichment_data( 
                        data_type = 'GO-TERM', 
                        p_value = 0.05, 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)

enr.to_csv("go.csv", index=False)


from cfi import encrichment_cell_heatmap

fig1 = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = None,
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)

fig1.savefig(
    "heatmap.svg",
    format="svg",
    dpi=300,
    bbox_inches="tight"
)



enr = instance1.get_enrichment_data( 
                        data_type = 'KEGG', 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)


enr.to_csv("kegg.csv", index=False)


fig2_1 = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = None,
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)

fig2_1.savefig(
    "heatmap2.svg",
    format="svg",
    dpi=300,
    bbox_inches="tight"
)


fig2_2 = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = instance1.get_included_cells(),
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)


fig2_2.savefig(
    "heatmap3.svg",
    format="svg",
    dpi=300,
    bbox_inches="tight"
)




enr = instance1.get_enrichment_data( 
                        data_type = 'REACTOME', 
                        p_value = 0.05, 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)

enr.to_csv("reactome.csv", index=False)

fig3 = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = instance1.get_included_cells(),
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)

fig3.savefig(
    "heatmap4.svg",
    format="svg",
    dpi=300,
    bbox_inches="tight"
)





enr = instance1.get_enrichment_data( 
                        data_type = 'specificity', 
                        p_value = 1, 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)

enr.to_csv("spec.csv", index=False)

fig4 = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = instance1.get_included_cells(),
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)

fig4.savefig(
    "heatmap5.svg",
    format="svg",
    dpi=300,
    bbox_inches="tight"
)



instance1.get_included_cells()


cell_int = instance1.get_gene_interactions('STRIATUM_1 # s1')

cell_int.to_csv("cell_int.csv", index=False)

from cfi import gene_interaction_network



fig5 = gene_interaction_network(idata = cell_int, min_con = 2)

from JVG import JVG

nt = JVG.NxEditor(fig5)
nt.edit()


instance1.calculate_cell_connections()

cell_con = instance1.get_cell_connections()

cell_con.to_csv("cell_con.csv", index=False)


from cfi import draw_cell_conections

fig6 = draw_cell_conections(cell_con)


from JVG import JVG

nt = JVG.NxEditor(fig6)
nt.edit()



jseq_object = COMPsc.project_dir(os.path.join(os.getcwd(), 'data'), ['s2'])
jseq_object.load_sparse_from_projects(normalized_data=True)



instance2 = CellFunCon(jseq_object)




instance2.calculate_cell_connections()






from cfi import compare_connections


instances_dict = {'healthy':instance2,
                  'disease':instance1}

comparison = compare_connections(instances_dict=instances_dict, 
                                 cells_compartment = None, 
                                 connection_type  = ['Adhesion-Adhesion',
                                                     'Gap-Gap',
                                                     'Ligand-Ligand',
                                                     'Ligand-Receptor',
                                                     'Receptor-Receptor',
                                                     'Undefined'])


comparison.to_csv("comparison.csv", index=False)


from cfi import volcano_plot_conections

fig7 = volcano_plot_conections(
    deg_data = comparison,
    p_adj = True,
    top = 25,
    top_rank = "p_value",
    p_val = 0.05,
    lfc = 0.25,
    rescale_adj = True,
    image_width = 12,
    image_high = 12,
)
            

fig7.savefig(
    "volcano.svg",
    format="svg",
    dpi=300,
    bbox_inches="tight"
)
          
            
            
            
        

