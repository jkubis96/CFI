import pandas as pd
from tqdm import tqdm
import sys
import os

_old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from gedspy import Enrichment
from gedspy import Analysis

sys.stdout.close()
sys.stdout = _old_stdout

import pickle
import copy


class CellFunCon:
    
    """
    A class to perform cell-type functional analysis and enrichment based on a JDTI/COMPsc dataset.

    This class provides methods to calculate marker genes for cell types, perform functional enrichment 
    (GO, KEGG, REACTOME, STRING, IntAct), and compute cell-cell interaction networks. 
    Projects can also be saved and loaded via pickle.

    Attributes
    ----------
    jdti : object
        JDTI or COMPsc object containing normalized single-cell data.
    cells_markers : pd.DataFrame or None
        DataFrame containing marker genes per cell type after calculation.
    enr_full_info : Enrichment
        Enrichment object containing all genes available for enrichment analysis.
    cells_enrichment : dict
        Dictionary storing enrichment results per cell type.
    cells_connection : pd.DataFrame or None
        DataFrame storing calculated cell-cell interaction information.
    mt_genes : bool
        Whether mitochondrial genes are included (default False).
    ribo_genes : bool
        Whether ribosomal genes are included (default False).
    """
    
    def __init__(self, jdti_object):
        
        """
        Initializes the CellFunCon object with a COMPsc/JDTI object.

        Parameters
        ----------
        jdti_object : object
            A COMPsc or JDTI object with normalized single-cell data.
        """

        self.jdti = jdti_object
        self.cells_markers = None
        self.enr_full_genes = None
        self.cells_connection = None
        self.mt_genes = False
        self.ribo_genes = False
        
                
        names =  self.jdti.normalized_data.loc[self.jdti.normalized_data.select_dtypes(include='number').sum(axis=1) > 0].index.tolist()
        names = list(set(names))
        
        if self.mt_genes is False:
            names = [x for x in names if 'MT-' not in x.upper()]
        if self.ribo_genes is False:
            names = [x for x in names if 'RPS' != x[:3].upper()]
            names = [x for x in names if 'RPL' != x[:3].upper()]
    
        enr=Enrichment()
        enr.select_features(names)
        
        self.enr_full_info = enr
        


    def save_project(self, filename):
        
        """
        Saves the current CellFunCon project as a pickle file.

        Parameters
        ----------
        filename : str
            Path to save the project (e.g., 'project_name.psc').

        Example
        -------
        >>> self.save_project('my_project.psc')
        """
        
        with open(f'{filename}.psc', "wb") as f:
            pickle.dump(self, f)
        print(f"Project saved as {filename}")


    @classmethod
    def load_project(cls, filename):
        
        """
        Loads a previously saved CellFunCon project from a pickle file.

        Parameters
        ----------
        filename : str
            Path to the saved pickle file.

        Returns
        -------
        CellFunCon
            Loaded CellFunCon self.

        Raises
        ------
        TypeError
            If the loaded object is not a CellFunCon self.
        ValueError
            If the file is not a valid CellFunCon project file.

        Example
        -------
        >>> self = CellFunCon.load_project('my_project.psc')
        """


        if '.psc' in filename:
            with open(filename, "rb") as f:
                obj = pickle.load(f)
            if not isinstance(obj, cls):
                raise TypeError("Plik nie zawiera obiektu Project")
            print(f"Projekt wczytany z {filename}")
            return obj
        else:
            raise ValueError('Project not belong to CellFunCon project data.')



    def calculate_cells_markers(self, 
                         min_exp = 0, 
                         min_pct = 0.05, 
                         n_proc=10):
        
        """
        Calculates marker genes for each cell type based on expression thresholds.

        Parameters
        ----------
        min_exp : float, optional
            Minimum expression level to consider a gene (default 0).
        min_pct : float, optional
            Minimum fraction of cells expressing a gene (default 0.05).
        n_proc : int, optional
            Number of parallel processes to use (default 10).

        Notes
        -----
        The results are stored in the `cells_markers` attribute.
        """
            
        self.jdti.calculate_difference_markers(min_exp=min_exp, 
                                             min_pct=min_pct, 
                                             n_proc=n_proc, 
                                             force=True)
    
        self.cells_markers = self.jdti.var_data
    

    def enrich_cells_fucntionality(self, 
                                   p_value = 0.05, 
                                   log_fc = 0.25, 
                                   top_max = 500):
        
        """
        Performs functional enrichment analysis for each cell type based on marker genes.

        Parameters
        ----------
        p_value : float, optional
            Maximum adjusted p-value for significant genes (default 0.05).
        log_fc : float, optional
            Minimum log fold-change threshold for marker genes (default 0.25).
        top_max : int, optional
            Maximum number of top marker genes per cell type to consider (default 500).

        Raises
        ------
        ValueError
            If `cells_markers` is not defined.

        Notes
        -----
        This method populates `cells_enrichment` with results for GO-TERM, KEGG, REACTOME, 
        STRING, IntAct, and specificity analyses.
        """
        
        
        if isinstance(self.cells_markers, pd.DataFrame):
            
            markers = self.cells_markers
            cells=set(markers['valid_group'])
            
            data_dict = {}
            
            max_c = len(cells)
            for n, c in enumerate(cells):
                print(f'\nAnalysis {n+1} of {max_c} cells --> {c} \n')
                tmp = markers[(markers['valid_group'] == c) & (markers['adj_pval'] <= p_value) & (markers['log(FC)'] > log_fc)]
                names = list(set(tmp['feature']))
                
                tmp = tmp[tmp['feature'].isin(names)]
    
                if len(tmp.index) < 10:
                    tmp = markers[(markers['valid_group'] == c) & (markers['p_val'] <= p_value) & (markers['log(FC)'] > log_fc)]
                    names = list(set(tmp['feature']))
                    
                    tmp = tmp[tmp['feature'].isin(names)]
                
                tmp = tmp.sort_values('esm', ascending=False).head(top_max)
    
               
                data_dict[c] = {}
                enr=copy.copy(self.enr_full_info)
                enr.genome = enr.genome[enr.genome['found_names'].isin(list(set(tmp['feature'])))].reset_index(drop = True)
                enr.enriche_specificiti()
                enr.enriche_KEGG()
                enr.enriche_GOTERM()
                enr.enriche_REACTOME()
                enr.enriche_IntAct()
                enr.enriche_STRING()
                enr.enriche_specificiti()

                data = enr.get_results()
                del enr
                
                ans = Analysis(data)
                ans.gene_interaction()
                ans.features_specificity()
                ans.REACTOME_overrepresentation()
                ans.KEGG_overrepresentation()
                ans.GO_overrepresentation()
                ans.features_specificity()
                
                data_dict[c] = ans.get_full_results()
            
            self.cells_enrichment = data_dict
            
        else:
            raise ValueError('`self.cells_markers` not defined. Use `self.cells_markers` to provide markers.')
            
            
 

    def get_enrichment_data(self, 
                            data_type = 'GO-TERM', 
                            p_value = 0.05, 
                            test = 'FISH', 
                            adj = 'BH', 
                            parent_inc = False, 
                            top_n = 50):
        
        """
        Retrieves enrichment results for all cells in a unified DataFrame.

        Parameters
        ----------
        data_type : str, optional
            Type of enrichment to retrieve ('GO-TERM', 'KEGG', 'REACTOME', 'specificity').
        p_value : float, optional
            Maximum p-value threshold (default 0.05).
        test : str, optional
            Name of the statistical test column to use (default 'FISH').
        adj : str, optional
            P-value adjustment method (default 'BH').
        parent_inc : bool, optional
            Whether to include parent terms in the results (default False).
        top_n : int, optional
            Maximum number of terms per cell type to include (default 50).

        Returns
        -------
        pd.DataFrame
            DataFrame containing filtered enrichment results with a 'cell' column indicating cell type.

        Raises
        ------
        ValueError
            If `data_type` is not one of the expected values.
        """
        
        if not any(x in data_type for x in ('GO-TERM', 'KEGG', 'REACTOME', 'specificity')):
            raise ValueError("Invalid value for 'data_type'. Expected: 'GO-TERM', 'KEGG', 'REACTOME' or 'specificity'.")
           
        if  data_type == 'GO-TERM':
            parent_col = 'parent'
    
        elif data_type == 'KEGG':
            parent_col = '2nd'
            
        elif data_type == 'REACTOME':
            parent_col = 'top_level'
            
        elif data_type == 'specificity':
            parent_col = 'None'
    
    
        pdl = []
        for i in self.cells_enrichment.keys():
            print(i) 
            if data_type == 'specificity':
                tmp_dict = self.cells_enrichment[i]['statistics'][data_type]
                tmp = []
                for k in tmp_dict.keys():
                    if k != 'HPA_subcellular_location':
                        tmp.append(pd.DataFrame(tmp_dict[k]))
                
                tmp = pd.concat(tmp)
                    
            else:
                tmp = pd.DataFrame(self.cells_enrichment[i]['statistics'][data_type])
    
                                 
            cols = [x for x in tmp.columns if test in x and adj in x]
            cols = sorted(cols, reverse=True)
            if parent_inc is False:
                cols = [x for x in cols if parent_col not in x.lower()]
            
            mask = (tmp[cols] <= p_value).all(axis=1)
            tmp = tmp.loc[mask]
            tmp['cell'] = i
            tmp = tmp.sort_values(by=['cell'] + cols, ascending=True)
    
            pdl.append(tmp.head(top_n))
            
            
        df = pd.concat(pdl)
        df['source'] = data_type
        df = df.reset_index(drop = True)
        
        return df
    

    def get_included_cells(self):
        
        """
        Returns the list of cell types included in the enrichment analysis.

        Returns
        -------
        list
            List of cell type names.

        Example
        -------
        >>> self.get_included_cells()
        ['CellType1', 'CellType2', ...]
        """
        
        cl = []
        for i in self.cells_enrichment.keys():
            print(i) 
            cl.append(i)
            
        return cl

    
    def get_gene_interactions(self, cell_name):
        
        """
        Retrieves gene or protein interaction data for a specific cell type.

        Parameters
        ----------
        cell_name : str
            Name of the cell type.

        Returns
        -------
        pd.DataFrame
            DataFrame containing interactions for the specified cell.

        Example
        -------
        >>> self.get_gene_interactions('CellType1')
        """
        
        tmp = pd.DataFrame(self.cells_enrichment[cell_name]['statistics']['interactions'])
        
        return tmp


    def calculate_cell_connections(self):
        
        """
        Calculates cell-cell interaction connections based on gene/protein co-expression.

        Notes
        -----
        Populates `cells_connection` with a DataFrame containing interactions between all pairs of cells.
        Each row represents an interaction between two cells and the involved genes/proteins.

        Raises
        ------
        ValueError
            If `normalized_data` is not defined in the JDTI object.
        """
        
        if isinstance(self.jdti.normalized_data, pd.DataFrame):
            
            
            cells=set(self.jdti.normalized_data.columns)
          
            data_dict = {}
            
            for c in tqdm(cells):
                
                  
                tmp = self.jdti.normalized_data.loc[:,c]
                names =  tmp.loc[tmp.select_dtypes(include='number').sum(axis=1) > 0].index.tolist()
                names = list(set(names))
            
               
                enr=copy.copy(self.enr_full_info)
                enr.genome = enr.genome[enr.genome['found_names'].isin(names)].reset_index(drop = True)
                enr.enriche_CellCon()
                data = enr.get_results()
                del enr
                
                data_dict[c] = data['CellConnections']
                
            full_data = []
            for c1 in tqdm(cells):
                for c2 in cells:
                    if c1 != c2:
                        c1_d = pd.DataFrame(data_dict[c1]['interactor2'])
                        c2_d = pd.DataFrame(data_dict[c2]['interactor1'])
                        
                        mutual_lr = c1_d['interaction'][c1_d['interaction'].isin(list(c2_d['interaction']))]
                        
                        to_ret = c1_d[
                            c1_d['interaction'].isin(list(mutual_lr))
                        ].drop(
                            ['Species', 'protein_id_1', 'protein_id_2', 'found_names_2'], 
                            axis=1
                        ).reset_index(drop=True)
                        
                        to_ret = to_ret.rename(columns={'found_names_1': 'interactor1'})
                        c2_subset = c2_d[['interaction', 'found_names_2']].rename(columns={'found_names_2': 'interactor2'})
                        
                        to_ret = to_ret.merge(c2_subset, on='interaction', how='left')
                        to_ret['cell1'] = c1
                        to_ret['cell2'] = c2
                        
                        full_data.append(to_ret)
                        
            self.cells_connection = pd.concat(full_data)
            
                        
            
        else:
            raise ValueError('`self.cells_markers` not defined. Use `self.cells_markers` to provide markers.')
            
            
    def get_cell_connections(self):
        
        """
        Returns the calculated cell-cell interaction connections.

        Returns
        -------
        pd.DataFrame
            DataFrame containing cell-cell interactions.

        Example
        -------
        >>> connections = self.get_cell_connections()
        """
        
        return self.cells_connection 





def compare_connections(instances_dict:dict, 
                        cells_compartment:dict | None = None, 
                        connection_type: list  = ['Adhesion-Adhesion',
                                                  'Gap-Gap',
                                                  'Ligand-Ligand',
                                                  'Ligand-Receptor',
                                                  'Receptor-Receptor',
                                                  'Undefined']):
    
    """
    Compare gene expression between two instances based on their cell connections.

    This function compares normalized gene expression data from exactly two
    instances stored in ``instances_dict``. Optionally, the comparison can be
    restricted to specific cell compartments for each instance. Differential
    expression analysis is performed using ``jdti.calc_DEG``.

    Parameters
    ----------
    instances_dict : dict
        Dictionary containing exactly two objects. Each object must have:
        
        - ``jdti.normalized_data`` : pandas.DataFrame
            Gene expression matrix with genes as rows and cells as columns.
        - ``cells_connection`` : pandas.DataFrame
            DataFrame containing at least the columns ``'interactor1'`` and
            ``'interactor2'``.

        The dictionary keys are used as group labels in the comparison.

    cells_compartment : dict or None, optional
        Dictionary mapping each key in ``instances_dict`` to a list of cell names
        to be used for the comparison. If ``None``, all cells are used and genes
        are filtered based on cell–cell connections.

    Returns
    -------
    pandas.DataFrame
        Differential expression results returned by ``calc_DEG``, filtered to
        include only rows where ``valid_group`` matches the first key in
        ``instances_dict``.

    Raises
    ------
    ValueError
        If any cell specified in ``cells_compartment`` is not present in the
        corresponding ``normalized_data`` columns.

    Notes
    -----
    - Only genes common to both instances are considered.
    - When ``cells_compartment`` is ``None``, genes are further restricted to
      those appearing in the cell–cell interaction networks of either instance.
    - The function assumes exactly two entries in ``instances_dict``.
    - Differential expression is computed with ``min_exp=0`` and ``min_pct=0.1``.

    See Also
    --------
    jdti.calc_DEG : Function used to compute differential expression.
    """
    
    import pandas as pd
    from jdti import calc_DEG
    
    if isinstance(cells_compartment, dict):
        
        keys_list = list(instances_dict.keys())
        tmp1 = instances_dict[keys_list[0]].jdti.normalized_data.copy()
        cells = cells_compartment[keys_list[0]]
        if any(cell not in tmp1.columns for cell in cells):
            raise ValueError('Any of {keys_list[0]} cells in dictionary "cells_compartment" do not occur!')
        tmp1 = tmp1.loc[:,cells]
        tmp1.columns = [keys_list[0]]*len(tmp1.columns)
        
        tmp2 = instances_dict[keys_list[1]].jdti.normalized_data.copy()
        cells = cells_compartment[keys_list[1]]
        if any(cell not in tmp2.columns for cell in cells):
            raise ValueError('Any of {keys_list[1]} cells in dictionary "cells_compartment" do not occur!')
        tmp2 = tmp2.loc[:,cells]
        tmp2.columns = [keys_list[1]]*len(tmp2.columns)
            
        common_idx = tmp1.index.intersection(tmp2.index)

        tmp1 = tmp1.loc[common_idx]
        tmp2 = tmp2.loc[common_idx]
        
        concat_df = pd.concat([tmp1, tmp2],  axis=1)
            
    else:
        
        keys_list = list(instances_dict.keys())
        tmp1 = instances_dict[keys_list[0]].jdti.normalized_data.copy()
        tmp1.columns = [keys_list[0]]*len(tmp1.columns)
        
        tmp2 = instances_dict[keys_list[1]].jdti.normalized_data.copy()
        tmp2.columns = [keys_list[1]]*len(tmp2.columns)
            
        common_idx = tmp1.index.intersection(tmp2.index)

        tmp1 = tmp1.loc[common_idx]
        tmp2 = tmp2.loc[common_idx]
        
        concat_df = pd.concat([tmp1, tmp2],  axis=1)
        
        
    tmp_df_1 = instances_dict[keys_list[0]].cells_connection
    tmp_df_2 = instances_dict[keys_list[1]].cells_connection

    
    tmp_df_1['directionality'] = [x if x is not None else 'Undefined' for x in tmp_df_1['directionality']]
    tmp_df_2['directionality'] = [x if x is not None else 'Undefined' for x in tmp_df_2['directionality']]

        
    tmp_df_1 = tmp_df_1[tmp_df_1['directionality'].isin(connection_type)]
    tmp_df_2 = tmp_df_2[tmp_df_2['directionality'].isin(connection_type)]
    

    
    tmp_con1 = list(set(list(tmp_df_1['interactor1']) + 
                        list(tmp_df_1['interactor2'])))
    
    tmp_con2 = list(set(list(tmp_df_2['interactor1']) + 
                        list(tmp_df_2['interactor2'])))
    
    genes = list(set(tmp_con1 + tmp_con2))
    
    genes2 = [x for x in genes if x in common_idx]
    
    concat_df = concat_df.loc[genes2,:]
    
    results = calc_DEG(
        data=concat_df,
        metadata_list  = None,
        entities = 'All',
        sets = None,
        min_exp = 0,
        min_pct = 0,
        n_proc = 10,
    )
    
    results = results[results['valid_group'] == keys_list[0]]
    
    return results
    
    






