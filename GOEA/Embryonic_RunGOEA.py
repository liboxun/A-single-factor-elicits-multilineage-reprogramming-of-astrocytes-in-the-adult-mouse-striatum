from __future__ import print_function

# coding: utf-8

# In[ ]:


# Run a Gene Ontology Enrichment Analysis (GOEA)

import sys
import pickle
print(sys.executable)
print(sys.path)
print(sys.version)

from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
# import scanpy as sc
# import harmonypy as hm

plt.rcParams['figure.figsize']=(8,8) #rescale figures
# sc.settings.verbosity = 3

print('numpy', np.__version__)
print('pandas', pd.__version__)
# print('scipy', scipy.__version__)
# print('sklearn', sklearn.__version__)
print('statsmodels', sm.__version__)
print('matplotlib', mpl.__version__)
print('seaborn', sns.__version__)

# sc.logging.print_versions()

# !date +%F

## 1. Download Ontologies and Associations

### 1a. Download Ontologies, if necessary

# Get http://geneontology.org/ontology/go-basic.obo
from goatools.base import download_go_basic_obo
obo_fname = download_go_basic_obo()

### 1b. Download Associations, if necessary

# Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
from goatools.base import download_ncbi_associations
fin_gene2go = download_ncbi_associations()

## 2. Load Ontologies, Associations and Background gene set 

### 2a. Load Ontologies

from goatools.obo_parser import GODag

obodag = GODag("go-basic.obo")

### 2b. Load Associations

from goatools.anno.genetogo_reader import Gene2GoReader

# Read NCBI's gene2go. Store annotations in a list of namedtuples
objanno = Gene2GoReader(fin_gene2go, taxids=[10090]) # Taxid 10090 is Mus musculus

# Get namespace2association where:
#    namespace is:
#        BP: biological_process               
#        MF: molecular_function
#        CC: cellular_component
#    assocation is a dict:
#        key: NCBI GeneID
#        value: A set of GO IDs associated with that gene
ns2assoc = objanno.get_ns2assc()

for nspc, id2gos in ns2assoc.items():
    print("{NS} {N:,} annotated mouse genes".format(NS=nspc, N=len(id2gos)))

### 2c. Load Background gene set

from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
# from goatools.test_data.genes_NCBI_10090_All import GENEID2NT as GeneID2nt_mus

len(GeneID2nt_mus.keys())

## 3. Initialize a GOEA object

from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_mus.keys(), # List of mouse protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method


## 4. Read study genes

#### Making a reference dictionary to convert symbols of my genes to gene IDs

GeneID2Symbol_dict = {}

for GeneID in list(GeneID2nt_mus.keys()):
    Symbols = []
    Symbols.append(GeneID2nt_mus[GeneID].Symbol)
    if not GeneID2nt_mus[GeneID].Aliases == '':
        for alias in GeneID2nt_mus[GeneID].Aliases.split(','):
            Symbols.append(alias.strip())
    GeneID2Symbol_dict[GeneID] = Symbols

##### Example search of a gene symbol

gene_search = 'Ascl1'

for GeneID, names in GeneID2Symbol_dict.items():
    if gene_search in names:
        print(GeneID)

#### Loading my genes and convert symbols to IDs

DE_dir = '/project/GCRB/Hon_lab/s418610/Projects/03.Invivo_neuronal_reprogramming/analysis/Cluster_with_10X_20k/data/10x_remapped/Revision/DE_analysis/'

input_path = DE_dir + 'embryonic_neurogenesis_test_list_collection_100iter'
with open(input_path, 'rb') as input_file:
 
    test_list_collection = pickle.load(input_file)

assert len(test_list_collection) == 100

## 5. Run Gene Ontology Enrichment Analysis (GOEA)
# You may choose to keep all results or just the significant results. In this example, we choose to keep only the significant results.

goea_collection_experiment = []
goea_collection_public = []

# test_list_all = test_list_collection[0]
for iteration, test_list_all in enumerate(test_list_collection): 
    print('\n\n\n\n')
    print('Current iteration: ', iteration, '\n')

    #### Examine genes enriched in experimental cells

    #### Remove DEGs that are common to all clusters (see DE analysis jupyter notebook)

    #### For each iteration of DE analysis, run GOEA for one time
    
    goea_experiment_up_all = []
    for cluster, test_cluster in test_list_all:
        experiment_up_genes = test_cluster.summary()[np.logical_and(test_cluster.summary()['qval'] < 1e-2, test_cluster.summary()['log2fc'] > 1)]
        geneids_study = []
        # Remove commonly differentially expressed genes (see DE analysis jupyter notebook)
        for gene_search in experiment_up_genes['gene']:
            if gene_search in ['Lars2', 'Gm42418', 'AY036118']: # The 7's and 6's
                continue
            for GeneID, names in GeneID2Symbol_dict.items():
                if gene_search in names:
        #            print(GeneID)
                    geneids_study.append(GeneID)
                    break

        goea_results_all = goeaobj.run_study(geneids_study)
        # goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
        goea_experiment_up_all.append((cluster, experiment_up_genes, geneids_study, goea_results_all))


    #### Examine genes enriched in experimental cells

    #### Remove DEGs that are common to all clusters (see DE analysis jupyter notebook)

    #### For each iteration of DE analysis, run GOEA for one time

    goea_public_up_all = []
    for cluster, test_cluster in test_list_all:
        public_up_genes = test_cluster.summary()[np.logical_and(test_cluster.summary()['qval'] < 1e-2, test_cluster.summary()['log2fc'] < -1)]
        geneids_study = []
        # Remove commonly differentially expressed genes (see DE analysis jupyter notebook)
        for gene_search in public_up_genes['gene']:
            if gene_search in ['Rpl10', 'Rps27', 'Rps12', 'Rpl9', 'Rpl27']:
                continue
            for GeneID, names in GeneID2Symbol_dict.items():
                if gene_search in names:
        #            print(GeneID)
                    geneids_study.append(GeneID)
                    break

        goea_results_all = goeaobj.run_study(geneids_study)
        # goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    #    goea_public_up_all.append((cluster, public_up_genes, geneids_study, goea_results_sig))
        goea_public_up_all.append((cluster, public_up_genes, geneids_study, goea_results_all)) # Store all test results, as GO terms of interest might not always be significant 
                                                                                                # in all cell types

    # len(goea_public_up_all)
    goea_collection_experiment.append((iteration, goea_experiment_up_all))
    goea_collection_public.append((iteration, goea_public_up_all))

## 6. Plot enriched GO terms of interest as a function of reprogramming progression

#### Extracting results

del test_list_collection

goea_DF_collection_public = []
for iteration, goea_public_up_all in goea_collection_public:

    goea_DF_list_public = []
    for cluster, _, _, goea_results_sig in goea_public_up_all:

        goID = []
        goName = []
        goNamespace = []
        goEnrichment = []
        goAdjP = []
        study_count = []
        study_n = []
        pop_count = []
        pop_n = []
        for goea_result in goea_results_sig:
            goID.append(goea_result.goterm.id)
            goName.append(goea_result.goterm.name)
            goNamespace.append(goea_result.goterm.namespace)
            goEnrichment.append(goea_result.enrichment)
            goAdjP.append(goea_result.p_fdr_bh)
            study_count.append(goea_result.study_count)
            study_n.append(goea_result.study_n)
            pop_count.append(goea_result.pop_count)
            pop_n.append(goea_result.pop_n)

        goea_summary = pd.DataFrame(data={'goID': goID,
                                         'goName': goName, 
                                         'goNamespace': goNamespace, 
                                         'goEnrichment': goEnrichment, 
                                         'goAdjP': goAdjP, 
                                         'study_count': study_count, 
                                         'study_n': study_n, 
                                         'pop_count': pop_count, 
                                         'pop_n': pop_n})

        # goea_summary.sort_values(by='goAdjP', ascending=True, inplace=True)

        # goea_summary = goea_summary[goea_summary['goEnrichment']=='e']

        goea_summary_MF = goea_summary[goea_summary['goNamespace']=='molecular_function']
        goea_summary_BP = goea_summary[goea_summary['goNamespace']=='biological_process']
        goea_summary_CC = goea_summary[goea_summary['goNamespace']=='cellular_component']

        goea_DF_list_public.append((cluster, goea_summary_MF, goea_summary_BP, goea_summary_CC))
        
    goea_DF_collection_public.append((iteration, goea_DF_list_public))

len(goea_DF_collection_public)

goea_DF_collection_experiment = []
for iteration, goea_experiment_up_all in goea_collection_experiment:

    goea_DF_list_experiment = []
    for cluster, _, _, goea_results_sig in goea_experiment_up_all:

        goID = []
        goName = []
        goNamespace = []
        goEnrichment = []
        goAdjP = []
        study_count = []
        study_n = []
        pop_count = []
        pop_n = []
        for goea_result in goea_results_sig:
            goID.append(goea_result.goterm.id)
            goName.append(goea_result.goterm.name)
            goNamespace.append(goea_result.goterm.namespace)
            goEnrichment.append(goea_result.enrichment)
            goAdjP.append(goea_result.p_fdr_bh)
            study_count.append(goea_result.study_count)
            study_n.append(goea_result.study_n)
            pop_count.append(goea_result.pop_count)
            pop_n.append(goea_result.pop_n)

        goea_summary = pd.DataFrame(data={'goID': goID,
                                         'goName': goName, 
                                         'goNamespace': goNamespace, 
                                         'goEnrichment': goEnrichment, 
                                         'goAdjP': goAdjP, 
                                         'study_count': study_count, 
                                         'study_n': study_n, 
                                         'pop_count': pop_count, 
                                         'pop_n': pop_n})

        # goea_summary.sort_values(by='goAdjP', ascending=True, inplace=True)

        # goea_summary = goea_summary[goea_summary['goEnrichment']=='e']

        goea_summary_MF = goea_summary[goea_summary['goNamespace']=='molecular_function']
        goea_summary_BP = goea_summary[goea_summary['goNamespace']=='biological_process']
        goea_summary_CC = goea_summary[goea_summary['goNamespace']=='cellular_component']

        goea_DF_list_experiment.append((cluster, goea_summary_MF, goea_summary_BP, goea_summary_CC))
        
    goea_DF_collection_experiment.append((iteration, goea_DF_list_experiment))

len(goea_DF_collection_experiment)

##### Write and read data

GOEA_dir = '/project/GCRB/Hon_lab/s418610/Projects/03.Invivo_neuronal_reprogramming/analysis/Cluster_with_10X_20k/data/10x_remapped/Revision/GOEA/'

output_path = GOEA_dir + 'embryonic_neurogenesis_goea_DF_collection_public_100iter'

with open(output_path, 'wb') as output_file:
 
    pickle.dump(goea_DF_collection_public, output_file)

output_path = GOEA_dir + 'embryonic_neurogenesis_goea_DF_collection_experiment_100iter'

with open(output_path, 'wb') as output_file:
 
    pickle.dump(goea_DF_collection_experiment, output_file)

