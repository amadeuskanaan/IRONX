import json
import urllib2
import os
from tables import open_file
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import preprocessing
#from alleninf.api import get_probes_from_genes
from alleninf.data import get_values_at_locations
from alleninf.api import get_mni_coordinates_from_wells#
from alleninf.analysis import fixed_effects, approximate_random_effects, bayesian_random_effects
import statsmodels.formula.api as smf
import math
#from variables.variables import *
from AHBA.genesets_iron import *
from sklearn.decomposition import TruncatedSVD, PCA
api_url        = "http://api.brain-map.org/api/v2/data/query.json"


ahba_dir     = '/Users/kanaaax/Google Drive/TS-EUROTRAIN/RESULTS_QSMv3/dataframes/AHBA/'


def get_probes_from_genes(gene_names):
    if not isinstance(gene_names, list):
        gene_names = [gene_names]
    # in case there are white spaces in gene names
    gene_names = ["'%s'" % gene_name for gene_name in gene_names]

    api_query = "?criteria=model::Probe"
    api_query += ",rma::criteria,[probe_type$eq'DNA']"
    api_query += ",products[abbreviation$eq'HumanMA']"
    api_query += ",gene[acronym$eq%s]" % (','.join(gene_names))
    api_query += ",rma::options[only$eq'probes.id','name']"

    data = json.load(urllib2.urlopen(api_url + api_query))

    d = {probe['id']: probe['name'] for probe in data['msg']}

    if not d:
        print 'Gene %s not available'%gene_name#
        # raise Exception("Could not find any probes for %s gene. Check "
        #                "http://help.brain-map.org/download/attachments/2818165/HBA_ISH_GeneList.pdf?version=1&modificationDate=1348783035873 "
        #                "for list of available genes." % gene_name)
        pass
    else:
        return d


def return_probe_expression(gene_probes_dict, geneset_name):
    dfs = []
    genes = gene_probes_dict.keys()

    if not os.path.isfile(os.path.join(ahba_dir, 'AHBA/AHBA_%s.csv' % geneset_name)):

        print 'Fetching normalized gene expression values for:', genes
        print ''
        for gene in genes:
            probe_ids = ["'%s'" % probe_id for probe_id in gene_probes_dict[gene].keys()]
            print 'Probe IDs for Gene %s = %s' % (gene, probe_ids)

            api_query = api_url + "?criteria=service::human_microarray_expression[probes$in%s]" % (','.join(probe_ids))
            data = json.load(urllib2.urlopen(api_query))
            print api_query

            cols = ['top_struct', 'struct_name', 'struct_id', 'donor_names', 'coords_native']
            probe_cols = ['%s_' % gene + str(i) for i in gene_probes_dict[gene].values()]
            cols = cols + probe_cols
            well_ids = [str(sample["sample"]["well"]) for sample in data["msg"]["samples"]]
            df = pd.DataFrame(index=well_ids, columns=cols)

            df['top_struct'] = [sample["top_level_structure"]["name"] for sample in data["msg"]["samples"]]
            df['struct_id'] = [sample["structure"]["id"] for sample in data["msg"]["samples"]]
            df['struct_name'] = [sample["structure"]["name"] for sample in data["msg"]["samples"]]
            df['donor_names'] = [sample["donor"]["name"] for sample in data["msg"]["samples"]]
            df['coords_native'] = [sample["sample"]["mri"] for sample in data["msg"]["samples"]]

            for i, probe_id in enumerate(gene_probes_dict[gene].values()):
                df['%s_%s' % (gene, probe_id)] = [float(expression_value) for expression_value in
                                                  data["msg"]["probes"][i]["expression_level"]]

            dfs.append(df)

        # concat all probe expression dataframes
        df = pd.concat(dfs, axis=1).T.groupby(level=0).first().T
        df.to_csv(os.path.join(ahba_dir, 'PROBES/PROBES_%s.csv' % geneset_name))

        # decompose probe expression values
        all_probes = ['%s_'%gene + str(i) for gene in gene_probes_dict.keys() for i in gene_probes_dict[gene].values()]

        # Need to simplify Dataframe by averaging all probes for every gene
        #split dataframe in probes and metadata
        df_probes     = df.iloc[:, :-5]
        df_metadata   = df.iloc[:, -5:]
        probes_unique = [s.split('_')[0] for s in df_probes.T.index.values]
        df_probes_mu  = df_probes.astype(float).T.groupby(probes_unique).mean().T

        df = pd.concat([df_probes_mu, df_metadata], axis=1)
        if len(set(probes_unique)) > 1:
            df['Mean'] = df[list(set(probes_unique))].mean(axis=1)
            df['Median'] = df[list(set(probes_unique))].median(axis=1)
            pca_genes = PCA()
            pca_genes.fit(np.array(np.asarray([df[gene] for gene in genes])))
            pca_probes = PCA()
            pca_probes.fit(np.array(np.asarray([df_probes[probe] for probe in all_probes])))
            try:
                df['SVD1g'] = pca_genes.components_[0, :]
                df['SVD1p'] = pca_probes.components_[0, :]
            except:
                print 'No SVD1'
            try:
                df['SVD2g'] = pca_genes.components_[1, :]
                df['SVD2p'] = pca_probes.components_[1, :]
            except:
                print 'No SVD2'
            try:
                df['SVD3g'] = pca_genes.components_[2, :]
                df['SVD3p'] = pca_probes.components_[2, :]
            except:
                print 'No SVD3'

            # print 'PC explained variance:', pca.explained_variance_ratio_
            #df['PC_EV'] = pca.explained_variance_ratio_[0]#, pca.explained_variance_ratio_[1], pca.explained_variance_ratio_[2],


            #package_directory = '/Users/kanaan/SCR/Github/alleninf/alleninf'
        package_directory = '/scr/malta1/Software/anaconda/envs/awesome/lib/python2.7/site-packages/alleninf'
        package_directory = '/Users/kanaaax/anaconda2/lib/python2.7/site-packages/alleninf'
        mni = pd.read_csv(os.path.join(package_directory, "data", "corrected_mni_coordinates.csv"),
                          header=0, index_col=0)
        mni.index = mni.index.map(unicode)
        df_concat = pd.concat([df, mni], axis=1).to_csv(os.path.join(ahba_dir, 'AHBA_%s.csv' % geneset_name))

        return df_concat
    return pd.read_csv(os.path.join(ahba_dir, 'AHBA/AHBA_%s.csv' % geneset_name), index_col=0)

def get_expression_df(genes, geneset_name):
    gene_probes = {}

    print '##########################################################################################################'
    print '##########################    Downloading AHBA for data for geneset %s      ##############################'%geneset_name
    print '##########################################################################################################'

    for gene in genes:
        probes  = get_probes_from_genes(gene)
        if probes:
            gene_probes[gene] = probes
        else:
            pass #print 'Gene %s has no probes' %gene

    return_probe_expression(gene_probes, geneset_name)
    #return df


genesets = ['IRON', 'IRON_D','ANMC', 'GLU', 'GABA', 'FTH', 'TF','FTH_ALL', 'FTL_ALL', 'FERRITIN']


# get_expression_df(IRON_HOMEOSTASIS.keys()  , 'IRON_HOMEOSTASIS')
# get_expression_df(IRON_D.keys()            , 'IRON_D')
# get_expression_df(IRON_TRANSPORT2          , 'IRON_TRANSPORT2')
# get_expression_df(FERRITIN              , 'FERRITIN')
get_expression_df(IRON_SULFUR, 'IRON_SULFUR')

