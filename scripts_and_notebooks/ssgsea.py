# -*- coding: utf-8 -*-

"""Perform ssGSEA on a given dataset using KEGG as a pathway database."""

import logging
import os

import gseapy
import pandas as pd

from pathway_forte.pathway_enrichment.functional_class import filter_gene_exp_data
from pathway_forte.constants import check_gmt_files

logger = logging.getLogger(__name__)


def run_ssgsea(filepath: str, output_dir: str):
    """Perform ssGSEA."""
    # Read data
    gene_exp = pd.read_csv(filepath, sep='\t')

    kegg_gene_set, _, _, _ = check_gmt_files()

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, merge_gene_set)

    logger.info(f'Running {filepath}')

    single_sample_gsea = gseapy.ssgsea(
        data=filtered_expression_data,
        gene_sets=kegg_gene_set,
        outdir=output_dir,  # do not write output to disk
        sample_norm_method='rank',  # choose 'custom' for your own rank list
        permutation_num=0,  # skip permutation procedure, because you don't need it
        no_plot=True,  # skip plotting to speed up
        format='png',
    )
    logger.info('Done with ssGSEA')

    single_sample_gsea.res2d.to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t')


if __name__ == '__main__':
    # TODO: change paths
    run_ssgsea('example/path', 'example/path')
    