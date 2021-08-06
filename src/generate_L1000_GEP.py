from cmapPy.pandasGEXpress.parse import parse
import argparse
import pandas as pd
import numpy as np
import os
import sys


def downloadFromGEO(filename, url):
    print('downloading {} from the GEO website...'.format(filename))
    from urllib import request
    with request.urlopen(url) as response, open(filename, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    if('gctx' in filename):
        import gzip
        import shutil
        with gzip.open(filename, 'rb') as f_in:
            with open(filename.rsplit('.',1)[0], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate bulk profiles of the specified cell type from the LINCS L1000 database.")

    parser.add_argument("-o", "--outdir", default='./LINCS/', help="path to output directory, default='./'")
    parser.add_argument("-c", "--celltype", default='MCF7', help='Cell Line name, default: MCF7. Options: A375|A549|HCC515|HEPG2|MCF7|PC3|VCAP|HT29')
    parser.add_argument("--inst", type=str, default='GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz', help="inst_info file (.txt.gz)")
    parser.add_argument("--gctx", type=str, default='GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx', help="LINCS L1000 level 3 GEPs (.gctx)")
    parser.add_argument("--gene", type=str, default='GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz', help="gene_info file (.txt.gz)")
    
    args = parser.parse_args()

    cell_types = ['A375','A549','HCC515','HEPG2','MCF7','PC3','VCAP','HT29']

    if not args.celltype in cell_types:
        sys.exit("Unacceptable cell type.")
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    if not (os.path.isfile(args.inst) and args.inst.endswith('.gz')):
        if args.inst == 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz':
            downloadFromGEO('GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5Finst%5Finfo%5F2017%2D03%2D06%2Etxt%2Egz')
        else:
            sys.exit("The inst_info file does not exist or is not .gz file.")
    if not (os.path.isfile(args.gctx) and args.gctx.endswith('.gctx')):
        if args.gctx == 'GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx':
            downloadFromGEO('GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5FLevel3%5FINF%5Fmlr12k%5Fn345976x12328%5F2017%2D03%2D06%2Egctx%2Egz')
            
        else:
            sys.exit("The gctx file does not exist or is not gz file.")
    if not (os.path.isfile(args.gene) and args.gene.endswith('.gz')):
        if args.gene == 'GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz':
            downloadFromGEO('GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138%5FBroad%5FLINCS%5Fgene%5Finfo%5F2017%2D03%2D06%2Etxt%2Egz')
        else:
            sys.exit("The gene_info file does not exist or is not .gz file.")


    # read inst_info and gene_info
    inst_info = pd.read_csv(args.inst, sep='\t', compression='gzip')
    sig_info = pd.read_csv(args.gene, sep='\t', usecols=['pr_gene_id','pr_gene_symbol'], compression='gzip')

    cell = args.celltype
    
    # select instance ids for a specific cell type
    inst_ids = inst_info['inst_id'][inst_info['cell_id'] == cell]
    # read gctx
    gctoo = parse(args.gctx, cid=inst_ids)
    gctoo.data_df.index = gctoo.data_df.index.astype(int)
    # covert rowids to gene names
    named_df = pd.merge(gctoo.data_df, sig_info, left_index=True, right_on=['pr_gene_id'], validate='1:1')
    average_df = named_df.groupby('pr_gene_symbol').mean().dropna().drop(labels='pr_gene_id',axis=1)
    # reverse to non-log for CIBERSORTx
    exp_df = 2**average_df
    exp_df.to_csv('{}/LINCS_L1000_GEP_{}.csv'.format(args.outdir, cell), sep='\t')