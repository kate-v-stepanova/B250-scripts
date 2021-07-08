import os
import json
import click
import redis
import pandas as pd
import glob
import math


BASE_DIR=os.getenv('BASE_DIR')
localhost="172.22.24.88" # cable
#localhost='172.24.140.135'# wifi
remote_host="172.22.54.5" # remote
remote_host="172.22.54.5"
#remote_host="b250-bioinfo"
#remote_host="172.22.25.0"
port = "6379"
host = localhost

# this script is to be run on the cluster!

@click.group()
def cli():
    # do nothing
    pass

@cli.command()
@click.option('--remote/--local', default=False)
def check_projects(remote):
    projects = []
    for project in os.listdir(BASE_DIR):
        analysis_path = os.path.join(BASE_DIR, project, "analysis/output")
        if os.path.isdir(analysis_path):
            # if analysis dir exists
            projects.append(project)
            # demultiplexed
            # umi extracted
            # clean (without rRNA and tRNA)
            # tophat_out
            # alignments (with star)
            # rpf_5p_density
            # subsequence_data
            # rrna_fragments
            # trna_fragments
            # ucsc_tracks
            # bc_split_stats, cutadapt_plot_stats, diricore_stats
            # figures
    if remote:
        rdb = redis.StrictRedis(host=remote_host, health_check_interval=30)
    else:
        rdb = redis.StrictRedis(host=localhost, health_check_interval=30)
    rdb_projects = rdb.smembers('projects')
    rdb_projects = [p.decode('utf-8') for p in rdb_projects]
    for project in projects:
        if project not in rdb_projects:
            rdb.sadd('projects', project)
            print(project)


def check_bc_stats(project, host):
        bc_file = os.path.join(BASE_DIR, project, "analysis/output/bc_split_stats.txt")
        if os.path.exists(bc_file):
            rdb = redis.StrictRedis(host=host, health_check_interval=30)
            df = pd.read_csv(bc_file, sep="\t")
            rdb.set("bc_split_{}".format(project), df.to_msgpack())
            df = df[['Barcode', 'Count']]
            df.columns = ['sample', 'reads']
            df = df.loc[df['sample'] != 'total']
            df = df.loc[df['sample'] != 'unmatched']
            df['reads'] = df['reads'].apply(lambda x: x/1000000).round(decimals=2) # to get milion reads
            rdb.set("sample_info_{}".format(project), json.dumps(df.to_dict('recods')))


@cli.command()
@click.option("--remote/--local", default=False)
@click.argument('project')
def bc_stats(project, remote):
    if remote:
        host = remote_host
    else:
        host = localhost
    check_bc_stats(project, host)

def check_cutadapt_stats(project, host):
        cutadapt_file = os.path.join(BASE_DIR, project, "analysis/output/cutadapt_plot_stats.txt")
        if os.path.exists(cutadapt_file):
            rdb = redis.StrictRedis(host=host, health_check_interval=30)
            df = pd.read_csv(cutadapt_file, sep="\t", header=None)
            print(df)
            df = df.transpose()
            print(df)
            df.columns = df.iloc[0]
            print(df.columns)
            #df.columns = [col.strip() for col in df.columns]
            print(df.columns)
            df = df.rename({"Total reads:": "total", "With adapters:": "with_adapter", "Too short:": "too_short", "Passed filters:": "passed"}, axis='columns')
#            print(df)
            
            #df = df.loc[1]
#            df = df.loc[df['total'] != "Total reads:"]
            #import pdb; pdb.set_trace()
#            print(df)
            #df.columns = ['total', 'with_adapter', 'too_short', 'passed']
            rdb.set('cutadapt_stats_{}'.format(project), json.dumps(df.to_dict('records')))

@cli.command()
@click.option("--remote--local", default=False)
@click.argument('project')
def cutadapt_stats(project, remote):
    if remote:
        host = remote_host
    else:
        host = localhost
    check_cutadapt_stats(project, host)

def check_transcript_regions(project, host):
    print("Transcript regions")
    input_file = os.path.join(BASE_DIR, project, "analysis/output/transcript_regions/reads_per_region.tsv".format(project))
    print(input_file)
    if os.path.exists(input_file):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(input_file, sep="\t", header=0)
#        rdb.set('transcript_regions_{}'.format(project), df.to_msgpack())
        rdb.set('transcript_regions_{}'.format(project), json.dumps(df.to_dict('records')))

@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project')
def transcript_regions(project, remote):
    if remote:
        host = remote_host
    else: 
        host = localhost
    check_transcript_regions(project, host)

def check_alignment_stats(project, host):
    print('Alignment stats')
    input_file = os.path.join(BASE_DIR, project, "analysis/output/all_stats.txt")
    if os.path.exists(input_file):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(input_file, sep="\t")
        rdb.set('alignment_stats_{}'.format(project), json.dumps(df.to_dict('records')))

def check_diricore_stats(project, host):
    print("Diricore stats")
    input_file = os.path.join(BASE_DIR, project, "analysis/output/alignment_stats/diricore_stats.txt")
    if os.path.exists(input_file):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(input_file, sep="\t", header=0)
        rdb.set('diricore_stats_{}'.format(project), json.dumps(df.to_dict('records')))


@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project')
def diricore_stats(project, remote):
    host = remote_host if remote else localhost
    check_diricore_stats(project, host)
    check_alignment_stats(project, host)

def check_ucsc(project, host=localhost):
    project_dir = os.path.join(BASE_DIR, project)
    rdb = redis.StrictRedis(host=host, health_check_interval=30)
    all_files = glob.glob("{}/analysis/output/ucsc_tracks/*".format(project_dir))
    bam_types = []
    for bam in all_files:
        if os.path.isdir(bam):
            bam_types.append(os.path.basename(bam))
    for bam_type in bam_types:
        ucsc_dir = os.path.join(project_dir, "analysis/output/ucsc_tracks/{}".format(bam_type))
        if os.path.isdir(ucsc_dir):
            ucsc_annotation = os.path.join(ucsc_dir, "ucsc_track_annotation.txt")
            if os.path.isfile(ucsc_annotation):
                print("File exists: {}. Updating redis".format(ucsc_annotation))
                ucsc_link = "http://genome.ucsc.edu/s/stephz/{}_{}".format(project, bam_type)
                try:
                    rdb.sadd("ucsc_link_{}".format(project), ucsc_link)
                except:
                    existing_link = rdb.get('ucsc_link_{}'.format(project))
                    rdb.delete('ucsc_link_{}'.format(project))
                    rdb.sadd('ucsc_link_{}'.format(project), existing_link)
                    rdb.sadd('uscsc_link_{}'.format(project), ucsc_link)
                print("Added UCSC link: {}".format(ucsc_link))
        else:
            print("Skipping: {}".format(ucsc_dir))

@cli.command()
@click.option("--remote/--local", default=False)
@click.argument('project')
def ucsc(project, remote):
    host = remote_host if remote == True else localhost
    check_ucsc(project, host)


def check_ma_plot(project, host):
    rdb = redis.StrictRedis(host=host, health_check_interval=30)
    ma_files = os.path.join(BASE_DIR, project, "analysis/output/ma_plot/ma_plot_all*.tsv")
    for ma_file in glob.glob(ma_files):
        df = pd.read_csv(ma_file, sep="\t")
        contrast = os.path.basename(ma_file).replace('ma_plot_all_', '').replace('.tsv', '')
#        df1['contrast'] = contrast
        rdb.sadd('contrasts_{}'.format(project), contrast)
        print(contrast)
        rdb.set('ma_plot_all_{}_{}'.format(project, contrast), df.to_msgpack())
    ma_files = os.path.join(BASE_DIR, project, "analysis/output/ma_plot/ma_plot_coding*.tsv")
    for ma_file in glob.glob(ma_files):
        df = pd.read_csv(ma_file, sep="\t")
        contrast = os.path.basename(ma_file).replace('ma_plot_coding_', '').replace('.tsv', '')
#        df1['contrast'] = contrast
        rdb.sadd('contrasts_{}'.format(project), contrast)
        print(contrast)
        rdb.set('ma_plot_coding_{}_{}'.format(project, contrast), df.to_msgpack())
    
    #print(df)

@cli.command()
@click.option("--remote/--local", default=False)
@click.argument("project")
def ma_plot(project, remote):
    host = remote_host if remote else localhost
    check_ma_plot(project, host)

def check_fold_change(project, host):
#/{}/14592/analysis/output/rna_seq/fc_data_for_contrasts.tsv
    fc_path = os.path.join(BASE_DIR, project, "analysis/output/fold_change/FC_all_genes.tsv")
    if os.path.isfile(fc_path):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(fc_path, sep="\t")
        rdb.set("fc_all_{}".format(project), df.to_msgpack())
    fc_path = os.path.join(BASE_DIR, project, "analysis/output/fold_change/FC_coding_genes.tsv")
    if os.path.isfile(fc_path):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(fc_path, sep="\t")
        rdb.set("fc_coding_{}".format(project), df.to_msgpack())

@cli.command()
@click.option("--remote/--local", default=False)
@click.argument("project")
def fold_change(project, remote):
    host = remote_host if remote else localhost
    check_fold_change(project, host)

def check_cpm_heatmap(project, host):
    cpm_path = os.path.join(BASE_DIR, project, "analysis/output/cpm_heatmap/cpm_coding_genes.tsv")
    if os.path.isfile(cpm_path):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(cpm_path, sep="\t")
        print("Writing coding genes: {}".format(cpm_path))
        rdb.set("cpm_coding_{}".format(project), df.to_msgpack())
    cpm_path = os.path.join(BASE_DIR, project, "analysis/output/cpm_heatmap/cpm_non_coding_genes.tsv")
    if os.path.isfile(cpm_path):
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        df = pd.read_csv(cpm_path, sep="\t")
        print("Writing non-coding genes: {}".format(cpm_path))
        rdb.set('cpm_non_coding_{}'.format(project), df.to_msgpack())

@cli.command()
@click.option('--remote/--local', default=False)
@click.argument("project")
def cpm_heatmap(project, remote):
    host = remote_host if remote else localhost
    check_cpm_heatmap(project, host)

def check_reads_per_position(project_id, genome, host):
    rdb = redis.StrictRedis(host, health_check_interval=30)
    rdb.sadd('projects', project_id)

    rrna_genes = '{}/static/{}/rRNA_genes.fasta'.format(BASE_DIR, genome)
    r_df = pd.read_csv(rrna_genes,  sep='>', header=None)
    r_df = r_df.loc[~r_df[1].isna()]
    r_df = r_df[1].str.split('|', expand=True)
    r_df.columns = ['gene', 'length']
    r_dict = {}
    for i, row in r_df.iterrows():
        gene = row['gene']
        length = row['length']
        r_dict[gene] = length
    rdb.set('{}_rrna_genes'.format(project_id), json.dumps(r_dict))

    rrna_positions_dir = "analysis/output/rrna/reads_per_position"
    path = os.path.join(BASE_DIR, project_id, rrna_positions_dir, '*.tsv')
    input_files = glob.glob(path)
    if not input_files:
        print("No input files found: {}".format(path))
        return
    full_df = None
    for input_file in input_files:
        df = pd.read_csv(input_file, sep='\t')
        filename = os.path.basename(input_file)
        sample_name = filename.replace('.tsv', '')
        df['sample'] = sample_name
        if full_df is None:
            full_df = df
        else:
            full_df = full_df.append(df, ignore_index=True)
    key = "{}_reads_per_position".format(project_id)

    # save a binary string
    print("Saving key to redis: {}".format(key))
    rdb.set(key, full_df.to_msgpack())

@cli.command()
@click.argument('project_id')
@click.argument('genome', default='hg19')
@click.option('--remote/--local', default=False)
def reads_per_position(project_id, genome, remote):
    """
    Reads the data from *_reads_per_position.txt files,
    aggregates all samples into one DataFrame,
    converts it to pandas binary string and saves to Redis.
    This will update an existing record, or create a new one if does not exist.
    """
    # default host: 127.0.0.1, default port: 6379.
    # to change, use redis.StrictRedis(host=HOST, port=PORT)
    # but we are not going to change this

    if remote:
        host="172.22.54.5"
    else:
        host = localhost
    check_reads_per_position(project_id, genome, host)


def check_periodicity(project_id, bam_type, host):
    rdb = redis.StrictRedis(host, health_check_interval=30)
    rdb.sadd('projects', project_id)
    if bam_type is None:
        path = os.path.join(BASE_DIR, project_id, 'analysis/output/periodicity/all_unique/*.heatmap.tsv')
    else:
        path = os.path.join(BASE_DIR, project_id, 'analysis/output/periodicity/{}/*.heatmap.tsv'.format(bam_type))
    input_files = glob.glob(path)
    if not input_files:
        print("No input files found: {}".format(path))
        return
    full_df = None
    for input_file in input_files:
        df = pd.read_csv(input_file, sep='\t')
        samplename = os.path.basename(input_file).replace('.heatmap.tsv', '').replace('.', '_')
        df['sample'] = samplename
        if full_df is None:
            full_df = df
        else:
            full_df = full_df.append(df, ignore_index=True)
    key = "{}_periodicity_heatmap".format(project_id)

    if rdb.exists(key):
        rdb.delete(key)
    rdb.set(key, full_df.to_msgpack())

@cli.command()
@click.argument('project_id')
@click.argument('bam_type', default=None)
@click.option('--remote/--local', default=False)
def periodicity(project_id, bam_type, remote):
    if remote:
        host = "172.22.54.5"
    else:
        host = "172.22.24.88"
    
    check_periodicity(project_id, bam_type, host)


@cli.command()
@click.argument('project_id')
def copy_periodicity(project_id):
    remote = redis.StrictRedis('172.22.54.5', health_check_interval=30)
    local = redis.StrictRedis('172.22.24.88', health_check_interval=30)
    key = "{}_periodicity_heatmap".format(project_id)
    dt = pd.read_msgpack(remote.get(key))
    local.set(key, pd.to_msgpack(dt))

def get_rrna_genes(project, host):
    indir = os.path.join(BASE_DIR, project, "analysis/output/rrna/reads_per_gene")
    df = None
    for f in glob.glob(os.path.join(indir, "*.tsv")):
         samplename = os.path.basename(f)
         samplename = samplename.replace(".tsv", "")
         df1 = pd.read_csv(f, sep="\t")
         df1["sample"] = samplename
         if df is None:
              df = df1
         else:
              df = df.append(df1, ignore_index=True)
    if df is not None:
        rdb = redis.StrictRedis(host, health_check_interval=30)
        print('Saving rRNA genes')
        rdb.set("rrna_genes_{}".format(project), df.to_msgpack())
     
         
@cli.command()
@click.argument("project_id")
@click.option("--remote/--local", default=False)
def rrna_genes(project_id, remote):
    if remote:
        host = "172.22.54.5"
    else:
        host = "172.22.24.88"    
    get_rrna_genes(project_id, host)
 
def get_te1(rna_id, rp_id, host):
    rna_file = os.path.join(BASE_DIR, rna_id, "analysis/output/translational_efficiency/{}_rpkm.tsv".format(rna_id))
    rp_file = os.path.join(BASE_DIR, rp_id, "analysis/output/translational_efficiency/{}_rpkm.tsv".format(rp_id))
    if not os.path.isfile(rna_file):
        print("No file found! {}".format(rna_file))
        exit()
    if not os.path.isfile(rp_file):
        print("No file found! {}".format(rp_file))
        exit()
    
    df1 = pd.read_csv(rna_file, sep="\t")
    df2 = pd.read_csv(rp_file, sep="\t")
    
    rdb = redis.StrictRedis(host, health_check_interval=30)
    rdb.set('{}_rpkm_rna'.format(rna_id), df1.to_msgpack())
    rdb.set('{}_rpkm_rna'.format(rp_id), df1.to_msgpack())
    
    rdb.set('{}_rpkm_rp'.format(rna_id), df2.to_msgpack())
    rdb.set('{}_rpkm_rp'.format(rp_id), df2.to_msgpack())

def get_te(rna_id, rp_id, host):
    project_id = rna_id
    for f in glob.glob(os.path.join(BASE_DIR, rna_id, "analysis/output/ribo_diff/te/all_*_plot_data.tsv")):
        df = pd.read_csv(f, sep="\t")
        rdb = redis.StrictRedis(host, health_check_interval=30)
        contrast = os.path.basename(f).replace('all_', '').replace('_plot_data.tsv', '')
        print("Saving data for contrast: {}".format(contrast))
        rdb.set('te_{}_{}'.format(rna_id, contrast), df.to_msgpack())
        if rp_id is not None:
            rdb.set('te_{}_{}'.format(rp_id, contrast), df.to_msgpack())


@cli.command()
@click.argument("rna_id")
@click.argument("rp_id", required=False)
@click.option("--remote/--local", default=False)
def te(rna_id, remote, rp_id=None):
    if remote:
        host = remote_host
    else:
        host = localhost
    get_te(rna_id, rp_id, host)

def check_snoRNAs(project_id, host):
    infile = os.path.join(BASE_DIR, project_id, "analysis/output/snoRNAs/barplot.txt")
    if os.path.isfile(infile):
        print("Writing snoRNA stats")
        df = pd.read_csv(infile, sep="\t")
        rdb = redis.StrictRedis(host=host, health_check_interval=30)
        rdb.set('{}_snoRNAs'.format(project_id), df.to_msgpack())

@cli.command()
@click.argument("project_id")
@click.option("--remote/--local", default=False)
def sno_rna(project_id, remote):
    host = remote_host if remote else localhost
    check_snoRNAs(project_id, host)

def get_all_stats(project, host=localhost):
    print("Getting all stats for {}".format(project))
    project_dir = os.path.join(BASE_DIR, project)
    analysis_dir = os.path.join(project_dir, "analysis/output")
    rdb = redis.StrictRedis(host=host, health_check_interval=30)
    if os.path.isdir(analysis_dir):
        rdb.sadd("projects", project)
        check_cutadapt_stats(project, host)
        check_bc_stats(project, host)
        check_transcript_regions(project, host)
        check_diricore_stats(project, host)
        check_ma_plot(project, host)
        check_ucsc(project, host)
#        check_reads_per_position(project, 'hg19', host)
        check_periodicity(project, None, host)
        get_rrna_genes(project, host)
        get_te(project, None, host)
        check_snoRNAs(project, host)
        #upload_volcano_plot(project, host)
        # upload_alignments(project, host)
        upload_psites(project, host)
        


@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project')
def all_stats(project, remote):
    host = remote_host if remote else localhost
    get_all_stats(project, host)

@cli.command()
@click.option('--remote/--local', default=False)
def all_projects(remote):
    if remote:
        host = remote_host
    else:
        host = localhost

    projects = os.listdir(BASE_DIR)
    rdb = redis.StrictRedis(host=host, health_check_interval=30)
    #projects = rdb.smembers('projects')
#    projects = [p.decode('utf-8') for p in projects]
    for project in projects:
        if os.path.isdir(os.path.join(BASE_DIR, project, "analysis/input")):
            rdb.sadd('projects', project)
            get_all_stats(project, host)

def clear_project_info(host, project):
    rdb = redis.StrictRedis(host=host, health_check_interval=30)
    
@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project_id')
def alignments(remote, project_id):
    if remote:
        host = remote_host
    else:
        host = localhost
    if project_id == 'all':
        indir = BASE_DIR
        for f in os.listdir(indir):
            if os.path.isdir(os.path.join(indir, f)):
                proj = os.path.basename(f)
                upload_alignments(proj, host)
    else:
        upload_alignments(project_id, host)
    

def upload_alignments(project_id, host):
    rdb = redis.StrictRedis(host=host, health_check_interval=30)
    indir = "{}/{}/analysis/output/ext_diricore".format(BASE_DIR, project_id)
    if not os.path.exists(indir):
        return
    bam_types = []
    for f in os.listdir(indir):
        if os.path.isdir(os.path.join(indir, f)):
            bam_types.append(os.path.basename(f))
    for bam_type in bam_types:
        indir = "{}/{}/analysis/output/ext_diricore/{}/tsv".format(BASE_DIR, project_id, bam_type)
        trans_file = "{}/static/hg19/transcriptLength.txt".format(BASE_DIR)
        trans = pd.read_csv(trans_file, sep="\t")
        rdb_trans = rdb.get('transcriptLength.txt')
        if rdb_trans is None:
            rdb.set('transcriptLength.txt', trans.to_msgpack())
        list_of_genes = []
        for f in glob.glob(os.path.join(indir, "*.tsv")):
            size = os.path.getsize(f)
                    
            print('Uploading alignments for : {}/{}/{}'.format(project_id, bam_type, os.path.basename(f)))
            sample = os.path.basename(f).replace('.tsv', '')
            df = pd.read_csv(f, sep="\t", header=None, names=["transcript", 'start', 'seq'])
            #df['end'] = df['start'] + df['seq'].str.len()
            df1 = pd.merge(df, trans, on="transcript", how="left")
            print("Writing {}".format(f))
            rdb.set('alignment__{}__{}__{}'.format(project_id, bam_type, sample), df.to_msgpack())
            rdb.sadd('bam_types_{}'.format(project_id), bam_type)
            list_of_genes += df1['gene_name'].unique().tolist()
        if list_of_genes:
            rdb.sadd('genes_{}'.format(project_id), *set(list_of_genes))

def upload_counts(project_id, bam_type):
    indir = "{}/{}/analysis/output/alignments/reads_per_gene/normalized/{}".format(BASE_DIR, project_id, bam_type)
    full_df = None
    pro1000 = "{}/static/hg19/pro_asp_leu/pro_trans_top1000.tsv".format(BASE_DIR)
    leu1000 = "{}/static/hg19/pro_asp_leu/leu_trans_top1000.tsv".format(BASE_DIR)
    asp1000 = "{}/static/hg19/pro_asp_leu/asp_trans_top1000.tsv".format(BASE_DIR)
    pro_df = pd.read_csv(pro1000)
    leu_df = pd.read_csv(leu1000)
    asp_df = pd.read_csv(asp1000)
    for f in glob.glob("{}/*.tsv".format(bam_type)):
        df = pd.read_csv(f, sep="\t")

@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project_id')
def counts(remote, project_id):
    indir = "{}/{}/analysis/output/alignments/reads_per_gene/normalized/".format(BASE_DIR, project_id)
    dir_content = glob.glob("{}/*".format(indir))
    bam_types = []
    for f in dir_content:
        if os.path.isdir(f):
             bam_types.append(os.path.basename(f))
    for bam_type in bam_types:
        upload_counts(project_id, bam_type)

def upload_volcano_plot(project_id, host):

    rdb = redis.StrictRedis(host, health_check_interval=30)
    indir = "{}/{}/analysis/output/diff_expr".format(BASE_DIR, project_id)
    contrasts = []
    for f in glob.glob(os.path.join(indir, "*.tsv")):
        contrast = os.path.basename(f).replace('diff_expr_', '').replace('.tsv', '')
        contrasts.append(contrast)
        df = pd.read_csv(f, sep="\t")
        df = df[['log2FoldChange', 'pvalue', 'padj', 'gene']]
        print('Uploading data for contrast: {}'.format(contrast))
        rdb.set('volcano_{}_{}'.format(project_id, contrast), json.dumps(df.to_dict('records')))
        rdb.sadd('contrasts_{}'.format(project_id), contrast)
    asp_top200 = "{}/static/hg19/pro_asp_leu/top200/asp_top200.tsv".format(BASE_DIR)
    pro_top200 = "{}/static/hg19/pro_asp_leu/top200/pro_top200.tsv".format(BASE_DIR)
    asp_df = pd.read_csv(asp_top200, sep="\t", header=None, names=['gene', 'ratio'])
    asp_df = asp_df.drop('ratio', axis=1)
    pro_df = pd.read_csv(pro_top200, sep="\t", header=None, names=['gene', 'ratio'])
    pro_df = pro_df.drop('ratio', axis=1)
    rdb.sadd('asp200', asp_df['gene'].tolist())
    rdb.sadd('pro200', pro_df['gene'].tolist())
        

@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project_id')
def volcano_plot(remote, project_id):
    host = remote_host if remote else localhost
    asp_top200 = "/{}/static/hg19/pro_asp_leu/top200/asp_top200.tsv".format(BASE_DIR)
    pro_top200 = "/{}/static/hg19/pro_asp_leu/top200/pro_top200.tsv".format(BASE_DIR)
    asp_df = pd.read_csv(asp_top200, sep="\t", header=None, names=['gene', 'ratio'])
    asp_df = asp_df.drop('ratio', axis=1)
    pro_df = pd.read_csv(pro_top200, sep="\t", header=None, names=['gene', 'ratio'])
    pro_df = pro_df.drop('ratio', axis=1)
    print(pro_df['gene'].tolist())
    print(asp_df['gene'].tolist()) 
    upload_volcano_plot(project_id, host)    

def upload_contrasts(project_id, host):
    contrast_file = "/{}/{}/analysis/input/metadata/rpf_density_contrasts.tsv".format(BASE_DIR, project_id)
    contrast_df = pd.read_csv(contrast_file, sep="\t", usecols=[0,1], names=['sample', 'contrast'])
    contrast_df['contrast'] = contrast_df['sample'] + '__vs__' + contrast_df['contrast']
    rdb = redis.StrictRedis(host, health_check_interval=30)
    contrasts = contrast_df['contrast'].tolist()
    for c in contrasts:
        rdb.sadd('contrasts_{}'.format(project_id), c)


def upload_psites(project_id, host, bam_type='all_unique'):
    psites_file = "/{}/{}/analysis/output/ext_diricore/{}/psites/psites.tsv".format(BASE_DIR, project_id, bam_type)
    if not os.path.isfile(psites_file):
        print('File not found: {}'.format(psites_file))
        return
    asites_file = "/{}/{}/analysis/output/ext_diricore/{}/psites/asites.tsv".format(BASE_DIR, project_id, bam_type)
    esites_file = "/{}/{}/analysis/output/ext_diricore/{}/psites/esites.tsv".format(BASE_DIR, project_id, bam_type)
    psites_df = pd.read_csv(psites_file, sep="\t")
    asites_df = pd.read_csv(asites_file, sep="\t")
    esites_df = pd.read_csv(esites_file, sep="\t")
    print('Redis: {}'.format(host))
    rdb = redis.StrictRedis(host, health_check_interval=30)
    print('Uploading psites: {}'.format(psites_file))
    rdb.set('psites_{}'.format(project_id), json.dumps(psites_df.to_dict('records')))
    print('Uploading asites: {}'.format(asites_file))
    rdb.set('asites_{}'.format(project_id), json.dumps(asites_df.to_dict('records')))
    print('Uploading esites: {}'.format(esites_file))
    rdb.set('esites_{}'.format(project_id), json.dumps(esites_df.to_dict('records')))
    
@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project_id')
@click.argument('bam_type')
def psites(remote, project_id, bam_type):
    host = remote_host if remote else localhost
    upload_psites(project_id, host, bam_type=bam_type)
    upload_contrasts(project_id, host)

@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project_id')
def contrasts(remote, project_id):
    host = remote_host if remote else localhost
    upload_contrasts(project_id, host)

def upload_psite_dotplots(project_id, host, bam_type='all_unique'):
    rdb = redis.StrictRedis(host, health_check_interval=30)
    psites_file = "/{}/{}//analysis/output/ext_diricore/{}/psites_dotplot/*/*psite_dotplot.tsv".format(BASE_DIR, project_id, bam_type)
    asites_file =  "/{}/{}//analysis/output/ext_diricore/{}/psites_dotplot/*/*asite_dotplot.tsv".format(BASE_DIR, project_id, bam_type)
    esites_file =  "/{}/{}//analysis/output/ext_diricore/{}/psites_dotplot/*/*esite_dotplot.tsv".format(BASE_DIR, project_id, bam_type)
    aas = [] 
    # P-sites
    for f in glob.glob(psites_file):
        print(f)
        aa = os.path.basename(f).replace('_psite_dotplot.tsv', '')
        aas.append(aa)
        df = pd.read_csv(f, sep="\t")
        cols = [c for c in df.columns if 'tpm_' not in c and 'rpkm_' not in c and c not in ['gene', 'Aa', 'codon']]
        df1 = df[cols]
        # take only codons with min 4 reads per codon in average per sample
        df1 = df1.loc[df1.sum(axis=1) >= 4 * len(cols)]
        df = df.loc[df1.index]
        key = 'psite_dotplot_{}_{}'.format(project_id, aa)
        rdb.set(key, json.dumps(df.to_dict('records')))
        # top genes for AA
        aa_file = "/{}/static/hg19/pro_asp_leu/top1000/{}_1000.tsv".format(BASE_DIR, aa)
        if os.path.isfile(aa_file):
             aa_df = pd.read_csv(aa_file, sep="\t", usecols=[0], names=['gene'], header=None)
             key = '{}_top1000'.format(aa)
             print(aa_df)
             rdb.set(key, json.dumps(aa_df.to_dict('list')))
    # A-sites
    for f in glob.glob(asites_file):
        print(f)
        aa = os.path.basename(f).replace('_asite_dotplot.tsv', '')
        df = pd.read_csv(f, sep="\t")
        cols = [c for c in df.columns if 'tpm_' not in c and 'rpkm_' not in c and c not in ['gene', 'Aa', 'codon']]
        df1 = df[cols]
        # take only codons with min 4 reads per codon in average per sample
        df1 = df1.loc[df1.sum(axis=1) >= 4 * len(cols)]
        df = df.loc[df1.index]
        key = 'asite_dotplot_{}_{}'.format(project_id, aa)
        rdb.set(key, json.dumps(df.to_dict('records')))
    # E-sites
    for f in glob.glob(esites_file):
        print(f)
        aa = os.path.basename(f).replace('_esite_dotplot.tsv', '')
        df = pd.read_csv(f, sep="\t")
        cols = [c for c in df.columns if 'tpm_' not in c and 'rpkm_' not in c and c not in ['gene', 'Aa', 'codon']]
        df1 = df[cols]
        # take only codons with min 4 reads per codon in average per sample
        df1 = df1.loc[df1.sum(axis=1) >= 4 * len(cols)]
        df = df.loc[df1.index]
        key = 'esite_dotplot_{}_{}'.format(project_id, aa)
        rdb.set(key, json.dumps(df.to_dict('records')))
    rdb.sadd('aa_dotplot_{}'.format(project_id), *set(aas))
    contrasts = "/{}/{}//analysis/input/metadata/rpf_density_contrasts.tsv".format(BASE_DIR, project_id)
    c_df = pd.read_csv(contrasts, sep="\t", header=None, usecols=[0,1], names=['sample', 'control'])
    c = []
    for i, row in c_df.iterrows():
        sample = row['sample']
        control = row['control']
        contrast = '{}__vs__{}'.format(sample, control)
        c.append(contrast)
    
    contrasts = rdb.sadd('contrasts_{}'.format(project_id), *set(c))
    



@cli.command()
@click.option('--remote/--local', default=False)
@click.argument('project_id')
@click.argument('bam_type')
def psite_dotplots(remote, project_id, bam_type):
    host = remote_host if remote else localhost
    upload_psite_dotplots(project_id, host, bam_type)

if __name__ == '__main__':
    cli()
