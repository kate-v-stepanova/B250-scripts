""" PSL - PSL format from the BLAT program """

def calcPercentIdentity(row):
    ''' Return number of non-identical matches. 
    Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
    qalisize = row['qspan']
    alisize = min(qalisize, row['tspan'])
    millibad = 0
    if alisize <= 0: return 0
    sizediff = alisize - row['tspan']
    if sizediff < 0:
        if ismrna:
            sizediff = 0
        else:
            sizediff = -sizediff
    insertfactor = row['qnuminsert']
    if not ismrna: insertfactor += row['tnuminsert']
    total = row['matches'] + row['repmatches'] + row['mismatches']
    if total != 0:
        millibad = (1000 * (row['mismatches'] + insertfactor + round(3*math.log(1 + sizediff)))) / total
    return 100.0 - millibad * 0.1
   

if __name__ == '__main__':
    import sys
    if len(sys.argv) <= 1:
        print("Path to .psl file is missing! ")
        print("Usage {} ~/analysis/14522_B6a/analysis/output/rrna_fragments/Ctrl_fraction_3.psl".format(sys.argv[0]))
        exit(1)
    psl_file = sys.argv[1]
    import pandas as pd
    print("Reading file: {}".format(psl_file))
    df = pd.read_csv(psl_file, sep="\t", names=["matches", "mismatches", "repmatches", "ncount", "qnuminsert", "qbaseinsert", "tnuminsert", "tbaseinsert", "strand", "qname", "qsize", "qstart", "qend", "tname", "tsize", "tstart", "tend", "blockcount", "blocksizes", "qstarts", "tstarts"], header=None)
    import pdb; pdb.set_trace()
    df['score'] = df['matches'] + (df['repmatches'] / 2) - df['mismatches'] - df['qnuminsert'] - df['tnuminsert']
    df['tspan'] = df['tend'] - df['tstart']
    df['qspan'] = df['qend'] - df['qstart']
    df['perc_ident'] = df.apply(calcPercentIdentity)
    df = df[['qname', 'score', 'qstart', 'quend', 'qsize', 'perc_ident', 'tname', 'strand', 'tstart', 'tend', 'tspan', 'qblockseqs', 'tblockseqs']]
    import pdb; pdb.set_trace()
    df.to_csv(psl_file.repalce('.psl', '_parsed.psl', sep="\t", header=False, index=False))
    # Replicate BLAT output from web
#    for line in open(psl_file):
#        p = Pslx(line)
#        print(p.qname, p.score(), p.qstart+1, p.qend, p.qsize,
#            "%.1f" % p.calcPercentIdentity(), p.tname, p.strand,
#            p.tstart+1, p.tend, p.tspan(), p.qblockseqs, p.tblockseqs)
    
