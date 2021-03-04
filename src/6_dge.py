#!/usr/bin/python

from __future__ import division
from numpy import cumsum
from collections import Counter, defaultdict
from itertools import chain, combinations, product
import pandas as pd
import datetime, os, operator, gzip, warnings, timeit

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib.pyplot as plt

__author__  = 'Shan Sabri'
__date__    = '2016/12/09'

wrdir  = '/path/to/working/directory/'

extend_knee = 5000
take_top_n_cells = 2500
barnyard_threshold = 1 # 70
barnyard_read_cutoff = 1 # 2500

input_format = 'uniq.stranded.bed.gz'

CHARS     = 'ACGT'
bc_flag   = 'BC:'
gene_flag = 'GENE:'
umi_flag  = 'UMI:'

VERBOSE = False
with_mixed = False
with_uncollapsed_dges = False

start = datetime.datetime.now()


def init_dir(files):
    print '{}\tCreating necessary directories'.format(datetime.datetime.now()-start)
    if not os.path.exists(os.path.join(wrdir, "output")):
        os.makedirs(os.path.join(wrdir, "output"))
    outdir = os.path.join(wrdir, "output")
    for f in files:
        if not os.path.exists(os.path.join(outdir, f)):
            os.makedirs(os.path.join(outdir, f))
    return outdir


def file_len(file):
    with gzip.open(file) as f:
        for i, l in enumerate(f, start=1):
            pass
    return i


def extract_barcodes(file):
    with gzip.open(file) as f:
        for idx, l in enumerate(f):
            l = l.strip().split('\t')
            obs_bc = [x for x in l if bc_flag in x][0].split(':')[1]
            yield obs_bc


def count_cell_barcodes(file):
    print datetime.datetime.now()-start, '\t- Creating barcode counter dictionary'
    barcode_counts = Counter(extract_barcodes(file))
    print('{}\t-- Identified {} barcodes corresponding to {} reads'.format(datetime.datetime.now()-start,
                                                                           len(barcode_counts), sum(barcode_counts.values())))
    return barcode_counts


def count_umis(file, barcodes):
    print '{}\t- Creating cell:[gene:transcipt] dictionary'.format(datetime.datetime.now() - start)
    cells = {}
    for idx, l in enumerate(gzip.open(file)):
        l = l.strip().split("\t")
        gene    = [x for x in l if gene_flag in x]
        obs_umi = [x for x in l if umi_flag in x]
        obs_bc  = [x for x in l if bc_flag in x]
        gene    = gene[0].split(":")[1].upper()
        obs_bc  = obs_bc[0].split(":")[1]
        obs_umi = obs_umi[0].split(":")[1]
        if obs_bc not in barcodes: continue
        if obs_bc not in cells: cells[obs_bc] = {}
        if gene not in cells[obs_bc]: cells[obs_bc][gene] = Counter()
        cells[obs_bc][gene].update([obs_umi])
    return cells


def hamming_circle(s, n):
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(CHARS) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == CHARS[r]:
                    cousin[p] = CHARS[-1]
                else:
                    cousin[p] = CHARS[r]
            yield ''.join(cousin)


def hamming_ball(s, n):
    return chain.from_iterable(hamming_circle(s, i) for i in range(n + 1))


def find_mergable(kmers, i, edits):
    if i in kmers:
        return i
    for s in hamming_ball(i, edits):
        if s in kmers:
            return s


def find_cell_barcodes_to_collapse(x, edits=1, cutoff=5):
    print '{}\t- Identifying cell barcodes to collapse'.format(datetime.datetime.now() - start)
    seen = set()
    collapsed = {}
    for k, v in x.most_common():
        if v < cutoff: continue
        true = find_mergable(seen, k, edits)
        if true:
            if VERBOSE:
                print '{}\t-- Collapsing {} ({} occurrences) with {} ({} occurrences)'.format(datetime.datetime.now() - start,
                                                                                              k, v, true, x[true])
            collapsed[k] = true
        else:
            seen.add(k)
    return collapsed


def collapse_cell_barcodes(f, dict, out):
    print '{}\t- Collapsing barcodes'.format(datetime.datetime.now() - start)
    with gzip.open(f) as i, gzip.open(out, 'w') as o:
        for idx, l in enumerate(i, start=1):
            l = l.strip().split('\t')
            obs  = [x for x in l if bc_flag in x]
            obs  = obs[0].split(':')[1]
            true = dict.get(obs)
            if true:
                if VERBOSE:
                    print '{}\t-- Replaceing less common {} with more common {} at line {}'.format(
                        datetime.datetime.now() - start, obs, true, idx)
                l[l.index(bc_flag + obs)] = bc_flag + true
            o.write('\t'.join(l) + '\n')


def collapse_umis(x, edits=1):
    print '{}\t- Figuring which UMIs to collapse'.format(datetime.datetime.now() - start)
    collapsed = {}
    for cell, gene_dict in x.iteritems():
        collapsed[cell] = {}
        for gene, umis in gene_dict.iteritems():
            collapsed[cell][gene] = len(umis)
            if len(umis) == 1: continue
            seen = set()
            for k, v in umis.most_common():
                true = find_mergable(seen, k, edits)
                if true:
                    if VERBOSE:
                        print '{}\t\t\tCan collapse {} ({} occurrences) with {} ({} occurrences) -- {}[{}]'.format(datetime.datetime.now()-start,
                                                                                                                   k, umis[k], true, umis[true], cell, gene)
                    umis[true] = umis[k] + umis[true]
                    del umis[k]
                else:
                    seen.add(k)
            collapsed[cell][gene] = len(umis)
    return collapsed


def plot_knee(x, outfile, xlim=(0,5000), h_line=0, close=True):
    print '{}\t- Plotting knee'.format(datetime.datetime.now()-start)
    total = sum(x.values())
    cum_perc = map(lambda x: round(x*100.0/total, 10), cumsum(sorted(x.values(), reverse=True)))

    if VERBOSE:
        for i in range(100, 1100, 100):
            print '{}\t-- {} cells capture {}% of total reads'.format(datetime.datetime.now()-start, i, cum_perc[i])

    # most_common = x.most_common()
    # with open(outfile, 'w') as o:
    #     o.write('Barcode\tReads\tCumulativePercent\n')
    #     for percent, kv in zip(cumPerc, most_common[:n]):
    #         key, value = kv
    #         o.write('%s\t%s\t%s\n' % (key, value, percent))

    plt.plot(cum_perc)
    plt.title(str(h_line) + " cells captured\n" + ' total.reads:' + str(total) + ', total.bcs:' + str(len(cum_perc)))
    plt.xlim(xlim)
    plt.axvline(x=h_line, ymin=0, ymax=1, color='gray', linestyle='--')
    plt.ylabel('Cumulative Fraction of Reads')
    plt.xlabel('Cellular Barcodes Sorted by Number of Reads [decending]')
    plt.savefig(outfile)
    if close: plt.close()


def get_top_cell_barcodes(x, n=1000):
    top = dict(x.most_common(n))
    return top.keys()


def split_by_species(f, species, out):
    print '{}\t- Splitting file by {} cells'.format(datetime.datetime.now()-start, species)
    with gzip.open(f) as i, gzip.open(out, 'w') as o:
        for l in i:
            l = l.strip().split('\t')
            chr = l[1]
            if species in chr:
                o.write('\t'.join(l)+'\n')


def split_by_mixed_cells(f, mixed_cells, out):
    print '{}\t- Splitting file by MIXED cells'.format(datetime.datetime.now() - start)
    with gzip.open(f) as i, gzip.open(out, 'w') as o:
        for idx, l in enumerate(i):
            l = l.strip().split('\t')
            bc = [x for x in l if bc_flag in x]
            bc = bc[0].split(':')[1]
            if bc in mixed_cells:
                o.write('\t'.join(l) + '\n')


def get_barnyard_info(file, top_barcodes):
    barnyard = {}
    with gzip.open(file) as f:
        for l in f:
            l = l.strip().split('\t')
            bc = [x for x in l if bc_flag in x]
            bc = bc[0].split(':')[1]
            if bc in top_barcodes:
                if bc in barnyard:
                    barnyard[bc] += 1
                else:
                    barnyard[bc] = 1
            else:
                barnyard[bc] = 0
    return barnyard


def identity_x_mapped_reads(f1, f2):
    print '{}\t- Identifing species X-mapped reads'.format(datetime.datetime.now() - start)
    reads, x_mapped_reads = set(), set()
    with gzip.open(f1) as m:
        for ml in m:
            ml = ml.strip().split('\t')
            reads.add(ml[0])
    with gzip.open(f2) as h:
        for hl in h:
            hl = hl.strip().split('\t')
            if hl[0] in reads:
                x_mapped_reads.add(hl[0])
    print '{}\t-- Found {} X-mapped reads'.format(datetime.datetime.now() - start, len(x_mapped_reads))
    return x_mapped_reads


def remove_x_mapped_reads(file, x_mapped_reads, out):
    with gzip.open(file) as f, gzip.open(out, 'w') as o:
        for l in f:
            l = l.strip().split('\t')
            if l[0] in x_mapped_reads: continue
            o.write('\t'.join(l) + '\n')


def plot_barnyard(barnyard, top_barcodes, outfile, rcfile=None, threshold=75, exclude_cells_below=10000):
    data = {'x':[], 'y':[], 'label':[]}

    if rcfile:
        with gzip.open(rcfile, 'w') as rc:
            rc.write('Barcode\tMouse\tHuman\n')
            for label, coord in barnyard.iteritems():
                if label not in top_barcodes: continue
                if len(coord) != 2: continue
                if coord[0] < exclude_cells_below and coord[1] < exclude_cells_below: continue
                data['x'].append(coord[0])
                data['y'].append(coord[1])
                data['label'].append(label)
                rc.write(str.join('\t', [label] + [str(c) for c in coord]) + '\n')
    else:
        for label, coord in barnyard.iteritems():
            if label not in top_barcodes: continue
            if len(coord) != 2: continue
            if coord[0] < exclude_cells_below and coord[1] < exclude_cells_below: continue
            data['x'].append(coord[0])
            data['y'].append(coord[1])
            data['label'].append(label)

    plt.figure(figsize=(10,8))
    plt.title('Barnyard Plot', fontsize=20)
    plt.xlabel('MOUSE Reads', fontsize=15)
    plt.ylabel('HUMAN Reads', fontsize=15)
    plt.xlim(0, max(max(data['x']), max(data['y'])))
    plt.ylim(0, max(max(data['x']), max(data['y'])))
    plt.axvline(x=exclude_cells_below, ymin=0, ymax=exclude_cells_below/max(max(data['x']), max(data['y'])),
                color='black', linestyle='-')
    plt.axhline(y=exclude_cells_below, xmin=0, xmax=exclude_cells_below/max(max(data['x']), max(data['y'])),
                color='black', linestyle='-')

    mouse_barcodes, human_barcodes, mixed_barcodes = [], [], []
    mouse_x, human_x, mixed_x = [], [], []
    mouse_y, human_y, mixed_y = [], [], []

    for label, x, y  in zip(data['label'], data['x'], data['y']):
        x, y = int(x), int(y)
        percent_mouse = int(round((x/(x+y))*100))
        percent_human = int(round((y/(x+y))*100))
        if percent_mouse >= threshold:
            mouse_barcodes.append(label)
            mouse_x.append(x)
            mouse_y.append(y)
        elif percent_human >= threshold:
            human_barcodes.append(label)
            human_x.append(x)
            human_y.append(y)
        else:
            mixed_barcodes.append(label)
            mixed_x.append(x)
            mixed_y.append(y)

    m = plt.scatter(x=mouse_x, y=mouse_y, label=mouse_barcodes, color="blue")
    h = plt.scatter(x=human_x, y=human_y, label=human_barcodes, color="red")
    x = plt.scatter(x=mixed_x, y=mixed_y, label=mixed_barcodes, color="gray")
    plt.legend((m, h, x), ('Mouse: ' + str(len(mouse_barcodes)),
                           'Human: ' + str(len(human_barcodes)),
                           'Mixed: ' + str(len(mixed_barcodes))),
               scatterpoints=1, loc=1, fontsize=8)
    plt.savefig(outfile)
    plt.close()

    return mouse_barcodes, human_barcodes, mixed_barcodes


def write_dge(dge, out):
    print "{}\t- Writing DGE".format(datetime.datetime.now()-start)
    df = pd.DataFrame(dge).fillna(0)
    df.to_csv(out, sep="\t", compression='gzip')


def write_lib_info(x, outfile):
    print "{}\t- Writing summary".format(datetime.datetime.now() - start)
    with open(outfile, "wb") as o:
        o.write('\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ("Total.reads", "Total.reads.collapsed",
                                                        "Total.reads.collapsed.mouse", "Total.reads.collapsed.human",
                                                        "Total.uncollapsed.UMIs", "Total.collapsed.UMIs",
                                                        "Total.reads.collapsed.mouse.cleaned", "Total.reads.collapsed.human.cleaned"))
        for k, v in x.iteritems():
            o.write('%s\t%s\n' % (k, "\t".join(str(i) for i in v)))

def clean_up(files):
    print "{}\t- Cleaning up".format(datetime.datetime.now() - start)
    for f in files:
        os.remove(f)


def pipeline():
    files = [os.path.splitext(f)[0].split(".")[0] for f in os.listdir(wrdir) if f.endswith(input_format)]

    outdir = init_dir(files)

    for f in files:
        print '{}\tProcessing {}'.format(datetime.datetime.now()-start, f)

        # lib_info = defaultdict(list)

        # STEP 0: Count the number of uncollapsed transcripts (# of lines in file)
        total_uncollapsed_transcripts = file_len(os.path.join(wrdir, '.'.join((f, input_format))))

        # STEP 1: Count the raw occurences of each barcode
        barcode_counts = count_cell_barcodes(os.path.join(wrdir, '.'.join((f, input_format))))
        # for k, v in barcode_counts.iteritems(): lib_info[k].append(v) # append num reads to lib info

        # STEP 2: Plot knee
        plot_knee(barcode_counts, outfile=os.path.join(outdir, f, '.'.join((f, "knee", "uncollapsed", "pdf"))),
                  xlim=(0, extend_knee), h_line=take_top_n_cells ,close=False)

        # STEP 3: Identify which cell barcodes to collapse
        barcodes_to_collapse = find_cell_barcodes_to_collapse(barcode_counts, edits=2, cutoff=50)

        # STEP 4: Collapse on the identified barcodes
        collapse_cell_barcodes(f=os.path.join(wrdir, '.'.join((f, input_format))),
                               dict=barcodes_to_collapse,
                               out=os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))))

        # STEP 5: sort | uniq on the collapsed 20-mer
        uniq_20mers = set()
        with gzip.open(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format)))) as collapsed_bed:
                for idx, l in enumerate(collapsed_bed):
                    l = l.strip().split('\t')
                    uniq_20mers.add(l[11])
        with open(os.path.join(outdir, f, '.'.join((f, "uniq.fraction.txt"))), 'w') as o:
            o.write('Uncollapsed_Transcripts\tCollapsed_Transcripts\tCollapse_Rate\n')
            o.write('{}\t{}\t{}\n'.format(str(total_uncollapsed_transcripts), str(len(uniq_20mers)),
                                          str(len(uniq_20mers)/total_uncollapsed_transcripts)))

        collapsed_barcode_counts = count_cell_barcodes(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))))
        # for k, v in collapsed_barcode_counts.iteritems(): lib_info[k].append(v)

        # STEP 5: Plot collapsed knee
        plot_knee(collapsed_barcode_counts, outfile=os.path.join(outdir, f, '.'.join((f, "knee", "collapsed", "pdf"))),
                  xlim=(0,extend_knee), h_line=take_top_n_cells, close=True)

        # STEP 6: Get most abundant barcodes from collapsed cell barcode counts dictionary
        top_collapsed_cell_barcodes = get_top_cell_barcodes(x=collapsed_barcode_counts, n=take_top_n_cells)


        # STEP 7: Split collapsed bed file by species and count number of reads within the top collapsed cell barcodes
        species = ["MOUSE", "HUMAN"]
        master_barnyard = {}
        for s in species:
            split_by_species(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))), species=s,
                             out=os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))))
            species_barcode_counts = count_cell_barcodes(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))))
            # for k, v in species_barcode_counts.iteritems(): lib_info[k].append(v)
            species_barnard_info = get_barnyard_info(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))), top_barcodes=top_collapsed_cell_barcodes)
            for k, v in sorted(species_barnard_info.iteritems(), key=operator.itemgetter(1), reverse=True):
                master_barnyard.setdefault(k, []).append(v)

        # STEP 8: Plot the barnyard
        plot_barnyard(master_barnyard,
                      top_barcodes=top_collapsed_cell_barcodes,
                      outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "pdf"))),
                      threshold=barnyard_threshold,
                      exclude_cells_below=0)

        mouse, human, mixed = plot_barnyard(master_barnyard,
                                            top_barcodes=top_collapsed_cell_barcodes,
                                            outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "cutoff", "pdf"))),
                                            rcfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "reads.counts", "txt.gz"))),
                                            threshold=barnyard_threshold,
                                            exclude_cells_below=barnyard_read_cutoff)

        # STEP 9: Write out bed files with MIXED cells
        if with_mixed:
            split_by_mixed_cells(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))), mixed_cells=mixed,
                                 out=os.path.join(outdir, f, '.'.join((f, "collapsed", "MIXED", input_format))))


        # STEP 9: Collapse UMIs and write out DGEs in a species-specific manner
        species = {'mouse': mouse, 'human': human}
        if with_mixed:
            species['mixed'] = mixed

        for s, v in species.iteritems():
            print '{}\t- Computing DGE for {} - {}'.format(datetime.datetime.now() - start, f, s)
            umi_counts = count_umis(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))), v)

            umi_counts_uncollapsed = {}
            for cell, gene_dict in umi_counts.iteritems():
                umi_counts_uncollapsed[cell] = {}
                for gene, umis in gene_dict.iteritems():
                    umi_counts_uncollapsed[cell][gene] = sum(umis.values())
            # for k, v in umi_counts_uncollapsed.iteritems():  lib_info[k].append(sum(v.values()))

            umi_counts_collapsed = collapse_umis(umi_counts, edits=1)
            # for k, v in umi_counts_collapsed.iteritems():  lib_info[k].append(sum(v.values()))

            write_dge(umi_counts_collapsed, os.path.join(outdir, f, '.'.join((f, "collapsed", s, "dge", "tsv.gz"))))

            if with_uncollapsed_dges:
                write_dge(umi_counts, os.path.join(outdir, f, '.'.join((f, "uncollapsed", s, "dge.showUMIs", "tsv.gz"))))
                write_dge(umi_counts_uncollapsed, os.path.join(outdir, f, '.'.join((f, "uncollapsed", s, "dge", "tsv.gz"))))


        # # STEP 10: Remove cross-mapped reads from species collapsed files
        x_mapped_reads = identity_x_mapped_reads(f1=os.path.join(outdir, f, '.'.join((f, "collapsed.MOUSE", input_format))),
                                                 f2=os.path.join(outdir, f, '.'.join((f, "collapsed.HUMAN", input_format))))
        cleaned_master_barnyard = {}
        for s in ["MOUSE", "HUMAN"]:
            remove_x_mapped_reads(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))),
                                  x_mapped_reads,
                                  out=os.path.join(outdir, f, '.'.join((f, "collapsed", s, "cleaned", input_format))))
            cleaned_species_barcode_counts = count_cell_barcodes(os.path.join(outdir, f, '.'.join((f, "collapsed", s, "cleaned", input_format))))
            # for k, v in cleaned_species_barcode_counts.iteritems(): lib_info[k].append(v)

            cleaned_species_barnard_info = get_barnyard_info(os.path.join(outdir, f, '.'.join((f, "collapsed", s, "cleaned", input_format))),
                                                             top_barcodes=top_collapsed_cell_barcodes)
            for k, v in sorted(cleaned_species_barnard_info.iteritems(), key=operator.itemgetter(1), reverse=True):
                cleaned_master_barnyard.setdefault(k, []).append(v)

        # STEP 11: Plot the barnyard
        plot_barnyard(cleaned_master_barnyard, top_barcodes=top_collapsed_cell_barcodes,
                      outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "cleaned", "pdf"))),
                      threshold=barnyard_threshold, exclude_cells_below=0)

        plot_barnyard(cleaned_master_barnyard, top_barcodes=top_collapsed_cell_barcodes,
                      outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "cleaned", "cutoff", "pdf"))),
                      threshold=barnyard_threshold, exclude_cells_below=barnyard_read_cutoff)


        # STEP 13: Write out library info
        # write_lib_info(lib_info, os.path.join(outdir, f, '.'.join((f, "summary", "tsv")))) # not sure how to have 0 place holders for rows with no entry?

        # STEP 14: garbage collect
        # delete_these_files = [os.path.join(outdir, f, '.'.join((f, "knee.uncollapsed.pdf"))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed", input_format)))]
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.MOUSE", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.HUMAN", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.MIXED", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.MOUSE.cleaned", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.HUMAN.cleaned", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "uncollapsed.mouse.dge.tsv.gz"))),
        #                       os.path.join(outdir, f, '.'.join((f, "uncollapsed.human.dge.tsv.gz"))),
        #                       os.path.join(outdir, f, '.'.join((f, "uncollapsed.mixed.dge.tsv.gz")))]
        # clean_up(delete_these_files)


if __name__ == '__main__':

    print 'START: ', datetime.datetime.now()
    pipeline()
    print 'FINISH: ', datetime.datetime.now()
