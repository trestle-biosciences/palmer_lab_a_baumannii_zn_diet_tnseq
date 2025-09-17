#!/usr/bin/env python
"""Convert BAM files to Transit Ready Wig file format
adapted from file located here:
Usage:
    bam_to_wig.py <BAM file> [<YAML config>]
    [--outfile=<output file name>
     --chrom=<chrom>
     --start=<start>
     --end=<end>
     --fasta=<fasta>
     --tn5=<tn5>
     --tempfile=<tempfile>]
     --ta-position=<ta-position>
chrom start and end are optional, in which case they default to everything.
chrom is optional. If not provided all chromosomes will be converted.
"""
import os
import sys
import tempfile
from optparse import OptionParser
from contextlib import contextmanager, closing
import pysam

def main(bam_file, chrom='all', start=0, end=None,
         outfile=None, normalize=False, use_tempfile=False,
         fasta_file=None, transposon_type='m', ta_position=None):

    if outfile is None:
        outfile = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if start > 0:
        start = int(start) - 1
    if end is not None:
        end = int(end)
    regions = [(chrom, start, end)]
    if os.path.abspath(bam_file) == os.path.abspath(outfile):
        sys.stderr.write("Bad arguments, input and output files are the same.\n")
        sys.exit(1)
    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
        if use_tempfile:
            #Use a temp file to avoid any possiblity of not having write permission
            out_handle = tempfile.NamedTemporaryFile(delete=False)
            wig_file = out_handle.name
        else:
            wig_file = "%s.wig" % os.path.splitext(outfile)[0]
            out_handle = open(wig_file, "w")
        with closing(out_handle):
            chr_sizes, wig_valid = write_bam_track(bam_file, regions, out_handle,
                                                   normalize, fasta_file, transposon_type, 
                                                   ta_position)

def get_dict_value(dict, key, default):
    try:
        return dict[key]
    except:
        return default

@contextmanager
def indexed_bam(bam_file):
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.AlignmentFile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()

def write_bam_track(bam_file, regions, out_handle, normalize, fasta_file, transposon_type, ta_position):
    out_handle.write("track %s\n" % " ".join(["type=wiggle_0",
        "name=%s" % os.path.splitext(os.path.split(bam_file)[-1])[0],
        "visibility=full",
        ]))
    normal_scale = 1e6
    is_valid = False
    fasta = pysam.FastaFile(fasta_file)
    with indexed_bam(bam_file) as work_bam:
        total = sum(1 for r in work_bam.fetch() if not r.is_unmapped) if normalize else None
        sizes = zip(work_bam.references, work_bam.lengths)
        if len(regions) == 1 and regions[0][0] == "all":
            regions = [(name, 0, length) for name, length in sizes]
        counts = {}
        for chrom, start, end in regions:
            if chrom in fasta.references:
                if end is None and chrom in work_bam.references:
                    end = work_bam.lengths[work_bam.references.index(chrom)]
                assert end is not None, "Could not find %s in header" % chrom
                out_handle.write("variableStep chrom=%s\n" % chrom)
                # get counts at each position
                counts[chrom] = {}  
                counts[chrom]['-'] = {}
                counts[chrom]['+'] = {}            
                if ta_position == 'internal_3prime':         
                    for read in work_bam.fetch(chrom, start, end):
                        if read.is_reverse == 1:
                            try:
                                counts[chrom]['-'][read.reference_start] +=1
                            except:
                                counts[chrom]['-'][read.reference_start] = 1
                        else:
                            try:
                                counts[chrom]['+'][read.reference_end-2] +=1  # minus 2:  -1 for python coords, -1 to align with positive strand T counts
                            except:
                                counts[chrom]['+'][read.reference_end-2] = 1
                elif ta_position == 'adjacent_5prime':
                    for read in work_bam.fetch(chrom, start, end):
                        if read.is_reverse == 1:
                            try:
                                counts[chrom]['-'][read.reference_end] +=1
                            except:
                                counts[chrom]['-'][read.reference_end] = 1
                        else:
                            try:
                                counts[chrom]['+'][read.reference_start-2] +=1  # minus 2:  to align with positive strand T site
                            except:
                                counts[chrom]['+'][read.reference_start-2] = 1


                
                print("gathered counts")

                # if transposon type is not tn5, only use TA site counts
                # and ensure ta sites with no counts are represented as zeros (per Transit model requirements)
                if transposon_type.lower() != 'tn5':
                    chrom_seq = fasta[chrom]
                    ta = {}
                    ta['+'] = [i for i in range(len(chrom_seq)) if chrom_seq[i:i+2] == 'TA']
                    # ta['-'] = [i+1 for i in range(len(chrom_seq)) if chrom_seq.lower()[i:i+2] == 'ta']

                    count_at_ta = sum([count for strand in counts[chrom].keys()
                                            for pos, count in counts[chrom][strand].items()
                                            if pos in ta['+']
                                            ])
                    count_not_at_ta = sum([count for strand in counts[chrom].keys()
                                            for pos, count in counts[chrom][strand].items()
                                            if pos not in ta['+']
                                            ])
                    total = count_at_ta
                    print(f"Total TA site count for region ({chrom}, {start}, {end}) in file {bam_file} is {count_at_ta}")
                    print(f"Total non-TA site count is {count_not_at_ta}")
                    print(f"The ratio of counts at a TA site is {count_at_ta/(count_not_at_ta + count_at_ta)}")

                    # ta_sites = {pos: 0 for pos in ta['+']}
                    ta_sites = {pos:sum(get_dict_value(counts[chrom][strand], pos, 0) 
                                            for strand in counts[chrom].keys())
                                        for pos in ta['+']}

                    count_of_ta_sites = len(ta_sites.keys())
                    count_of_ta_sites_with_counts = sum([1 for i,v in ta_sites.items() if v > 0])
                    percent_of_ta_sites_with_counts = count_of_ta_sites_with_counts/count_of_ta_sites

                    print(f"Total count of TA sites in region ({chrom}, {start}, {end}) in file {bam_file} is {count_of_ta_sites}")
                    print(f"Of these,  {count_of_ta_sites_with_counts} are have non-zero counts, giving an occupancy ratio of {percent_of_ta_sites_with_counts}")
    
                    for pos, count in ta_sites.items():
                        if normalize:
                            n = float(count) / total * normal_scale
                        else:
                            n = count
                        out_handle.write("%s %.1f\n" % (pos+1, n))
                    is_valid = True
                else:
                    print("NOT YET SUPPORTED")
        return sizes, is_valid



if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-o", "--outfile", dest="outfile")
    parser.add_option("-c", "--chrom", dest="chrom")
    parser.add_option("-s", "--start", dest="start")
    parser.add_option("-e", "--end", dest="end")
    parser.add_option("-n", "--normalize", dest="normalize",
                      action="store_true", default=False)
    parser.add_option("-t", "--tempfile", dest="use_tempfile",
                      action="store_true", default=False)   
    parser.add_option("-f", "--fasta", dest="fasta_file")
    parser.add_option("-y", "--transposon", dest="transposon_type", default='m')
    parser.add_option("-p", "--ta-position", dest="ta_position")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        print("Incorrect arguments")
        print(__doc__)
        sys.exit()
    kwargs = dict(
        outfile=options.outfile,
        chrom=options.chrom or 'all',
        start=options.start or 0,
        end=options.end,
        use_tempfile=options.use_tempfile,
        fasta_file=options.fasta_file,
        transposon_type=options.transposon_type,
        ta_position = options.ta_position)
    main(*args, **kwargs)
