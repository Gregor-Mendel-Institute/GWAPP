#!/usr/local/bin/python2.7
# encoding: utf-8
'''
plotter -- shortdesc

plotter is a wrapper script for plotting GWAS results

It defines classes_and_methods

@author:     Ümit Seren
        
@copyright:  2012 Gregor Mendel Institute. All rights reserved.
        
@license:    MIT

@contact:    uemit.seren@gmail.com
@deffield    updated: Updated
'''

import sys
import os

import numpy
import tables
import itertools
import math
import gwaResults
__all__ = []
__version__ = 0.1
__date__ = '2012-12-06'
__updated__ = '2012-12-06'




def main(argv=None):
    '''Command line options.'''
    from argparse import ArgumentParser
    from argparse import RawDescriptionHelpFormatter

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by Ümit Seren on %s.
  Copyright 2012 Gregor Mendel Institute. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-r", "--result", dest="result", help="result for which to plot [PHENOTYPE1/Fullset/raw/lm]",required=True)
        parser.add_argument("-o", "--output_folder", dest="outputfolder", help="Location where to output the plots")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('--png',dest='png',action='store_true',default=False)
        parser.add_argument('--pdf',dest='pdf',action='store_true',default=False)
        parser.add_argument(dest="file", help="HDF5 GWAS result file", metavar="FILE")
        # Process arguments
        args = parser.parse_args()
        _plot_gwas_result(args.file,args.outputfolder,args.result,args.png,args.pdf)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def _plot_gwas_result(hdf5_file,outputfolder,result,png=True,pdf=False,mode='r',mac=15):
    #ipdb.set_trace()
    data,threshold = _get_gwas_result(hdf5_file,result,mode)
    pdf_file = None
    png_file = None
    base_folder = './'
    filename = '%s_mac%s' % (result.replace('/','_'),int(mac))
    if outputfolder is not None:
        base_folder=outputfolder
    image_file = ''
    if png:
        png_file = '%s/%s.png' % (base_folder,filename)
        image_file = png_file
    if pdf:
        pdf_file= '%s/%s.pdf' % (base_folder,filename)
        image_file= pdf_file

    res = gwaResults.Result(chromosomes=data['chromosome'],positions=data['position'],scores=data['score'],mafs=data['maf'],macs=data['mac'])
    res.filter_attr('macs',mac)
    res.plot_manhattan(png_file=png_file,pdf_file=pdf_file, percentile=90, type="score",
                            ylab="$-$log$_{10}(p)$", plot_bonferroni=True,b_threshold=-math.log10(threshold),
                            neg_log_transform=False)
    return image_file
    #_plot_manhattan_chart(data,threshold=threshold,png_file=png_file,pdf_file=pdf_file)

def _get_chromosome_ends(data,chr_list):
    sorted_ix = numpy.argsort(data,order=('chromosome','position'))
    sorted_by_chr_pos = data[sorted_ix]
    #data.sort(order=['chromosome','position'])
    chromosome_ends = []
    for chr in chr_list:
        if chr == 1:
            continue
        chromosome_ends.append(sorted_by_chr_pos['position'][numpy.where(sorted_by_chr_pos['chromosome'] == chr)[0][0]-1])

    chromosome_ends.append(sorted_by_chr_pos['position'][len(sorted_by_chr_pos['chromosome'])-1])
    return chromosome_ends

def _get_chromosome_splits(data,chr_list):
    """
        Returns list of indices (and prev chromosome), for the when the chromosomes
        change in the scores, positions indices.
        """
    chromosome_splits = []
    for chr in chr_list:
        chromosome_splits.append(numpy.where(data['chromosome'] == chr)[0][0])
    chromosome_splits.append(len(data['chromosome'])-1)
    return chromosome_splits

def _plot_manhattan_chart(data,b_threshold=None, percentile=100,highlight_markers = None,
                          highlight_loci=None,plot_xaxis=True,chrom_col_map=None,
                          markersize= 6,max_score=None,min_score=None,plot_bonferroni=False,threshold=None,pdf_file=None,png_file=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    num_scores = len(data['score'])

    "Plotting a Manhattan-style plot with %i markers." % num_scores
    chrom_set = set(data['chromosome'])
    chromosomes = list(chrom_set)
    chromosomes.sort()
    chromosome_ends = _get_chromosome_ends(data,chromosomes)

    if len(chrom_set) == 1:
        percentile = 0.0

    if percentile != 0.0:
        data = data[:num_scores*percentile]

    if highlight_markers:
        new_h_markers = []
        if len(highlight_markers[0]) == 2:
            indices = result.get_indices(highlight_markers)
            scores = [result.snp_results['scores'][i] for i in indices]
            for i, (c, p) in enumerate(highlight_markers):
                new_h_markers.append((c, p, scores[i]))
        else:
            new_h_markers = [(cps[0], cps[1], cps[2]) for cps in highlight_markers]
        highlight_markers = new_h_markers

    if not max_score:
        max_score = max(data['score'])
        if highlight_markers:
            h_scores = [s for c, p, s in highlight_markers]
            max_score = max(max_score, max(h_scores))
    if not min_score:
        min_score = min(data['score'])


    if highlight_loci:
        hl_dict = {}
        for chrom in chromosomes:
            hl_dict[chrom] = []
        for c, p in highlight_loci:
            hl_dict[c].append(p)



    scoreRange = max_score - min_score
    offset = 0
    chromosome_splits = _get_chromosome_splits(data,chromosomes)

    ticksList1 = []
    ticksList2 = []
    plt.figure(figsize=(11, 2.8))
    plt.axes([0.045, 0.15, 0.95, 0.71])
    starPoints = [[], [], []]
    chr_offsets = []
    for i, chromosome_end in enumerate(chromosome_ends):
        chr_offsets.append(offset)
        print i, chromosome_splits
        index1 = chromosome_splits[i]
        index2 = chromosome_splits[i + 1]
        scoreList = data['score'][index1:index2]
        posList = data['position'][index1:index2]
        chrom = chromosomes[i]
        newPosList = [offset + pos for pos in posList]

        for s_i, (score, pos) in enumerate(itertools.izip(scoreList, newPosList)):
            if score > max_score:
                starPoints[0].append(pos)
                starPoints[1].append(max_score)
                starPoints[2].append(score)
                score = max_score
            scoreList[s_i] = score

        if not chrom_col_map:
            plt.plot(newPosList, scoreList, ".", markersize=markersize, alpha=0.7, mew=0)
        else:
            color = chrom_col_map[chrom]
            plt.plot(newPosList, scoreList, ".", markersize=markersize, alpha=0.7, color=color, mew=0)

        if highlight_loci:
            for p in hl_dict[chrom]:
                plt.axvline(offset + p, color='#1166FF', alpha=0.6)

        oldOffset = offset
        #            textPos.append(offset + chromosome_end / 2 - 2000000)
        offset += chromosome_end
        if plot_xaxis:
            if len(chromosome_ends) == 23: #This is probably Human!
                ticksList1.append(oldOffset + chromosome_end / 2)
                if chrom == 23:
                    ticksList2.append('X')
                else:
                    ticksList2.append(chrom)
            else:
                for j in range(oldOffset, offset, 4000000):
                    ticksList1.append(j)
                for j in range(0, chromosome_end, 4000000):
                    if j % 8000000 == 0 and j < chromosome_end - 4000000 :
                        ticksList2.append(j / 1000000)
                    else:
                        ticksList2.append("")


    plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=markersize + 2, mew=0)
    if len(starPoints[0]) > 0:
        i = 0
        while i < len(starPoints[0]):
            max_point = i
            cur_pos = starPoints[0][i]
            while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 3000000:
                if starPoints[2][i] > starPoints[2][max_point]:
                    max_point = i
                i += 1
            plt.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")


    if highlight_markers:
        ys = []
        xs = []
        for c, p, score in highlight_markers:
            x = chr_offsets[c - 1] + p
            xs.append(x)
            if score > max_score:
                plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
                ys.append(max_score)
            else:
                ys.append(score)
        plt.plot(xs, ys, ".", color="#ff9944", markersize=markersize + 4, alpha=0.8, mew=0)


    if plot_bonferroni:
        if not b_threshold:
            b_threshold = -math.log10(1.0 / (num_scores * 20.0))
        plt.plot([0, sum(chromosome_ends)], [b_threshold, b_threshold], color='#000000', linestyle="-.")
    if threshold :
        threshold = -math.log10(threshold)
        plt.plot([0, sum(chromosome_ends)], [threshold, threshold], color='#6495ed', linestyle='-.')


    x_range = sum(chromosome_ends)
    plt.axis([-x_range * 0.01, x_range * 1.01, min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
    if plot_xaxis:
        plt.xticks(ticksList1, ticksList2, fontsize='x-small')
        #pdb.set_trace()
    plt.ylabel('$-$log$_{10}(p)$')

    if plot_xaxis:
        if len(chromosome_ends) == 23: #This is probably Human!
            plt.xlabel("Chromosome")
        else:
            plt.xlabel("Mb")
    else:
        plt.xlabel("bases")

    image_file = ''
    if pdf_file:
        plt.savefig(pdf_file, format="pdf")
        image_file = pdf_file
    if png_file:
        plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
        image_file = png_file
    if not (pdf_file or png_file):
        plt.show()
    plt.clf()
    plt.close()
    return image_file

def _get_gwas_result(hdf5_file,result,mode='r'):
    f = tables.open_file(hdf5_file,mode)
    result_table = f.get_node('/phenotypes/%s' % result)
    data =  result_table[:]
    threshold = result_table._v_attrs['pval_threshold']
    f.close()
    return data,threshold


if __name__ == "__main__":
    sys.exit(main())
