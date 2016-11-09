import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import argparse
import subprocess
import collections
import math
from pybedtools import BedTool

def generate_gene_dict(beds, output_folder):
    """
    Generate a combined BED file and fill to get the genes
    :param beds: a list of all the BED files to be merged
    :param output_folder: folder to save any output files
    :return: file path to a BED file for the coordinates of each gene (as produced in the pipeline
    """

    if output_folder.endswith('/'):
        output_genes = output_folder + 'genes.bed'
        output_regions = output_folder + 'merged_regions.bed'
    else:
        output_genes = output_folder + '/genes.bed'
        output_regions = output_folder + '/merged_regions.bed'

    all_lines = []

    for bed in beds:
        f = open('/results/Analysis/projects/NBS/' + bed, 'r')
        lines = [line.strip('\n') for line in f.readlines()]
        f.close()

        for line in lines:
            if line.startswith('#'):
                continue
            else:
                all_lines.append(line)

    raw_regions = '\n'.join(all_lines)
    f = open ('/results/Analysis/projects/NBS/raw_regions.bed', 'w')
    f.write(raw_regions)
    f.close()

    whole = BedTool(raw_regions, from_string=True)
    whole_sorted =  whole.sort()
    whole_merged = whole_sorted.merge(nms=True)
    whole_merged.saveas(output_regions)

    command = 'python /results/Pipeline/SDGSPipeline/scripts/fill_bed.py --bed ' + output_regions + ' --out ' + output_genes

    try:
        subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as e:
        print(command)
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    f = open(output_genes)
    lines = [line.strip('\n') for line in f.readlines()]
    f.close()

    gene_dict = {}

    for line in lines:
        fields = line.split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        if chrom not in gene_dict:
            gene_dict[chrom] = collections.OrderedDict()
            gene_dict[chrom][name] = {'start': start, 'end': end}
        elif name not in gene_dict[chrom]:
            gene_dict[chrom][name] = {'start': start, 'end': end}
        else:
            print('gene duplicate ' + line + ':' + name)


    return gene_dict


def generate_bed_dict(beds, output_folder):
    """
    Generate a dictionary of bed files using their abbreviation as the key.
    The dictionary include the regions within the BED file and the coverage values will be added to this subsequently.
    :param beds: list of bed files
    :param output_folder: folder to save any output files
    :return: dictionary of BED files with the abbreviation as the key
    """
    bed_dict = {}
    for bed in beds:
        # find the bed abbreviation from the file on the server to name the result files
        #command = "grep " + bed + " /results/Analysis/MiSeq/MasterBED/abbreviated_bed_names.txt | cut -f2"
        try:
            # abv = subprocess.check_output(command, shell=True).replace('\n', '')
            # print(abv)
            # if not abv:
            #     print('No abv found for ' + bed)
            #     abv = bed
            abv = 'NBS1'
            bed_dict[abv] = {'name':bed, 'regions':{}}

            f = open('/results/Analysis/projects/NBS/' + bed, 'r')
            lines = [line.strip('\n') for line in f.readlines()]
            f.close()

            for line in lines:
                #if line == lines[0]:
                    #continue
                fields = line.split('\t')
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                #print name
                if chrom not in bed_dict[abv]['regions']:
                    bed_dict[abv]['regions'][chrom] = collections.OrderedDict()
                    bed_dict[abv]['regions'][chrom][name] = {'start':start, 'end':end, 'coverage':{}}
                elif name not in bed_dict[abv]['regions'][chrom]:
                    bed_dict[abv]['regions'][chrom][name] = {'start':start, 'end': end, 'coverage': {}}
                else:
                    print('region duplicate ' + line + ':' + bed)
        except subprocess.CalledProcessError as e:
            print('Error executing command: ' + str(e.returncode))
            exit(1)

    return bed_dict

def get_coverage(coverage_input, bed_dict):
    """
    Add the coverage values for each region to the BED dictionary
    :param coverage_input: list of files containing the coverage at each position in the broad panel
    :param bed_dict: dictionary of BED files and regions for generation of coverage graphs
    :return: bed_dict containing the coverage for each patient at each position
    """

    #generate coverage dictionary (used to populate BED dict)
    coverage = {}

    for bed in coverage_input:
        f = open(bed, 'r')
        lines = [line.strip('\n') for line in f.readlines()]
        f.close()

        for line in lines:
            #ignore header line
            #if line == lines[0] and bed != 'genes':
                #continue
            #print line
            fields = line.split('\t')
            chrom = fields[0]
            start = int(fields[1])
            #end position (fields[2]) not relevant
            cov = int(fields[3])
            try:
                if chrom not in coverage:
                    coverage[chrom] = {}
                    coverage[chrom][start] = {'cov':[cov,], 'name':''}
                elif start not in coverage[chrom]:
                    coverage[chrom][start] = {'cov':[cov,], 'name':''}
                else:
                    #add coverage value to existing list
                    covs = coverage[chrom][start]['cov']
                    covs.append(cov)
                    coverage[chrom][start]['cov'] = covs
            except KeyError: #start and end coordinates for the genes panel extend beyond the ROI for coverage
                pass

    #iterate through BED dictionary and add coverage values to relevant regions
    for abv in bed_dict.keys():
        for chrom in bed_dict[abv]['regions'].keys():
            for name in bed_dict[abv]['regions'][chrom].keys():
                start = bed_dict[abv]['regions'][chrom][name]['start']
                end = bed_dict[abv]['regions'][chrom][name]['end']
                #calculate average coverage and portions for region plots
                count_less_18 = 0
                count_18_30 = 0
                count_30_100 = 0
                count_more_100 = 0

                #list start<=i<end
                for i in range(start, end):
                    try:
                        cov = coverage[chrom][i]['cov']
                        coverage[chrom][i]['name'] = name
                        avg = sum(cov) / float(len(cov))
                        if avg < 18:
                            count_less_18 += 1
                        elif 18<= avg < 30:
                            count_18_30 += 1
                        elif 30 <= avg < 100:
                            count_30_100 += 1
                        else:
                            count_more_100 += 1

                        if i not in bed_dict[abv]['regions'][chrom][name]['coverage']:
                            bed_dict[abv]['regions'][chrom][name]['coverage'][i] = cov
                        else:
                            print('duplicate')

                        length = end - start
                        p_less_18 = count_less_18 / float(length)
                        p_18_30 = count_18_30 / float(length)
                        p_30_100 = count_30_100 / float(length)
                        p_more_100 = count_more_100 / float(length)

                        bed_dict[abv]['regions'][chrom][name]['p_less_18'] = p_less_18
                        bed_dict[abv]['regions'][chrom][name]['p_18_30'] = p_18_30
                        bed_dict[abv]['regions'][chrom][name]['p_30_100'] = p_30_100
                        bed_dict[abv]['regions'][chrom][name]['p_more_100'] = p_more_100

                    except KeyError:
                        pass

    return coverage, bed_dict

def plot_region_coverage(panel_name, panel_dict, output_folder):
    """
    Generate a plot for each panel showing the proportion of positions that have an arverage coverage in each of the
    pre-defined bins (less than 18X, between 18 and 30X, between 30 and 100X and greater than 100X)
    :param panel_name: abbreviation for the panel used in the pipeline
    :param panel_dict: dictionary for the specific sub-panel
    :param output_folder: location for graphs to be saved
    :return:
    """

    if output_folder.endswith('/'):
        output_file = output_folder + panel_name + '_coverage_plot.png'
    else:
        output_file = output_folder + '/' + panel_name + '_coverage_plot.png'

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
              'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
              'chrY']
    regions = panel_dict['regions']

    #fill arrays with the proportions in region order to be added to the graph
    names = []
    less_18 = []
    between_18_30 = []
    between_30_100 = []
    more_100 = []

    #get regions in chromosome order; regions should be in the order they were entered into the dictionary
    for chrom in chroms:
        try:#chromosome may not be in the panel
            for name in regions[chrom].keys():
                names.append(name)
                less_18.append(regions[chrom][name]['p_less_18'])
                between_18_30.append(regions[chrom][name]['p_18_30'])
                between_30_100.append(regions[chrom][name]['p_30_100'])
                more_100.append(regions[chrom][name]['p_more_100'])
        except KeyError:
            pass

    # split panels with more than 70 regions onto multiple lines
    if len(names) > 70:
        print(panel_name)
        graph_no = int(math.ceil(len(names) / 70.0))
        spacing = [i + 1 for i in range(70)]
        width = 0.75
        tick_pos = [i + (width / 2.) for i in spacing]
        f, subplots = plt.subplots(graph_no, 1, figsize=(19, 10 * graph_no))

        subs = range(graph_no)

        sub_less_18 = numpy.split(numpy.asarray(less_18), [i * 70 for i in subs])
        sub_18_30 = numpy.split(numpy.asarray(between_18_30), [i * 70 for i in subs])
        sub_30_100 = numpy.split(numpy.asarray(between_30_100), [i * 70 for i in subs])
        sub_more_100 = numpy.split(numpy.asarray(more_100), [i * 70 for i in subs])
        sub_regions = numpy.split(numpy.asarray(names), [i * 70 for i in subs])

        plot_no = 1
        for ax in subplots.flatten():
            if len(sub_less_18[plot_no]) < 70:
                a = sub_less_18[plot_no]
                b = sub_18_30[plot_no]
                c = sub_30_100[plot_no]
                d = sub_more_100[plot_no]
                # extend the array to make 70
                zeros = [0] * (70 - len(sub_less_18[plot_no]))
                print(len(a))
                print(zeros)
                a = numpy.append(a, zeros)
                b = numpy.append(b, zeros)
                c = numpy.append(c, zeros)
                d = numpy.append(d, zeros)
                print(len(a))
            else:
                a = sub_less_18[plot_no]
                b = sub_18_30[plot_no]
                c = sub_30_100[plot_no]
                d = sub_more_100[plot_no]

            ax.bar(spacing, a, width, color='r', label='Less than 18X')
            ax.bar(spacing, b, width, bottom=a, color='y', label='Between 18X and 30X')
            ax.bar(spacing, c, width, bottom=[i + j for i, j in zip(a, b)],
                   color='g', label='Between 30X and 100X')
            ax.bar(spacing, d, width,
                   bottom=[i + j + k for i, j, k in zip(a, b, c)], color='m',
                   label='Over 100X')

            ax.set_xticks(tick_pos)
            ax.set_xticklabels(sub_regions[plot_no], rotation='vertical')
            ax.set_yticks(numpy.arange(0, 1.1, 0.1))
            ax.set_xlim([min(tick_pos) - width, max(tick_pos) + width])
            if plot_no == 1:
                ax.set_title(panel_name)
            plot_no += 1

    else:
        print(panel_name)
        spacing = [i + 1 for i in range(len(names))]
        width = 0.75
        tick_pos = [i + (width / 2.) for i in spacing]
        f, ax1 = plt.subplots(1, figsize=(0.2 * len(names) + 5, 10))

        ax1.bar(spacing, less_18, width, color='r', label='Less than 18X')
        ax1.bar(spacing, between_18_30, width, bottom=less_18, color='y', label='Between 18X and 30X')
        ax1.bar(spacing, between_30_100, width, bottom=[i + j for i, j in zip(less_18, between_18_30)],
                color='g', label='Between 30X and 100X')
        ax1.bar(spacing, more_100, width,
                bottom=[i + j + k for i, j, k in zip(less_18, between_18_30, between_30_100)], color='m',
                label='Over 100X')
        print(tick_pos)
        ax1.set_xlim([min(tick_pos) - width, max(tick_pos) + width])
        plt.title(panel_name)
        plt.xticks(tick_pos, names, rotation='vertical')
        plt.yticks(numpy.arange(0, 1.1, 0.1))

    f.text(0.001, 0.5, 'Proportion of region at given depth', verticalalignment='center', rotation='vertical')
    f.text(0.5, 0.01, 'Region Name', verticalalignment='center')
    plt.legend(loc='upper right')
    f.tight_layout()
    plt.savefig(output_file)

def generate_gene_plots(gene_dict, coverage, output_folder):
    """
    Generate a figure for each gene showing the coverage stats for each position in each exon.
    The figures are split into a plot for each exon and show the average coverage at each position (red line), the
    maximum coverage at each position (green line) and the minimum coverage at each position (blue line)
    :param gene_dict: dictionary for the coverage values for each position in each gene
    :param coverage: dictionary containing the coverage information for the validation patients
    :param output_folder: location for graphs to be saved
    :return:
    """

#print coverage['chr1']['76190447']['cov']

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
              'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
              'chrY']
    for chrom in chroms:
        try:
            for gene in gene_dict[chrom].keys():
                #print gene
                if output_folder.endswith('/'):
                    output_file = output_folder + gene + '_plot.png'
                else:
                    output_file = output_folder + '/' + gene + '_plot.png'

                x = []
                maximums = []
                avg = []
                minimums = []

                start = gene_dict[chrom][gene]['start']
                end = gene_dict[chrom][gene]['end']
                # print start
                # print end
                plot_list = {}
                plot_number = 0
                new_plot = False
                intron_number = 0
                exon = ''
                index = 0
                for i in range(start, end - 1):
                    try:
                        cov = coverage[chrom][i]['cov']
                        if new_plot:
                            print('plot added')
                            region_max = numpy.array(maximums)
                            region_min = numpy.array(minimums)
                            region_avg = numpy.array(avg)
                            region_length = numpy.array(x)
                            maximums = []
                            minimums = []
                            avg = []
                            x = []
                            index = 0
                            plot_list[plot_number] = {'title': exon, 'length': max(region_length), 'x': region_length,
                                                      'max': region_max,
                                                      'min': region_min, 'avg': region_avg}
                            plot_number += 1
                            new_plot = False
                        #print coverage[chrom][i]['name']
                        region_name = coverage[chrom][i]['name']
                        #print region_name
                        exon = region_name.split('_')[1]
                        intron_number += 1
                        cov_mean = numpy.mean(cov)
                        cov_max = max(cov)
                        cov_min = min(cov)
                        maximums.append(cov_max)
                        minimums.append(cov_min)
                        avg.append(cov_mean)
                        x.append(index)
                        index += 1
                    except KeyError:
                        if intron_number > 0:
                            new_plot = True

                if new_plot:
                    print('plot added')
                    region_max = numpy.array(maximums)
                    region_min = numpy.array(minimums)
                    region_avg = numpy.array(avg)
                    region_length = numpy.array(x)
                    plot_list[plot_number] = {'title': exon, 'length': max(region_length), 'x': region_length,
                                              'max': region_max,
                                              'min': region_min, 'avg': region_avg}
                    plot_number += 1

                if plot_number > 5:
                    row_number = int(math.ceil(plot_number / 5.))
                    f, subplots = plt.subplots(row_number, 5, sharey='all')
                    f.set_size_inches(30, 6 * row_number)
                else:
                    print(plot_number)
                    f, subplots = plt.subplots(1, plot_number, sharey='all')
                    f.set_size_inches(30, 6)
                f.text(0.001, 0.5, 'Coverage', verticalalignment='center', rotation='vertical')
                if len(plot_list.keys()) == 1:
                    p = plot_list[0]
                    subplots.plot(p['x'], p['max'], 'go', p['x'], p['avg'], 'ro', p['x'], p['min'], 'bo')
                    subplots.text(0.01, 0.99, p['title'], horizontalalignment='right', verticalalignment='top',
                                  transform=subplots.transAxes)
                else:
                    reverse = False
                    if float(plot_list[0]['title'].replace('Ex', '')) > float(plot_list[1]['title'].replace('Ex', '')):
                        reverse = True
                    i = 0
                    for ax in subplots.flatten():
                        if reverse:
                            index = len(plot_list) - (i + 1)
                        else:
                            index = i
                        try:
                            p = plot_list[index]
                            ax.plot(p['x'], p['max'], 'go', p['x'], p['avg'], 'ro', p['x'], p['min'], 'bo')
                            ax.text(1.0, 0.99, p['title'], horizontalalignment='right', verticalalignment='top',
                                    transform=ax.transAxes)
                            ax.set_aspect('auto')
                        except KeyError:
                            ax.axis('off')
                        i += 1
                plt.tight_layout()

                plt.savefig(output_file)
                plt.close()
        except KeyError:
            pass

def main():
    """
    This script creates graphs showing a breakdown of coverage across the list of patients provided in the arguments
    (-d) for the BED files listed (-b).
    Two sets of graphs are generated: proportion of regions covered at pre-defined depths and a breakdown of each exon
    in each gene to help highlight any areas that may cause problems when the panel goes into service.
    The graphs are saved in the specified folder (-o) and named either by the panel or gene.
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', help='Comma separated list of coverage files (full path)', required=True)
    parser.add_argument('-b', help='Comma separated list of BED files (file names only)', required=True)
    parser.add_argument('-o', help='Directory for graph output', required=True)

    args = parser.parse_args()

    coverage_input = args.d.split(',')
    beds = args.b.split(',')
    output_folder = args.o

    #add BED tools to the path so the python library can be used
    command = 'export PATH=${PATH}:/results/Pipeline/program/bedtools-2.17.0/bin'
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print(command)
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    #generate dictionaries
    bed_dict = generate_bed_dict(beds, output_folder)
    coverage, bed_dict = get_coverage(coverage_input, bed_dict)
    gene_dict = generate_gene_dict(beds,output_folder)

    #generate graphs
    for bed in bed_dict.keys():
        plot_region_coverage(bed, bed_dict[bed], output_folder)

    generate_gene_plots(gene_dict, coverage, output_folder)

if __name__ == '__main__':
    main()