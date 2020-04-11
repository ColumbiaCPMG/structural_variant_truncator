import sys
import csv
import re
import portion as P
import argparse
import intervaltree


def main():

    args = parseArgs(sys.argv)

    regions = loadGenomicCoordinatesFile(args.regions)
    sys.stderr.write('Regions File Processed\n')

    processPennCNVfile(args.variants, regions)

    params = {
        # 'disj': ' or ',
        'sep': '\t',
        'left_closed': '',
        'right_closed': '',
        # 'left_open': '..',
        # 'right_open': '..',
        # 'conv': lambda v: '"{}"'.format(v),
    }

    # for i in p:
    #     print("1\t", P.to_string(i, **params))

    # p = getMatchingIntervals('2', test_cnv2, regions)
    # print (2, p)
    # for i in p:
    #     print("2\t", P.to_string(i, **params))

    # for key in sorted(regions):
    #     print (key + "\t" + ";".join(sorted(regions[key])))


def loadGenomicCoordinatesFile(regions_file_path):
    regions = {}
    with open(regions_file_path, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')

        # skip headers
        row = next(reader)
        while row[0][0] == '#':
            row = next(reader)

        for i, row in enumerate(reader):

            # prepend 'chr' to the chromosome if it isn't already there
            chromosome = row[0]

            if not re.match('chr', chromosome): 
                chromosome = 'chr' + row[0] 
                

            

            pos = P.closed(int(row[1]), int(row[2]))

            if not chromosome in regions:
                # regions[chromosome] = [pos]
                regions[chromosome] = pos
            else:
                # regions[chromosome].append(pos)
                regions[chromosome] = regions[chromosome] | pos

            print('Line:' + str(i), end='\r')
        print('')

    return regions


def processPennCNVfile(variants_file_path, regions):
    try:
        variants_file = open(variants_file_path, 'r')
    except OSError:
        sys.stderr.write("Could not open/read file:" + variants_file)
        sys.exit()

    with variants_file:
        reader = csv.reader(variants_file, delimiter='\t')
        writer = csv.writer(sys.stdout, delimiter="\t")
        
        # skip header
        header = next(reader)
        header.append("Truncated")
        writer.writerow(header)

        for row in reader:
            chromosome, loc = row[0].split(':')
            start, end = loc.split('-')

            p = getMatchingIntervals(chromosome, P.closed(int(start), int(end)), regions)


            original_row = row.copy()
            original_row.append('original')
            writer.writerow(original_row)
            # row[0] = chr + ':' + start + '-' + end

            if p != P.empty():
                print('chromosome:', chromosome, p)


def getMatchingIntervals(chromosome, structural_variant, regions):
    intersections = []

    if chromosome not in regions:
        sys.stderr.write(
            "\n** Warninng Chr " + chromosome + " not found in Regions file.\n")
        return P.empty()

    for region in regions[chromosome]:
        # check to see if there is an overlap between the
        portion = region & structural_variant
        if not portion == P.empty():
            intersections = intersections + [portion]

    p = P.empty()
    for i in intersections:
        p = p | i

    return p


def parseArgs(args):
    parser = argparse.ArgumentParser(
        description='Truncate Structural Variants (SVs) to specified genomic regions.')
    parser.add_argument('-v', '--variants',
                        help='File containing variants in PennCNV format.')
    parser.add_argument('-r', '--regions',
                        help='File containing genomic regions in format "chromosome,start,end" (tab delimited also ok).')

    return parser.parse_args()


if __name__ == "__main__":
    main()
