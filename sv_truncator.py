import sys
import csv
import re
import argparse
from intervaltree import Interval, IntervalTree


def main():

    args = parseArgs(sys.argv)

    trees = loadGenomicCoordinatesFile(args.regions)
    sys.stderr.write('Regions File Processed\n')

    processStructuralVariantsFile(args.variants, trees)
    sys.stderr.write('Conversion Complete\n')


# Load the file which contains the genomic regions to restrict the variants to.
# Typicall this would be an list of exons.
def loadGenomicCoordinatesFile(regions_file_path):
    # regions = {}
    trees = {}

    try:
        regions_file = open(regions_file_path, 'r')
    except OSError:
        sys.stderr.write("Could not open/read file:" + regions_file_path)
        sys.exit()
    except TypeError:
        sys.stderr.write("Error: You must specify a genomic regions file using -r.\n\n")
        sys.exit()

    with regions_file:
        reader = csv.reader(regions_file, delimiter='\t')

        # skip headers
        # row = next(reader)
        # while row[0][0] == '#':
        #     row = next(reader)

        for i, row in enumerate(reader):

            # prepend 'chr' to the chromosome if it isn't already there
            chromosome = row[0]
            if not re.match('chr', chromosome): 
                chromosome = 'chr' + row[0] 
            
            # if for some reason the start and end coordinates are not ints, 
            # skip to next line
            if not str.isdigit(row[1]) or not str.isdigit(row[2]):
                continue

            # Get the start and end coordinates.
            # Add 1 to the end coordinate if they are the same, since if they are the
            # same, that's an empty coordinate and we want an Interval of at least 1 basepair
            start, end = int(row[1]), int(row[2])
            if start == end:
                end += 1
            
            pos = Interval(start, end)

            # add a tree for each chromosome to the trees dictionary
            # if not chromosome in regions:
            if not chromosome in trees:
                # regions[chromosome] = [pos] 
                # regions[chromosome] = pos
                t = IntervalTree()
                t[start:end] = ''
                trees[chromosome] = t
            # if the key for that chromosome already exist, append the interval
            # the data it contains doesn't matter, so use and empty string
            else:
                # regions[chromosome].append(pos)
                trees[chromosome][start:end] = ''
                # regions[chromosome] = regions[chromosome] | pos

            # Show wich line has been processed in the terminal
            sys.stderr.write('Line:' + str(i) + '\r')

        sys.stderr.write('\n')

    return (trees)


def processStructuralVariantsFile(variants_file_path, trees):
    try:
        variants_file = open(variants_file_path, 'r')
    except OSError:
        sys.stderr.write("Could not open/read file:" + variants_file)
        sys.exit()
    except TypeError:
        sys.stderr.write("Error: You must specify a Variants file using -v.\n\n")
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

            # skip the line if start and end aren't digits
            if not str.isdigit(start) or not str.isdigit(end):
                continue

            q = getMatchingIntervalsFromTree(chromosome, Interval(int(start), int(end), 0), trees)

            original_row = row.copy()
            original_row.append('original')
            writer.writerow(original_row)
        
            if q and not q.is_empty():
                row[0] = chromosome + ':' + str(q.begin()) + '-' + str(q.end())
                row.append('truncated')
                writer.writerow(row)


def getMatchingIntervalsFromTree(chromosome, structural_variant, trees):
    it = IntervalTree()

    if chromosome not in trees:
        # sys.stderr.write(
            # "** Warning Chr " + chromosome + " not found in Regions file.\n")
        return it

    overlap = IntervalTree(trees[chromosome].overlap(structural_variant.begin, structural_variant.end))

    if not overlap.is_empty():
        lower = max(structural_variant.begin, overlap.begin())
        upper = min(structural_variant.end, overlap.end())
        i = Interval(lower, upper)
        it.add(i)

    return it


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
