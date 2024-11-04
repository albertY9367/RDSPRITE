import pysam
import re
import gzip
import os
import sys
from collections import defaultdict, Counter
from tqdm import tqdm


class Position:
    """This class represents a genomic position, with type of nucleic acid (RNA or DNA)

    Methods:
    - to_string(): Returns a string representation of this position in the form
      "R/DPM(feature)_chrX:1000"
    """
    def __init__(self, read_type, feature, chromosome, start_coordinate, end_coordinate):
        self._type = read_type
        self._feature = feature
        self._chromosome = chromosome
        self._start_coordinate = start_coordinate
        self._end_coordinate = end_coordinate

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return (self._type == other._type and
                self._feature == other._feature and
                self._chromosome == other._chromosome and
                self._start_coordinate == other._start_coordinate and
                self._end_coordinate == other._end_coordinate)

    def __hash__(self):
        return hash((self._type, self._feature, self._chromosome, 
                     self._start_coordinate, self._end_coordinate))

    def to_string(self):
        try:
            out = self._type + "[" + self._feature + "]" + "_" + \
                self._chromosome + ":" + \
                str(self._start_coordinate) + "-" + str(self._end_coordinate)
        except:
            print(self._type, self._feature, self._chromosome)
            print('Elements are not as expect!')
            sys.exit()
        return out

    def score(self):
        if self._type == 'RPM':
            if self._chromosome.startswith('chr'):
                return 1
            else:
                return 2
        else:
             return 3
                
           
  

class UMIs:
    """This class hold the UMI

    Methods:
        to_string: Returns a sting of the UMI
    """

    def __init__(self):
        self._umi = []

    def add_umi(self, umi):
        self._umi.append(umi)

    def unique(self):
        return '\t'.join(set(self._umi))

    def to_string(self):
        return '\t'.join(self._umi)

class Cluster:
    """This class represents a barcoding cluster as a collection of genomic
    positions.

    The underlying data structure is a set, so duplicate positions are
    discarded.

    Methods:
    - add_position(position): Adds a genomic position to this cluster

    - size(): Returns the number of reads or positions in this cluster

    - to_string(): Returns a string representation of this cluster as a
      tab-delimtited series of positions. See Position#to_string for how
      positions are represented as strings.

    - to_list(): Returns the Position class as a list (similar to_string())
    """

    def __init__(self):
        self._positions = set()

    def __iter__(self):
        return iter(self._positions)

    def add_position(self, position):
        self._positions.add(position)

    def size(self, read_type=None):
        if read_type == None:
            return len(self._positions)
        else:
            return sum([1 if pos._type == read_type else 0 for pos in self._positions])

    def to_string(self):
        positions_sorted = sorted(list(self._positions), key = lambda x: x.score())
        position_strings = [position.to_string() for position in positions_sorted]
        #position_ordered = sorted(position_strings, key = lambda x: x.startswith('DPM'))
        return "\t".join(position_strings)

    def to_list(self):
        position_strings = [position.to_string() for position in self._positions]
        position_ordered = sorted(position_strings, key = lambda x: x.startswith('DPM'))
        return position_ordered

    def count_type(self):
        rna_dna_count = Counter([position._type for position in self._positions])
        return rna_dna_count

class Clusters:
    """This class represents a collection of barcoding clusters.

    Methods:
    - get_cluster(barcode): Returns the cluster that corresponds to the given
      barcode. If the cluster does not exist, it is initialized (with zero
      positions), and this empty cluster is returned.

    - get_items(): iterate over clusters dictionary yielding keys and values

    - add_position(barcode, position): Adds the position to the cluster
      that corresponds with the given barcodes

    - to_strings(): Returns an iterator over the string representations of all
      of the contained clusters.

    - remove_cluster(barcode): Removes a cluster with the specified barcode

    - add_umi(barcode, umi): Used for fastq to cluster format, will used umi
      instead of position coordinates

    - unique(): keep only unique cluster entries

    - make_lookup(): make a lookup table for converting cluster back into bam
    """
    def __init__(self):
        self._clusters = {}

    def __iter__(self):
        return iter(self._clusters.values())

    def __getitem__(self, barcode):
        return self._clusters[barcode]

    def get_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = Cluster()
        return self._clusters[barcode]

    def get_items(self):
        return self._clusters.items()

    def add_position(self, barcode, position):
        self.get_cluster(barcode).add_position(position)

    def to_strings(self):
        for barcode, cluster in self._clusters.items():
            #yield barcode + '\t' + cluster.to_string()
            yield barcode + '\t' + cluster.to_string()

    def remove_cluster(self, barcode):
        del self._clusters[barcode]

    def get_umi_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = UMIs()
        return self._clusters[barcode]

    def add_umi(self, barcode, umi):
        self.get_umi_cluster(barcode).add_umi(umi)

    def unique(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.unique()

    def make_lookup(self):
        lookup = defaultdict(set)
        for barcode, cluster in self._clusters.items():
            lookup[barcode].update(cluster.to_list())
        return lookup

    def make_stripped_lookup(self):
        lookup = defaultdict(set)
        for barcode, cluster in self._clusters.items():
            barcode_strip = barcode.split('.')[:-1]
            barcode_merge = '.'.join(barcode_strip)
            lookup[barcode_merge].update(cluster.to_list())
        return lookup




def anno_type(element):
    '''
    Extract annotation type
    '''
    anno_type = element.split('.')[-1]
    if anno_type == 'exon':
        return 1
    elif anno_type == 'intron':
        return 2
    elif anno_type == 'repeat':
        return 3
    elif anno_type == 'none':
        return 4
    else:
        return 0


def order_annotation(anno):
    '''
    Enforce exon, intron, repeat, none order
    '''
    anno_list = anno.split(';')
    sorted_anno = sorted(anno_list, key=anno_type)
    return ';'.join(sorted_anno)

def fix_annotation(anno):
    d = {'exon':[], 'intron':[], 'repeat':[], 'none':[]}
    anno_list = anno.split(';')
    for a in anno_list:
        aname, atype = a.rsplit('.', 1)
        d[atype].append(aname)
    expanded_list = []
    for key,value in d.items():
        if len(value) == 0:
            expanded_list.append('*.' + key)
        elif len(value) == 1:
            expanded_list.append(value[0] + '.' + key)
        else:
            expanded_list.append(','.join([name + '.' + key for name in value]))
            #expanded_list.append('Ambiguous.none.' + key)
    sorted_anno = sorted(expanded_list, key=anno_type)[:-1]
    if len(sorted_anno) != 3:
        print('Too many annotations')
        sys.exit()
    anno_join = ';'.join(sorted_anno)
    anno_clean = anno_join.replace('*.exon', '').replace('*.intron', '').replace('*.repeat', '')
    return anno_clean
    

def expand_annotation(anno):
    print(anno)
    anno_list = anno.split(';')
    if '.exon' not in anno:
        anno_list.append('*.exon')
    if '.intron' not in anno:
        anno_list.append('*.intron')
    if '.repeat' not in anno:
        anno_list.append('*.repeat')
    sorted_anno = sorted(anno_list, key=anno_type)
    if '.none' in anno:
        sorted_anno = sorted_anno[:-1]
    if len(sorted_anno) != 3:
        print('Too many annotations')
        sys.exit()
    anno_join = ';'.join(sorted_anno)
    anno_clean = anno_join.replace('*.exon', '').replace('*.intron', '').replace('*.repeat', '')
    return anno_clean

def get_clusters(bamfile, num_tags):
    """Parses a BAM file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each BAM record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3][DPM/RPM]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.

    Output cluster barcodes will be in the format of:
    Tag1.Tag2.Tag3.SampleName
    """
    
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    rpm_counts = 0
    dpm_counts = 0
    annotation_flag = True

    for bam in bamfile:
        #get sample name from bamfile
        file_name = os.path.basename(bam)
        sample_name = file_name.split('.')[0]
        try:
            with pysam.AlignmentFile(bam, "rb") as f:
                for read in f.fetch(until_eof = True):
                    name = read.query_name
                    match = pattern.search(name)
                    barcode = list(match.groups())
                    strand = '+' if not read.is_reverse else '-'
                    #strip RPM DPM from barcode
                    if 'RPM' in barcode:
                        rpm_counts = rpm_counts + 1
                        #get featureCounts annotation
                        if read.has_tag('XT'):
                            gene_anno = read.get_tag('XT')
                            anno = fix_annotation(gene_anno) + ";" + strand
                        elif read.has_tag('XS'):
                            gene_anno = read.get_tag('XS')
                            anno = fix_annotation(gene_anno) + ";" + strand
                        else:
                            anno = ";;;" + strand
                            annotation_flag = False
                        position = Position('RPM', anno, read.reference_name,
                                            read.reference_start, read.reference_end)
                        barcode.remove('RPM')
                    elif 'DPM' in barcode:
                        dpm_counts = dpm_counts + 1
                        position = Position('DPM', strand, read.reference_name,
                                            read.reference_start, read.reference_end)
                        barcode.remove('DPM')
                    else:
                        raise Exception('RPM or DPM not present in full barcode')
                    barcode.append(sample_name)
                    barcode_str = ".".join(barcode)
                    clusters.add_position(barcode_str, position)
        except ValueError:
            print('BAM file provided is not a BAM or is empty!')
    
    if annotation_flag == False:
        print("No annotations found for rpm. Writing clusters with only strand information")
    
    print("Total DPM: ", dpm_counts)
    print("Total RPM: ", rpm_counts)
    return clusters


def get_clusters_fastq(fastqfile, num_tags, orientation, umi_len):
    """Parses a fastq file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each fastq record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """
    #strip RPM DPM from barcode
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    count = 0
    with file_open(fastqfile) as f:
        for name, seq, thrd, qual in fastq_parse(f):
            #extract UMI and place as position data
            match = pattern.search(name)
            barcode = list(match.groups())
            if orientation == 'r1': #read 1 orientation, UMI on left
                umi = seq[:umi_len]
                # umi = seq
            elif orientation == 'r2': #read 2 orientation, UMI on right
                umi = seq[-umi_len:]
                # umi = seq
            else:
                raise Exception('Incorrect orientation specified')

            barcode_str = "|".join(barcode[:num_tags])
            clusters.add_umi(barcode_str, umi)

            count += 1
    print('Reads parsed:', count)
    return clusters


def write_clusters_to_file(clusters, outfile, unique=False):
    """Writes a Clusters object to a file"""

    count = 0
    #filename = os.path.basename(outfile)
    #sample_name = filename.split('.')[0]
    with open(outfile, 'w') as f:

        if unique:
            for cluster_string in clusters.unique():
                f.write(cluster_string)
                f.write("\n")
                count += 1
        else:
            for cluster_string in clusters.to_strings():
                #cluster_string = cluster_string.replace('BARCODE', str(count) + '.' + sample_name)
                f.write(cluster_string)
                f.write("\n")
                count += 1
    print('Number of clusters written: ',count)


def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:

        linecount += 1
        if linecount % 4 == 1:
            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()
            assert name.startswith('@'),\
                   "ERROR: The 1st line in fastq element does not start with '@'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 2:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()
        elif linecount % 4 == 3:
            try:
                thrd = line.decode('UTF-8').rstrip()
            except AttributeError:
                thrd = line.rstrip()
            assert thrd.startswith('+'),\
                   "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual),\
                    "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)

            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4



def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f

def parse_cluster(c_file):
    '''
    Parse cluster file

    Args:
        c_file(str): input path of cluster file
    '''

    total_reads = 0
    clusters = Clusters()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')

    with file_open(c_file) as c:
        for line in tqdm(c):
            barcode, *reads = line.decode('utf-8').rstrip('\n').split('\t')
            for read in reads:
                total_reads += 1
                try:
                    match = pattern.search(read)
                    read_type, feature, chrom, start, end = match.groups()
                    position = Position(read_type, feature, chrom, start, end)
                    clusters.add_position(barcode, position)
                except:
                    print(read)
                    raise Exception('Pattern did not match above printed string')
    print('Total cluster reads:', total_reads)
    return(clusters)



