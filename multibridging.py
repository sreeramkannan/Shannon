import time, sys
from mbgraph import *
import pdb
K_value = 24

def extract_reads(filename):
    """Extract all reads from fasta file FILENAME
    and return them as a list of strings.
    """
    with open(filename) as f:
        return [l[:-1].upper() for l in f if l[0] != '>']

def set_read_length(filename):
    with open(filename) as f:
        for l in f:
            if l[0] != '>':
                Read.L = len(l) - 1
                return

def load_reads(filename, double_stranded, record_kmers):
    """Load FASTA reads.
    """
    print("Loading reads...")
    print(time.asctime())
    for line in extract_reads(filename):
        Read.add_read(line, double_stranded, record_kmers)

def load_mated_reads(file_1, file_2, double_stranded):
    def pair(r1, r2):
        r1.mate_pair = 1
        r2.mate_pair = 2
        r1.mate = r2
        r2.mate = r1

    #There must be a newline at the end of the file!
    print("Building mate-paired graph...")
    Read.MATED_READS = True
    with open(file_1) as f1, \
        open(file_2) as f2:
        gen1 = (l[:-1].upper() for l in f1 if l[0] != '>')
        gen2 = (l[:-1].upper() for l in f2 if l[0] != '>')
        for line1 in gen1:
            line2 = next(gen2)
            r1 = Read.add_read(line1, double_stranded)
            r2 = Read.add_read(Read.reverse_complement(line2), double_stranded)
	    #r2 = Read.add_read(line2,double_stranded)
            if double_stranded:
                r1a, r1b = r1
                r2a, r2b = r2
                pair(r1a, r2a)
                pair(r2b, r1b)
            else:
                pair(r1, r2)

def load_cpp(files, nodes_file, edges_file, double_stranded):
    """Loads condensed files from C++ implementation.
    FILES are the read files, NODES/EDGES are the node/edge files.
    """
    setup(files)
    clear()

    print(time.asctime())
    print("Loading graph from C++ output.")

    if len(files) == 2:
        load_mated_reads(files[0], files[1], double_stranded)
    else:
        load_reads(files[0], double_stranded, False)
    Node.kmer_dict = {}
    Node.nodes = []

    nodes = {}

    with open(nodes_file) as f:
        for line in f:
            node_id, bases, prevalence = line.split()
            n = Node(bases)
            n.id = int(node_id)
            n.prevalence = int(round(float(prevalence)))
            nodes[node_id] = n

    with open(edges_file) as f:
        for line in f:
            id1, id2 = line.split()
            weight = Read.K - 1 #FIXME: Read this from edges file
            nodes[id1].link_to(nodes[id2], int(weight))

def load_two_kmers(two_kmer_file):
    two_kmers = {}

    with open(two_kmer_file) as f:
        for line in f:
            two_kmer, count = line.split()
            two_kmers[two_kmer] = int(count)

    return two_kmers



def load_kplus3mers(kplus3mers_file):
    kplus3mers = set()
    with open(kplus3mers_file) as f:
        for line in f:
            kplus3mers.add(line.split()[0])
    return kplus3mers


def setup(files):
    """Perform some setup functions.
    """
    set_read_length(files[0])
    Read.K = K_value
    Node.SIZE_THRESHOLD = Read.L

def condense(files, double_stranded = False):
    """Runs the algorithm up to condensing on input files FILES (1 file for single-ended reads,
    2 for paired end reads), writing the results to OUTPUT_DIR.
    """
    setup(files)
    clear()
    print(time.asctime())
    if len(files) == 2:
        load_mated_reads(files[0], files[1], double_stranded)
    else:
        load_reads(files[0], double_stranded, True)

    Node.kmer_dict = {}

    print(time.asctime())
    print("Condensing graph...")
    print(str(len(Node.nodes)) + " initial nodes.")

    Node.condense_all()
    print(str(len(Node.nodes)) + " nodes after condensing.")

def log(text):
    print("%s: %s" % (time.asctime(), text))

def run(output_dir, error_correction = False, compute_fringes = False,PE = False):
    """Continues the algorithm from condensing, writing the results to OUTPUT_DIR.
    """
    if error_correction:
        Node.destroy_suspicious()
        print(str(len(Node.nodes)) + " nodes after destroying suspicious nodes.")
        Node.collapse_all()
        print(str(len(Node.nodes)) + " nodes after collapsing similar nodes.")
    if 0: #Replace by flag
	two_kmers_file = '/data/sreeramk/Full_Assembler/S15_Quorum_algo_input/2Kmer.dict'	
        #two_kmers_file = '/data/sreeramk/Full_Assembler/Snyder_group2/2Kmer_algo_input/2Kmer.dict'
	two_kmers = load_two_kmers(two_kmers_file)
	kplus3mers_file = '/data/sreeramk/Full_Assembler/S15_Quorum_algo_input/kplus3mer.dict' 
        #kplus3mers_file = '/data/sreeramk/Full_Assembler/Snyder_group2/2Kmer_algo_input/kplus3mer.dict'
	kplus3mers = load_kplus3mers(kplus3mers_file)
	#pdb.set_trace()
	info = Node.remove_edges_erroneous(Read.K, two_kmers, kplus3mers)
        print("No of edges, no < K+3 deleted, No < K+3 accepted, No > K+3 accepted, No. > 2K accepted, No.> 2K  destroyed, No < 2K and not supported by contig: " + str(info[0]) + ", " + str(info[1]) + ", " + str(info[2])+', '+ str(info[3]) + ', '+str(info[4]) + ', ' + str(info[5]) + ', ' + str(info[6]))
	#raw_input()
    log("Bridging graph.")
    Read.find_bridging_reads()
    Node.bridge_all()
    Node.condense_all()

    print(str(len(Node.nodes)) + " nodes after bridging.")
    log("Finding copy counts.")
    Node.find_approximate_copy_counts()
    #Node.disregard_loops()

    log("Breaking cycles.")
    #Node.condense_all()
    Node.remove_destroyed()
    #Node.aggressive_break_cycles()
    Node.break_cycles(False)
    Node.condense_all()

    log("Finding approximate copy counts.")
    #Node.find_approximate_copy_counts()
    log("Finding copy counts.")
    Node.find_copy_counts(False, False)
    #Node.find_approximate_copy_counts()
    known_paths()
    if PE:
    	log("Finding mate pairs.")
    	Read.find_mate_pairs()

    print(str(len(Node.nodes)) + " final nodes.")
    log("Exporting graph.")

    output_components(output_dir)  # Be careful - this will destroy the graph
    log("Done.")

def output_components(output_dir):
    """Outputs a textual representation of the graph as components, one
    component per set of node/edge/known-paths files. Destroys the graph.
    """
    if output_dir is None:
        return

    component = 0
    paths_by_start = {}
    for path in Read.known_paths:
        if path[0] not in paths_by_start:
            paths_by_start[path[0]] = []
        paths_by_start[path[0]].append(path)

    with open(output_dir + '/single_nodes.txt', 'w', 0) as single_file:
        single_file.write("ID\tBases\tCopycount\tNormalization\n")

        for source_node in Node.nodes:
            if hasattr(source_node, 'destroyed'):
                continue
            with open(output_dir + '/nodes'+str(component)+'.txt', 'w', 0) as nodefile, \
                open(output_dir + '/edges'+str(component)+'.txt', 'w', 0) as edgefile, \
                open(output_dir + '/paths'+str(component)+'.txt', 'w', 0) as pathfile:
                component_nodes, component_edges = source_node.add_component()
                component_nodes = (Node.topological_sort(component_nodes))
		component_node_set = set(component_nodes)

                if len(component_nodes) == 1:
                    source_node.hash = -1
                    single_file.write(source_node.to_string())
                    source_node.destroyed = True
                    continue

                node_hash = 0
                nodefile.write("ID\tBases\tCopycount\tNormalization\n")
                pathfile.write("ID1\tID2\tEtc.\n")
                for node in component_nodes:
                    node.hash = node_hash
                    node_hash += 1
                    nodefile.write(node.to_string())
                    node.destroyed = True

		for node in component_nodes:
		    if not paths_by_start.get(node):
			continue;
                    for path in paths_by_start[node]:
                        path = [str(n.hash) for n in path]
                        pathfile.write("\t".join(path) + "\n")

                edgefile.write("InID\tOutID\tWeight\tCopycount\tNormalization\n")
                for edge in component_edges:
                    edgefile.write(edge.to_string())

                component += 1


def main():
    if len(sys.argv) == 1:
        arguments = ['', '-c', 'input/nodes.txt', 'input/edges.txt',
                     'input/reads.fasta', 'output_aggressive']
    else:
        arguments = sys.argv
    sys.setrecursionlimit(100000)
    names = []
    double_stranded = False
    error_correction = False
    compute_fringes = False
    cpp = False
    corr_2kmer = True
    for arg in arguments[1:]:
        if arg == '-d':
            double_stranded = True
        elif arg == '-e':
            error_correction = True
        elif arg == '-f':
            compute_fringes = True
        elif arg == '-c':
            cpp = True
	elif arg == '-2k':
	    corr_2kmer = True
	elif len(arg)>=8 and arg[:6]=='--kmer':
	    K_value = int(arg[7:])
        else:
            names.append(arg)

    names, output_dir = names[:-1], names[-1]

    if cpp:
        node_file, edge_file = names[:2]
        names = names[2:]
	PE = True if len(names)==2 else False
        load_cpp(names, node_file, edge_file, double_stranded)
    else:
	PE = True if len(names)==2 else False
        condense(names,double_stranded)

    run(output_dir, error_correction, compute_fringes,PE)

if __name__ == '__main__':
    main()
