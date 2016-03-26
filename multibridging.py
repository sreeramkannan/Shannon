import time, sys, pdb
from mbgraph import *

def extract_reads(filename, weighted):
    """Extract all reads from fasta file FILENAME
    and return them as a list of (weight, string) tuples.
    """
    last_weight = None
    reads = []
    with open(filename) as f:
        for line in f:
            if line[0] == '>':
                data = line.split()
                if weighted:
                    last_weight = float(data[-1])
                else:
                    last_weight = 1.0
            else:
                reads.append((last_weight, line[:-1].upper()))
    return reads

def load_reads(filename, double_stranded, weighted):
    """Load FASTA reads.
    """
    log("Loading reads.")
    for read in extract_reads(filename, weighted):
        Read.add_read(read, double_stranded)

def load_mated_reads(file_1, file_2, double_stranded):
    #Does not allow weighted reads.
    def pair(r1, r2):
        r1.mate_pair = 1
        r2.mate_pair = 2
        r1.mate = r2
        r2.mate = r1

    #There must be a newline at the end of the file!
    log("Loading mated reads.")
    Read.MATED_READS = True
    with open(file_1) as f1, open(file_2) as f2:
        gen1 = (l[:-1].upper().strip() for l in f1 if l[0] != '>')
        gen2 = (l[:-1].upper().strip() for l in f2 if l[0] != '>')
        for line1 in gen1:
            line2 = next(gen2)
            r1 = Read.add_read((1.0, line1), double_stranded)
            r2 = Read.add_read((1.0, Read.reverse_complement(line2)), double_stranded)

            if double_stranded:
                r1a, r1b = r1
                r2a, r2b = r2
                pair(r1a, r2a)
                pair(r2b, r1b)
            else:
                pair(r1, r2)

def load_cpp(node_file, edge_file):
    """Loads condensed files from C++ implementation.
    FILES are the read files, NODES/EDGES are the node/edge files.
    Will implicitly set
    """
    log("Loading nodes from C++ output.")

    nodes = {}
    with open(node_file) as f:
        for line in f:
            node_id, bases, prevalence = line.split()
            n = Node(bases)
            n.id = int(node_id)
            n.prevalence = round(float(prevalence))
            nodes[node_id] = n

    with open(edge_file) as f:
        for line in f:
            id1, id2 = line.split()
            weight = Read.K - 1
            nodes[id1].link_to(nodes[id2], int(weight))
def load_jellyfish(node_file, edge_file):
    """Loads condensed files from Jelyfish.
    """
    log("Loading nodes from K-mer files.")

    nodes = {}
    with open(node_file) as f:
        for line in f:
            bases, prevalence = line.split()
            if bases in nodes:
                n = nodes[bases]
                n.prevalence += round(float(prevalence))
            else:
                n = Node(bases)
                n.prevalence = round(float(prevalence))
                nodes[bases] = n

    with open(edge_file) as f:
        for line in f:
            bases, prevalence = line.split()
            k1, k2 = bases[:-1], bases[1:]
            weight = Read.K - 1
            e = nodes[k1].link_to(nodes[k2], int(weight))
            e.copy_count = round(float(prevalence))

def load_single_jellyfish(edge_file):
    """Loads condensed files from Jellyfish.
    """
    log("Loading nodes from K-mer files.")

    nodes = {}
    # with open(node_file) as f:
    #     for line in f:
    #         bases, prevalence = line.split()
    #         n = Node(bases)
    #         n.prevalence = round(float(prevalence))
    #         nodes[bases] = n

    with open(edge_file) as f:
        for line in f:
            bases, prevalence = line.split()
            assert Read.K == len(bases) -1 
            k1, k2 = bases[:-1], bases[1:]
            if k1 not in nodes:
                n1 = Node(k1);  nodes[k1] = n1
            if k2 not in nodes:
                n2 = Node(k2);  nodes[k2] = n2
            weight = Read.K - 1
            e = nodes[k1].link_to(nodes[k2], int(weight))
            e.copy_count = round(float(prevalence))

    for node in Node.nodes:
        node.prevalence = sum(e.weight for e in node.out_edges)



def setup(read_file):
    """Perform some setup functions.
    """
    with open(read_file) as f:
        f.readline()
        Read.L = len(f.readline()) - 1
    Node.SIZE_THRESHOLD = Read.L
    clear()

def log(text):
    print("%s: %s" % (time.asctime(), text))

def run(output_dir, error_correction = False, compute_fringes = False):
    """Continues the algorithm starting with condensing, writing the results to OUTPUT_DIR.
    """
    log("Condensing.")
    Node.condense_all()
    no_nodes = max(len(Node.nodes),1)
    total_out_edges =sum(len(n.out_edges) for n in Node.nodes)
    total_no_bases = sum(len(n.bases) for n in Node.nodes)
    log("After condensing: No of nodes:" + str(no_nodes))
    log("Total out edges per node:" + str(float(total_out_edges)/no_nodes))
    log("Total no of bases per node:" + str(float(total_no_bases)/no_nodes))

    log(str(len(Node.nodes)) + " nodes after condensing.")

    if error_correction:
        Node.destroy_suspicious()
        log(str(len(Node.nodes)) + " nodes after destroying suspicious nodes.")
        Node.collapse_all()
        log(str(len(Node.nodes)) + " nodes after collapsing similar nodes.")

    Read.find_bridging_reads()
    log("Found bridging reads.")
    Node.bridge_all()
    Node.condense_all()

    log(str(len(Node.nodes)) + " nodes after bridging.")
    log("Finding copy counts.")
    Node.find_approximate_copy_counts()
    Node.disregard_loops()

    log("Breaking cycles.")
    Node.condense_all()
    Node.remove_destroyed()
    Node.break_cycles(False)

    log("Finding approximate copy counts.")
    Node.find_approximate_copy_counts()
    log("Finding known paths.")
    known_paths()
    log("Finding mate pairs.")
    Read.find_mate_pairs()

    #construct_reads()
    #log("Rebridging graph.")
    #Read.find_bridging_reads()
    #Node.bridge_all()
    #Node.condense_all()

    #log("Finding copy counts.")
    #Node.find_approximate_copy_counts()
    #Node.disregard_loops()
    #log("Finding known paths.")
    #known_paths()
    #Read.find_mate_pairs()

    log(str(len(Node.nodes)) + " final nodes.")
    log("Exporting graph.")

    output_components(output_dir)
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
                component_nodes = Node.topological_sort(component_nodes)

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
                    if node not in paths_by_start: continue
                    paths = paths_by_start[node]
                    for path in paths_by_start[node]:
                        path = [str(n.hash) for n in path]
                        pathfile.write("\t".join(path) + "\n")

                edgefile.write("InID\tOutID\tWeight\tCopycount\tNormalization\n")
                for edge in component_edges:
                    #Update edge copy count based on exact copy count
                    edge.copy_count = Read.known_edges[tuple([edge.in_node,edge.out_node])]/max(Read.L - edge.weight - 1, 1)
                    if Read.known_edges[tuple([edge.in_node,edge.out_node])] >= 2: #EDGE_TH
                        edgefile.write(edge.to_string())
                component += 1

def main():
    sys.setrecursionlimit(10000)
    if len(sys.argv) == 1:
        arguments = ['', '-f', '-e', '-d',
                     'input/kmer.dict', 'input/k1mer.dict', '--kmer=24',
                     'input/r1_100k.fa', 'input/r2_100k.fa', 'output']
    else:
        arguments = sys.argv
    names = []
    double_stranded = False
    error_correction = False
    compute_fringes = False
    cpp = False
    weighted = False
    Read.K = 24
    use_only_k1mer = False
    for arg in arguments[1:]:
        if arg == '-d':
            double_stranded = True
        elif arg == '-e':
            error_correction = True
        elif arg == '-f':
            compute_fringes = True
        elif arg == '-c':
            cpp = True
        elif arg == '-w':
            weighted = True
        elif arg[:7] == "--kmer=":
            Read.K = int(arg[7:])
        elif arg == '--only_k1':
            use_only_k1mer = True
        else:
            names.append(arg)

    node_file, edge_file, output_dir = names[0], names[1], names[-1]
    read_files = names[2:-1]

    setup(read_files[0])
    log("Starting Multibridging.")
    if len(read_files) == 1:
        load_reads(read_files[0], double_stranded, weighted)
    elif len(read_files) == 2:
        load_mated_reads(read_files[0], read_files[1], double_stranded)
    else:
        load_reads(read_files[0], double_stranded)
        load_mated_reads(read_files[1], read_files[2], double_stranded)

    if cpp:
        load_cpp(node_file, edge_file)
    else:
        if use_only_k1mer:
            load_single_jellyfish(edge_file)
        else:
            load_jellyfish(node_file, edge_file)
        
 
    run(output_dir, error_correction, compute_fringes)

if __name__ == '__main__':
    main()
