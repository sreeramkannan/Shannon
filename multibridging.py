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

def load_reads(filename, double_stranded, weighted, no_reads_cutoff):
    """Load FASTA reads.
    """
    log("Loading reads.")
    for (curr_no,read) in enumerate(extract_reads(filename, weighted)):
        if curr_no <= no_reads_cutoff:
            Read.add_read(read, double_stranded)
        else:
            break

def load_reads_inMem(rps, double_stranded, no_reads_cutoff):
    rl = 0; nreads = 0
    for i in range(min(len(rps),no_reads_cutoff)):
        #print(rps[i])
        Read.add_read((1.0,rps[i].strip()),double_stranded)
        rl+=len(rps[i]); nreads+=1
    Read.L=int(float(rl)/nreads)
    Node.SIZE_THRESHOLD = Read.L
    print("No of reads:"+str(nreads))


def load_mated_reads_inMem(rps1,rps2, double_stranded, no_reads_cutoff):
    log("Loading mated reads.")
    assert len(rps1)==len(rps2)
    Read.MATED_READS = True
    def pair(r1, r2):
        r1.mate_pair = 1
        r2.mate_pair = 2
        r1.mate = r2
        r2.mate = r1
    nr=0; lr = 0;
    for i in range(min(len(rps1),no_reads_cutoff)):
        r1 = Read.add_read((1.0, rps1[i].strip()), double_stranded)
        r2 = Read.add_read((1.0, Read.reverse_complement(rps2[i].strip())), double_stranded)
        nr+=1; lr+= len(rps1[i])+len(rps2[i])
        if double_stranded:
                r1a, r1b = r1
                r2a, r2b = r2
                pair(r1a, r2a)
                pair(r2b, r1b)
        else:
                pair(r1, r2)
    Read.L = int(lr/(2.0*nr))
    Node.SIZE_THRESHOLD = Read.L
    print('No of paired reads:' + str(nr))

def load_mated_reads(file_1, file_2, double_stranded, no_reads_cutoff):
    #Does not allow weighted reads. 
    #Only read reads till no_reads_cutoff 
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
        for (curr_no,line1) in enumerate(gen1):
            if curr_no > no_reads_cutoff:
                break
            #Continue otherwise
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

def load_kmers_inMem(contigs,weights):
    """Loads condensed files from Jellyfish.
    """
    log("Loading nodes from K-mer files.")

    nodes = {}    
    for i in range(len(weights)):
        for (j,prevalence) in enumerate(weights[i]):
            bases=contigs[i][j:j+Read.K+1]
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
    print('No of kmers:' + str(len(Node.nodes)))


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
    Node.find_copy_counts()
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
                    #np = tuple([edge.in_node,edge.out_node]) #node-pair
                    if edge.copy_count > 0: #either the edge has a copy count or edge weight >= Read.K
                        #edge.copy_count = max(Read.known_edges.get(np,0),1)/max(Read.L - edge.weight - 1, 1)
                        edgefile.write(edge.to_string())
                component += 1  

def main(arguments,inMem=False,contigs=[],weights=[],rps=[]):
    arguments=arguments.strip().split()
    sys.setrecursionlimit(10000)
    names = []
    double_stranded = False
    error_correction = False
    compute_fringes = False
    cpp = False
    weighted = False
    Read.K = 24
    use_only_k1mer = False
    for arg in arguments:
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

    
    log("Starting Multibridging.")
    print(inMem)
    if inMem:
        load_kmers_inMem(contigs,weights);
        del contigs[:], weights[:]
        no_kmers = len(Node.nodes);
        FACTOR = 10
        no_reads_cutoff = no_kmers*FACTOR;

        if len(rps)==1:
            load_reads_inMem(rps[0],double_stranded,no_reads_cutoff)
            del rps[0][:]
        elif len(rps)==2:
            load_mated_reads_inMem(rps[0],rps[1],double_stranded,no_reads_cutoff)
            del rps[0][:], rps[1][:]
    else:    
        setup(read_files[0])
        if cpp:
            load_cpp(node_file, edge_file)
        else:
            if use_only_k1mer:
                load_single_jellyfish(edge_file)
            else:
                load_jellyfish(node_file, edge_file)
            
        no_kmers = len(Node.nodes);
        FACTOR = 10
        no_reads_cutoff = no_kmers*FACTOR;

        if len(read_files) == 1:
            #Single-end
            load_reads(read_files[0], double_stranded, weighted, no_reads_cutoff)
        elif len(read_files) == 2:
            #Paired-end
            load_mated_reads(read_files[0], read_files[1], double_stranded, no_reads_cutoff)
        else: 
            #There is an unpaired read file and 2 paired_read files
            load_reads(read_files[0], double_stranded, no_reads_cutoff)
            load_mated_reads(read_files[1], read_files[2], double_stranded, no_reads_cutoff)
    print('output dir is ' + output_dir)
    run(output_dir, error_correction, compute_fringes)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        arguments = ['', '-f', '-e', '-d',
                     'input/kmer.dict', 'input/k1mer.dict', '--kmer=24',
                     'input/r1_100k.fa', 'input/r2_100k.fa', 'output']
    else:
        arguments = sys.argv; arguments = arguments[1:]
    main('\t'.join(arguments))
