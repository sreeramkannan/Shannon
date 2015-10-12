import time, doctest

CYCLE_DESTROY = True

class Read(object):
    """A read is a segment of DNA of length L that has been read.

    It contains information about the base pairs that were read (seq) as well
    as the nodes that were constructed from the read (nodes).
    """

    #Dictionary of all reads as bases:Read object
    reads = {}
    known_paths = set()
    mate_paths = set()

    MATE_PAIR_LENGTH = 500
    MATE_PAIR_MIN_LENGTH = 0
    MATED_READS = False

    def __init__(self, bases, copy_count):
        assert bases not in Read.reads
        
        self.bases = bases
        Read.reads[bases] = self
        self.copy_count = copy_count

        #If this read is in a mate pair,
        #the other read in the pair.
        self.mate = None
        #1 if this read is the first of a mate pair; 2 if it is the second;
        #otherwise None.
        self.mate_pair = None
        #The nodes that this read traverses.
        self.nodes = None

        #The base from which this read starts, on the first node traversed.
        self.start_base = None
        #The base on which this read ends, on the last node traversed.
        self.end_base = None

    @staticmethod
    def add_read(read, double = False):
        """Given a read as a weight/string tuple, construct and return a read.

        If DOUBLE, also add the reverse complement of this read.
        """
        copy_count, bases = read
        #If already in reads, only change copy count.
        if bases in Read.reads:
            read = Read.reads[bases]
            read.copy_count += copy_count
        else:
            read = Read(bases, copy_count)

        if double:
            rc = (copy_count, Read.reverse_complement(bases))
            return (read, Read.add_read(rc, False))

        return read

    @staticmethod
    def reverse_complement(bases):
        """Return the reverse complement of BASES. Assumes BASES is
        all uppercase.

        >>> Read.reverse_complement("ATCGGGG")
        'CCCCGAT'
        """
        replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
        for ch1, ch2 in replacements:
            bases = bases.replace(ch1, ch2)
        return bases[::-1].upper()

    def bridges(self, node, index):
        """Return True if your bases bridge NODE at INDEX, False otherwise.
        """
        if index <= 0:
            return False
        if len(self.bases) <= index + len(node.bases):
            return False

        return (self.bases[index:index+len(node.bases)] == node.bases)

    @staticmethod
    def find_bridging_reads():
        """Find all reads which bridge an X-node and link them to their
        nodes.
        """
        #Construct set of first K bases of all X-nodes
        start_sequences = {}
        for node in Node.all_xnodes():
            bases = node.bases[:Read.K]
            if bases not in start_sequences:
                start_sequences[bases] = []
            start_sequences[bases].append(node)

        #Look for reads which contain those bases
        for read_bases, read in Read.reads.items():
            for start in range(1, len(read.bases)-Read.K):
                node_bases = read_bases[start:start+Read.K]
                #If read contains a start sequence and bridges the
                #corresponding X-node(s), link them.
                if node_bases in start_sequences:
                    xnodes = start_sequences[node_bases]
                    #Check for matching X-nodes
                    matches = (x for x in xnodes if read.bridges(x, start))
                    for xnode in matches:
                        xnode.link_read(read, start)

    @staticmethod
    def find_mate_pairs():
        """Find any certain paths in the graph and add them to the list of
        known paths.

        A path is certain if, for the two reads in a mate pair, there
        is one unique path of length <= MATE_PAIR_LENGTH
        connecting the start base on the first node of the first read to the
        end base on the last node of the second read.
        """
        #Set of all read node-pairs (n1, n2)
        mate_nodes = list()

        for b, r in Read.reads.items():
            node_pair = r.mate_pair_check()
            if node_pair is not None:
                mate_nodes.append(node_pair)

        print 'Total pairs: {}'.format(len(mate_nodes))
        mate_nodes = set(mate_nodes)
        print 'Unique pairs: {}'.format(len(mate_nodes))
        print 'Unique sources: {}'.format(len(set(n1 for n1, n2 in mate_nodes)))
        for node1, node2 in mate_nodes:
            paths = node1.find_mate_path(len(node1.bases) - 1, node2, 0)
            if len(paths) == 1 and len(paths[0]) > 2:
                Read.known_paths.add(tuple(paths[0]))
                Read.mate_paths.add(tuple(paths[0]))

    def mate_pair_check(self):
        """Iff this read is the first of a mate pair, and both reads have
        known paths, return the first node/start base of this read, and the
        last node/end base of the second read. Otherwise return null.
        """
        if self.mate_pair == 1 and self.nodes and self.mate.nodes:
            link = (self.nodes[0], self.mate.nodes[-1])
            if link[0] is link[1]: return
            if link[1] in link[0].successors(): return
            return link

    def __str__(self):
        return self.bases

class Edge(object):
    """An edge between nodes.
    """

    def __init__(self, weight, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.weight = weight
        self.in_node.out_edges.append(self)
        self.out_node.in_edges.append(self)
        self.traversed = False
        self.copy_count = 0.0

    def destroy(self):
        self.in_node.out_edges.remove(self)
        self.out_node.in_edges.remove(self)
        self.in_node = None
        self.out_node = None

    def condense(self, pop):
        """Given an unambiguous edge (in_node has only
        this outgoing edge, a nd out_node has only this incoming edge), condense the
        two nodes into one node, inheriting edges and read parentage.

        This method also works on self-loops but will disrupt the graph. For example, if
        there are edges (A => A, A => B), A can be condensed to itself, resulting in
        AA => B.
        """
        source = self.in_node
        destination = self.out_node
        new_bases = source.bases + destination.bases[self.weight:]

        condensed = Node(new_bases)

        if source is not destination:
            condensed.count = source.count + destination.count
            condensed.prevalence = source.prevalence + destination.prevalence
        else:
            condensed.count = source.count
            condensed.prevalence = source.prevalence

        condensed.norm = source.norm + destination.norm
        if condensed.norm == 0:
            condensed.copy_count = source.copy_count + destination.copy_count
        else:
            condensed.copy_count = ((source.copy_count * source.norm + 
                destination.copy_count * destination.norm) /
                condensed.norm)

        if source is destination:
            self.destroy()
            for edge in list(source.out_edges):
                new_e = condensed.link_to(edge.out_node, edge.weight)
                new_e.copy_count = edge.copy_count
                edge.destroy()
            for edge in list(source.in_edges):
                new_e = condensed.link_from(edge.in_node, edge.weight)
                new_e.copy_count = edge.copy_count
                edge.destroy()

            condensed.copy_count = source.copy_count / 2.0
            condensed.norm = source.norm

            condensed.reads = source.reads
            source.reads = []
            source.destroy(pop)
            return condensed

        #Update edges
        for edge in list(source.in_edges):
            previous_node = edge.in_node
            edge.in_node.link_to(condensed, edge.weight)
            edge.destroy()
        for edge in list(destination.out_edges):
            edge.out_node.link_from(condensed, edge.weight)
            edge.destroy()

        #Update reads
        source_reads = set(source.reads)
        destination_reads = [(r, i-(len(source.bases)-self.weight))
            for (r, i)
            in destination.reads]
        destination_reads = [(r, i) for (r, i) in destination_reads
            if (r, i) not in source_reads]
        condensed.reads = list(source_reads) + destination_reads
        source.reads = []
        destination.reads = []

        #Clean up database
        self.destroy()
        source.destroy(pop)
        destination.destroy(pop)
        return condensed

    def local_condense(self):
        """Given that SELF may or may not have already been destroyed, try to
        condense SELF if valid. If condensing succeeds, try to condense
        all in- and out-edges of the resulting node.
        """
        if self.in_node is None:
            return
        if len(self.in_node.out_edges) > 1 or len(self.out_node.in_edges) > 1:
            return

        condensed = self.condense(False)
        for edge in condensed.in_edges + condensed.out_edges:
            edge.local_condense()

    def extend_forward(self):
        """If SELF is an edge v => q, construct a node w as v extended forward
        one base along q. Then construct w => q and return w.
        """
        v, q = self.in_node, self.out_node
        extension = q.bases[self.weight]
        w = Node(v.bases + extension)
        w.link_to(q, self.weight+1)
        for (read, index) in v.reads:
            if read.bridges(w, index+1):
                w.link_read(read, index+1)
        w.bridged = False
        return w

    def extend_back(self):
        """If SELF is an edge p => v, construct a node u as v extended back
        one base along p. Then construct p => u and return u.
        """
        p, v = self.in_node, self.out_node
        extension = p.bases[-self.weight-1]
        u = Node(extension + v.bases)
        p.link_to(u, self.weight+1)
        for (read, index) in v.reads:
            if read.bridges(u, index-1):
                u.link_read(read, index-1)
        u.bridged = False
        return u

    def to_string(self):
        return (str(self.in_node.hash) + "\t"
            + str(self.out_node.hash) + "\t"
            + str(self.weight) + "\t"
            + str(self.copy_count) + "\t"
            + str(max(Read.L - self.weight - 1, 0)) + "\n")

    def __str__(self):
        return self.in_node.bases[:500] + '=>' + str(self.weight) + self.out_node.bases[:500]

class Node(object):
    """A Node is a sequence of base pairs that has been pulled from a read.

    It contains information about the sequence itself (bases), any reads that
    bridge it if it is an X-node, and the connecting edges in the graph
    (in_edges and out_edges).
    """

    #List of nodes
    nodes = []
    #Dictionary of K-mers at the beginning
    kmer_dict = {}
    #Average prevalence required for a K-mer to not be destroyed
    PREVALENCE_THRESHOLD = 2
    #Hamming distance required for two nodes to be collapsed
    HAMMING_FRACTION = 0.1

    def __init__(self, bases):
        self.bases = bases
        self.reads = []
        self.in_edges = []
        self.out_edges = []
        self.norm = 1.0
        self.copy_count = 0.0
        #The number of times a K-mer in this node appeared in a read
        self.prevalence = 0.0
        #The number of K-mers condensed together to form this node
        self.count = 1.0

        #Add to database
        Node.nodes.append(self)

    def in_degree(self):
        """Return SELF's in-degree.
        """
        return len(self.in_edges)
    def out_degree(self):
        """Return SELF's out-degree.
        """
        return len(self.out_edges)

    def predecessors(self):
        """Return a generator for SELF's predecessors.
        """
        return (e.in_node for e in self.in_edges)
    def successors(self):
        """Return a generator for SELF's successors.
        """
        return (e.out_node for e in self.out_edges)

    def precedes(self, node):
        """Return True if SELF is a predecessor of NODE, false otherwise.
        """
        for edge in self.out_edges:
            if edge.out_node is node:
                return True
        return False
    def succeeds(self, node):
        """Return True if NODE is a predecessor of SELF, false otherwise.
        """
        return node.precedes(self)

    def link_to(self, next_node, weight):
        """Make a link from SELF to NEXT_NODE. Return the resulting edge.
        """
        edge = Edge(weight, self, next_node)
        return edge
    def link_from(self, previous_node, weight):
        """Make a link from PREVIOUS_NODE to SELF. Return the resulting edge.
        """
        return previous_node.link_to(self, weight)

    def max_in(self):
        """Return the maximum weight of any in-edge.
        """
        return max([e.weight for e in self.in_edges]+[0])
    def max_out(self):
        """Return the maximum weight of any out-edge.
        """
        return max([e.weight for e in self.out_edges]+[0])
        
    @staticmethod
    def add_node(bases, previous_node):
        """Add a K-mer from a read with bases BASES and linked
        to previous node PREVIOUS_NODE to the database, or add the link to an
        existing node if necessary. Return the created/existing node.
        """
        #Make new node/use existing one
        if bases in Node.kmer_dict:
            n = Node.kmer_dict[bases]
        else:
            if bases in Node.kmer_dict:
                n = Node.kmer_dict[bases]
            else:
                n = Node(bases)
                Node.kmer_dict[bases] = n

        n.prevalence += 1

        #Update previous_node links
        if previous_node and not previous_node.precedes(n):
            n.link_from(previous_node, Read.K-1)

        return n

    def link_read(self, read, index):
        """Given a read READ that bridges this X-node (contains it from INDEX
        to INDEX+self.bases-1 inclusive), link READ to this X-node.
        """
        self.reads.append((read, index))
    @staticmethod
    def remove_destroyed():
        Node.nodes = [n for n in Node.nodes if not hasattr(n, 'destroyed')]
    def destroy(self, pop):
        """Remove SELF from list of nodes.

        Iff POP is false, only mark self for destruction.
        """
        assert len(self.in_edges) == 0, "Remaining links: %s" % self
        assert len(self.out_edges) == 0, "Remaining links: %s" % self
        assert len(self.reads) == 0, "Remaining reads: %s" % self.reads
        if pop:
            Node.nodes.remove(self)
        else:
            self.destroyed = True
    def full_destroy(self):
        """Remove all edges and reads and remove SELF from list of nodes.

        Iff POP is false, only mark self for destruction.
        """
        for edge in list(self.in_edges):
            edge.destroy()
        for edge in list(self.out_edges):
            edge.destroy()
        self.reads = []
        self.destroy(False)

    def is_xnode(self):
        return len(self.in_edges) >= 2 and len(self.out_edges) >= 2
    def refresh_bridging_reads(self):
        """Check whether bridging reads are still bridging reads, and
        remove if not. Reads can be invalid if they do not fully cover
        the node + 1 base at each edge, or if the bases covered at the
        edge do not correspond to nodes in the graph (e.g. if the
        original nodes have been removed).
        """
        reads = [(read, index) for (read, index) in self.reads
            if index > 0 and
            len(read.bases) > index+len(self.bases) and
            read.bases[index:index + len(self.bases)] == self.bases]
        reads = list(set(reads))
        real_reads = []

        for read, index in reads:
            bridge_in, bridge_out = False, False
            for edge in self.in_edges:
                if read.bases[index-1] == \
                    edge.in_node.bases[len(edge.in_node.bases) - edge.weight - 1]:
                    bridge_in = True
            for edge in self.out_edges:
                if read.bases[index + len(self.bases)] == \
                    edge.out_node.bases[edge.weight]:
                    bridge_out = True
            if bridge_in and bridge_out:
                real_reads.append((read, index))
        self.reads = real_reads

    @staticmethod
    def condense_all():
        """Condense edges until no more unambiguous edges remain.

        For efficiency reasons, nodes are marked as destroyed while in
        the condensing process, then all destroyed nodes are removed
        at the end.
        """
        condensed = 0
        #Find nodes with only one outgoing edge
        for source_node in (n for n in Node.nodes if len(n.out_edges) == 1):
            assert not hasattr(source_node, 'destroyed')
            #If next node also has only one incoming edge, collapse
            out_edge = source_node.out_edges[0]
            if len(out_edge.out_node.in_edges) == 1 \
                and source_node is not out_edge.out_node:
                out_edge.condense(False)
                condensed += 1
                if condensed % 1000 == 0:
                    pass
        Node.remove_destroyed()

    @staticmethod
    def bridged_xnodes():
        """Return a generator for all bridged X-nodes.
        """
        for n in Node.all_xnodes():
            if n.is_bridged_xnode():
                yield n
    def is_bridged_xnode(self):
        """Return True iff SELF is a bridged X-node.
        """
        self.refresh_bridging_reads()
        #Find if all in/out-edges are bridged
        in_bases, out_bases = set(), set()
        for read, index in self.reads:
            in_bases.add(read.bases[index-1])
            out_bases.add(read.bases[index+len(self.bases)])
        bridged_in = len(self.in_edges) - len(in_bases)
        bridged_out = len(self.out_edges) - len(out_bases)
        if bridged_in == 0 and bridged_out == 0:
            return True
        elif bridged_in == 1 and bridged_out == 1:
            return True
        return False

    @staticmethod
    def all_xnodes():
        """Return a generator for all existing X-nodes.
        """
        return (n for n in Node.nodes if n.is_xnode())
    @staticmethod
    def bridge_all():
        """Continue bridging until no bridged X-nodes remain.
        """
        while True:
            bridged = 0
            to_bridge = list(Node.bridged_xnodes())
            #print("{} to bridge".format(len(to_bridge)))
            for node in to_bridge:
                node.bridging_step()
                bridged += 1
            print("%s: Bridged %s nodes" % (time.asctime(), bridged))
            Node.remove_destroyed()
            if bridged == 0:
                return

    def bridging_step(node):
        """Perform the bridging step for a single node.
        """
        node.refresh_bridging_reads()
        assert len(node.reads) > 0
        assert len(node.in_edges) >= 2 and len(node.out_edges) >= 2
        #print('  {}: Bridging length {} with {} bridges'.format(time.asctime(), len(node.bases), len(node.reads)))

        u_list, w_list = [], []
        v_back, v_forward, self_loop_weight = None, None, None
        in_edges, out_edges = list(node.in_edges), list(node.out_edges)
        for in_edge in in_edges:
            u = in_edge.extend_back()
            if in_edge.in_node is node:
                v_back = u
                self_loop_weight = in_edge.weight
            u_list.append(u)
        for out_edge in out_edges:
            w = out_edge.extend_forward()
            if out_edge.out_node is node:
                v_forward = w
            w_list.append(w)

        for edge in list(node.in_edges):
            edge.destroy()
        for edge in list(node.out_edges):
            edge.destroy()

        if v_back:
            assert v_forward is not None
            v_forward.link_to(v_back, self_loop_weight+2)

        links = {}
        for n in u_list + w_list:
            links[n] = 0

        for read, index in list(node.reads):
            bases = read.bases[index-1:index+len(node.bases)]
            matched_u = [u for u in u_list if u.bases == bases]
            bases = read.bases[index:index+len(node.bases)+1]
            matched_w = [w for w in w_list if w.bases == bases]

            if len(matched_u) != 1 or len(matched_w) != 1:
                continue

            u, w = matched_u[0], matched_w[0]
            u.link_read(read, index - 1)
            w.link_read(read, index)
            if not u.precedes(w):
                u.link_to(w, len(node.bases))
                u.bridged = True
                w.bridged = True
                links[u] += 1
                links[w] += 1
        node.reads = []

        unbridged_u = [u for u in u_list if not u.bridged]
        unbridged_w = [w for w in w_list if not w.bridged]
        if len(unbridged_u) == 1 and len(unbridged_w) == 1:
            u, w = unbridged_u[0], unbridged_w[0]
            u.link_to(w, len(node.bases))
            links[u] += 1
            links[w] += 1
        else:
            assert len(unbridged_u) + len(unbridged_w) == 0

        link_count = sum([c for n, c in links.items()])
        for n in u_list + w_list:
            n.prevalence = (float(links[n]) / link_count) * node.prevalence

        #Condense edges if necessary
        for n in u_list + w_list:
            for edge in n.in_edges + n.out_edges:
                edge.local_condense()

        #Destroy node
        node.destroy(False)
    def connected(self):
        return len(self.in_edges) > 0 and len(self.out_edges) > 0

    @staticmethod
    def sequence():
        """Return the DNA base sequence obtained by traversing an Eulerian
        cycle of the graph.
        """
        start_node, edges = Node.eulerian_cycle()
        sequence = start_node.bases
        for edge in edges:
            sequence += edge.out_node.bases[edge.weight:]
        return sequence
    @staticmethod
    def eulerian_cycle():
        """Return an Eulerian circuit of the node graph as a list of edges.
        """
        for node in Node.nodes:
            assert len(node.in_edges) == len(node.out_edges)

        start_node = Node.nodes[0] #Start of Eulerian cycle
        edges = []

        def path_nodes():
            yield start_node
            for edge in edges:
                yield edge.out_node

        while True:
            #Look for untraversed edge
            splice_point = None
            for index, node in enumerate(path_nodes()):
                untraversed = [e for e in node.out_edges \
                               if not e.traversed]
                if len(untraversed) > 0:
                    splice_point = index
                    splice_node = node
                    next_edge = untraversed[0]
                    next_edge.traversed = True
                    splice_edges = [next_edge]
                    current_node = next_edge.out_node
                    break

            #All edges traversed
            if splice_point is None:
                break

            #Randomly traverse until arriving back at start
            while current_node is not splice_node:
                untraversed = [e for e in current_node.out_edges \
                               if not e.traversed]
                assert len(untraversed) > 0
                next_edge = untraversed[0]
                next_edge.traversed = True
                splice_edges.append(next_edge)
                current_node = next_edge.out_node

            #Splice your path in
            edges = edges[:splice_point] + splice_edges + edges[splice_point:]

        return start_node, edges

    def add_component(self):
        """Return two lists of nodes/edges in the current connected component.
        """
        nodes, edges = set(), set()
        queue = [self]

        while len(queue) > 0:
            n = queue.pop()
            if n in nodes: continue
            nodes.add(n)

            for e in n.out_edges:
                edges.add(e)

            for e in n.out_edges:
                queue.append(e.out_node)
            for e in n.in_edges:
                queue.append(e.in_node)
        return nodes, edges
    @staticmethod
    def topological_sort(nodes):
        """Topologically sort the set NODES and return the resulting list.
        """
        added = set()
        added_list = []
        fringe = [n for n in nodes if len(n.in_edges) == 0]

        while len(fringe) > 0:
            v = fringe.pop()
            if v in added:
                continue

            added.add(v)
            added_list.append(v)
            successors = (e.out_node for e in v.out_edges)
            for n in successors:
                if set(n.predecessors()).issubset(added):
                    fringe.append(n)

        if len(added) != len(nodes):
            print("Could not topological sort: only sorted %s nodes" % len(added))
        return added_list

    @staticmethod
    def find_approximate_copy_counts():
        """Find approximate copy counts based on kmer count and node length.
        """
        Read.known_paths = set()
        for node in Node.nodes:
            node.norm = len(node.bases) - Read.K + 1
            node.copy_count = float(node.prevalence) / node.norm

        for node in Node.nodes:
            for edge in node.out_edges:
                norm = max(Read.L - edge.weight - 1, 0)
                in_n, out_n = edge.in_node, edge.out_node
                nodes_count = in_n.copy_count * in_n.norm
                nodes_count += out_n.copy_count * out_n.norm
                if norm == 0:
                    edge.copy_count = 0
                else:
                    edge.copy_count = 0.5 * nodes_count / norm

    @staticmethod
    def set_pre_post_length(compute_fringes):
        """Set the pre- and post-length of each node. If COMPUTE_FRINGES,
        try to compute the true minimal length. Otherwise, set edges to
        0 and inner nodes to a sufficiently large number.
        """
        LOTS = Read.MATE_PAIR_LENGTH * 3 + Read.L * 3
        if not compute_fringes:
            for node in Node.nodes:
                if len(node.in_edges) == 0:
                    node.pre_len = 0
                else:
                    node.pre_len = LOTS
                if len(node.out_edges) == 0:
                    node.post_len = 0
                else:
                    node.post_len = LOTS
        else:
            Node.find_pre_post_length()
            for node in Node.nodes:
                if node.pre_len is None:
                    node.pre_len = LOTS
                if node.post_len is None:
                    node.post_len = LOTS

    @staticmethod
    def find_pre_post_length():
        """For each node in the graph, find the length of the minimal
        path from a start node to that node, and set node.pre_len to
        that value. Similarly for post_len.

        The lengths will be None for any node which is connected to a
        cycle.
        """
        for node in Node.nodes:
            node.pre_len = None
            node.post_len = None

        fringe = [n for n in Node.nodes if len(n.in_edges) == 0]
        while len(fringe) > 0:
            node = fringe.pop()
            if node.pre_len is not None:
                continue

            pred_lens = [e.in_node.pre_len + len(e.in_node.bases) - e.weight
                for e in node.in_edges]
            node.pre_len = max([0] + pred_lens)

            successors = (e.out_node for e in node.out_edges)
            for s in successors:
                dependencies = list(e.in_node for e in s.in_edges)
                if len([n for n in dependencies if n.pre_len is None]) == 0:
                    fringe.append(s)

        fringe = [n for n in Node.nodes if len(n.out_edges) == 0]
        while len(fringe) > 0:
            node = fringe.pop()
            if node.post_len is not None:
                continue

            succ_lens = [e.out_node.post_len + len(e.out_node.bases) - e.weight
                for e in node.out_edges]
            node.post_len = max([0] + succ_lens)

            predecessors = (e.in_node for e in node.in_edges)
            for p in predecessors:
                dependencies = list(e.out_node for e in p.out_edges)
                if len([n for n in dependencies if n.post_len is None]) == 0:
                    fringe.append(p)

    def find_mate_path(self, start_base, goal, end_base):
        """Return paths from base START_BASE on SELF to base END_BASE
        on GOAL of length less than Read.MATE_PAIR_LENGTH.

        A path consists of a list of nodes from SELF to GOAL inclusive. It includes
        both the START_BASE and END_BASE bases.

        Paths which consist of one node are ignored.
        """
        paths = []
        fringe_length = end_base + (len(self.bases) - start_base)
        min_length = Read.MATE_PAIR_MIN_LENGTH - fringe_length
        max_length = Read.MATE_PAIR_LENGTH - fringe_length
        for e in self.out_edges:
            paths += ([self] + path for path in
                e.out_node.mate_search(goal,
                    max_length + e.weight, min_length + e.weight))
        return paths

    def mate_search(self, goal, max_length, min_length):
        """Return paths from SELF to GOAL of length less
        than MAX_LENGTH and at least MIN_LENGTH.

        A path consists of a list of nodes from SELF to GOAL inclusive.
        The length of a path is the length of the base sequence
        containing SELF in its entirety, and traversing to one base in
        GOAL.
        """
        if max_length <= 0:
            return []
        if self is goal and min_length <= 1:
            return [[goal]]
        paths = []
        for edge in self.out_edges:
            new_length = len(self.bases) - edge.weight
            e_paths = edge.out_node.mate_search(goal,
                max_length - new_length, min_length - new_length)
            paths += [[self] + list(p) for p in e_paths]
        return paths

    def all_sequences_from(self, offset, traversed):
        """Return all sequences of nodes starting from this
        node excluding the first OFFSET bases. Do not traverse any node 3
        times, given that TRAVERSED is a dictionary of how many times each
        node has already been traversed.
        """
        if self not in traversed:
            traversed[self] = 0

        if traversed[self] == 2:
            return []
        traversed = dict(traversed)
        traversed[self] += 1

        if len(self.out_edges) == 0:
            return [self.bases[offset:]]
        seqs = []
        for e in self.out_edges:
            seqs += e.out_node.all_sequences_from(e.weight, traversed)
        return [self.bases[offset:] + seq for seq in seqs]

    def reachable_cycle(self, no_cycles, traversed = []):
        """Return a cycle that contains this node, or None if no
        such cycle exists.

        Given that the nodes in TRAVERSED form a path to this node.
        And that the nodes in NO_CYCLES are guaranteed to not be in
        any cycles.

        Also update NO_CYCLES.

        The cycle should be a list of nodes, where the first and last
        nodes are the same.
        """
        traversed = list(traversed)
        traversed.append(self)
        for node in (e.out_node for e in self.out_edges):
            if node in traversed:
                contains_cycle = traversed + [node]
                return contains_cycle[contains_cycle.index(node):]
            if node in no_cycles:
                continue

            possible_cycle = node.reachable_cycle(no_cycles, traversed)
            if possible_cycle:
                return possible_cycle
        no_cycles.add(self)

    def fast_dfs(self, start, traversed, limit):
        """Return a cycle that contains START, or None if no
        such cycle exists, given that TRAVERSED nodes have already been
        traversed.

        The cycle should be a list of nodes, where the first and last
        nodes are the same.
        """
        if limit == 0: return
        traversed.add(self)
        for node in (e.out_node for e in self.out_edges):
            if node is start:
                return [self, start]
            if node in traversed: continue
            ret = node.fast_dfs(start, traversed, limit - 1)
            if ret is not None:
                return [self] + ret
        return

    @staticmethod
    def find_cycle(no_cycles):
        """Return a cycle in this graph, or None if one does not exist.

        The cycle should be a list of nodes, where the first and last
        nodes are the same.
        """
        for node in Node.nodes:
            if node not in no_cycles:
                possible_cycle = node.reachable_cycle(no_cycles)
                if possible_cycle:
                    return possible_cycle

    def has_self_loop(self):
        """If SELF has a self-loop, return True. Otherwise
        return False.
        """
        self_edges = [e for e in self.out_edges \
            if e.out_node is self]
        if len(self_edges) > 0:
            return True
        else:
            return False

    def break_self_loop(self):
        """Assumes that SELF has a self-loop. Break said loop.
        """
        print("Breaking self loop %s" % self.bases)
        self_edge = [e for e in self.out_edges \
            if e.out_node is self][0]
        self_edge.condense(False)

    @staticmethod
    def break_long_cycle(cycle):
        """Break a cycle in this graph which is not a self-loop.
        If the cycle contains a node which has a self-loop, break
        that instead.
        """
        xnodes = [n for n in cycle if len(n.in_edges) >= 2
            and len(n.out_edges) >= 2]
        if len(xnodes) == 0:
            xnodes = [n for n in cycle if len(n.in_edges) >= 2
                or len(n.out_edges) >= 2]
        if len(xnodes) == 0: return

        min_disrupt = None
        best_split_node, best_c_in, best_c_out, best_o_in, best_o_out = \
            None, None, None, None, None

        for node in xnodes:
            if node.has_self_loop():
                node.break_self_loop()
                return
            prev_node = cycle[cycle.index(node) - 1]
            next_node = cycle[(cycle.index(node) + 1) % len(cycle)]
            c_in = [e for e in node.in_edges if e.in_node is prev_node]
            c_out = [e for e in node.out_edges if e.out_node is next_node]
            o_in = [e for e in node.in_edges if e.in_node is not prev_node]
            o_out = [e for e in node.out_edges if e.out_node is not next_node]

            disrupt = abs(sum(e.copy_count for e in o_in) - c_out[0].copy_count) + \
                abs(sum(e.copy_count for e in o_out) - c_in[0].copy_count)
            if min_disrupt is None or disrupt < min_disrupt:
                min_disrupt = disrupt
                best_split_node = node
                best_c_in, best_c_out, best_o_in, best_o_out = \
                    c_in, c_out, o_in, o_out

        node = best_split_node
        print("Breaking cycle on node %s" % node.bases)
        clone = Node(node.bases)
        clone.norm = 0
        node.norm = 0

        edge_cc = sum(e.copy_count for e in best_c_in + best_c_out + best_o_in + best_o_out)
        if edge_cc == 0:
            proportion = 0.5
        else:
            proportion = float(sum(e.copy_count for e in best_c_in + best_o_out))
            proportion /= edge_cc
        clone.copy_count = node.copy_count * proportion
        node.copy_count = node.copy_count * (1 - proportion)

        for e in best_c_in:
            clone.link_from(e.in_node, e.weight)
            e.destroy()
        for e in best_o_out:
            clone.link_to(e.out_node, e.weight)
            e.destroy()

    @staticmethod
    def break_cycle(cycle):
        """Break a cycle in this graph, minimizing the disruption to
        flow. Cycle is denoted by CYCLE, a list of nodes, with the first
        node and last node the same.
        """
        cycle = cycle[1:]

        if CYCLE_DESTROY:
            node = cycle[0]
            node.full_destroy()
        elif len(cycle) == 1:
            node = cycle[0]
            node.break_self_loop()
        else:
            Node.break_long_cycle(cycle)

    @staticmethod
    def aggressive_break_cycles():
        import heapq

        N = len(Node.nodes)
        print("Topologically sorting %s nodes." % N)
        processed = set()
        delta_heap = []
        sources, sinks = [], []
        next_sources, next_sinks = [], []
        in_count, out_count = [0] * N, [0] * N

        for i in range(N):
            node = Node.nodes[i]
            node.id = i
            in_count[i] = node.in_degree()
            out_count[i] = node.out_degree()
            if in_count[i] == 0:
                next_sources.append(i)
            elif out_count[i] == 0:
                next_sinks.append(i)
            else:
                heapq.heappush(delta_heap, (in_count[i] - out_count[i], i))

        def add_source(src):
            processed.add(src)
            sources.append(src)
            for node in (n.id for n in Node.nodes[src].successors()):
                in_count[node] -= 1
                if in_count[node] == 0:
                    next_sources.append(node)
                else:
                    heapq.heappush(delta_heap, (in_count[node] - out_count[node], node))

        def add_sink(sink):
            processed.add(sink)
            sinks.append(sink)
            for node in (n.id for n in Node.nodes[sink].predecessors()):
                out_count[node] -= 1
                if out_count[node] == 0:
                    next_sinks.append(node)
                else:
                    heapq.heappush(delta_heap, (in_count[node] - out_count[node], node))

        while len(delta_heap) > 0:
            if len(next_sources) > 0:
                s = next_sources.pop()
                if s not in processed:
                    add_source(s)
            elif len(next_sinks) > 0:
                s = next_sinks.pop()
                if s not in processed:
                    add_sink(s)
            else:
                delta, n = heapq.heappop(delta_heap)
                if n in processed: continue
                if delta != in_count[n] - out_count[n]: continue
                add_source(n)
        print("%s processed." % len(processed))

        nodes = sources + sinks[::-1]
        edges_removed = 0
        for i in range(N):
            Node.nodes[nodes[i]].id = i

        for node in Node.nodes:
            for out_edge in list(node.out_edges):
                if out_edge.out_node.id <= node.id:
                    out_edge.destroy()
                    edges_removed += 1
        print("%s edges removed." % edges_removed)

        for node in Node.nodes:
            for out_edge in list(node.out_edges):
                assert out_edge.out_node.id > node.id, "Whoa, problem."

    @staticmethod
    def break_cycles(dfs=True):
        """Break cycles in this graph, minimizing the disruption to
        flow.
        """
        if dfs:
            limit = 1
            while limit < 100:
                for i in range(len(Node.nodes)):
                    n = Node.nodes[i]
                    if (hasattr(n, 'destroyed')): continue
                    cycle = n.fast_dfs(n, set(), limit)
                    if cycle:
                        print("%s-cycle found at %s: %s" % (len(cycle) - 1, i, time.time()))
                        Node.break_cycle(cycle)
                Node.remove_destroyed()
                Node.condense_all()
                print("Found all cycles with limit %s, %s nodes remain: %s" % (limit, len(Node.nodes), time.time()))
                limit += 1
        
        no_cycles = set()

        cycle = Node.find_cycle(no_cycles)
        while cycle is not None:
            print("Found cycle of length %s at %s" % (len(cycle), time.time()))
            Node.break_cycle(cycle)
            cycle = Node.find_cycle(no_cycles)
        Node.remove_destroyed()
        Node.condense_all()
        assert Node.find_cycle(set()) is None

    def average_prevalence(self):
        """Return the average prevalence of K-mers in this node.
        """
        return self.prevalence / self.count

    @staticmethod
    def destroy_suspicious():
        """After the condensing stage, destroy all suspicious nodes.
        """
        while Node.destroy_some_suspicious(): pass

    @staticmethod
    def check_suspicious():
        """Report statistics on the distribution of suspicious nodes.
        """
        Node.nodes.sort(key = lambda n: n.average_prevalence())
        buckets = 100
        for b in range(buckets):
            start = int((float(b)/buckets) * len(Node.nodes))
            print("Bucket %d: %.2f" % (b, Node.nodes[start].average_prevalence()))
        print("Max: %.2f" % Node.nodes[-1].average_prevalence())

    def is_suspicious(self):
        """Return True iff SELF is a suspicious node.

        Suspicious nodes meet the following four conditions:
        1. The average prevalence is less than Node.PREVALENCE_THRESHOLD.
        2. The node size is less than Node.SIZE_THRESHOLD.
        3. This node has predecessors and their average out-degree
        is at least 2.
        4. This node has successors and their average in-degree
        is at least 2.
        """
        if len(self.bases) <= Node.SIZE_THRESHOLD and (self.in_degree() == 0 or self.out_degree() == 0):
            return True

        if self.average_prevalence() >= Node.PREVALENCE_THRESHOLD:
            return False
        #if len(self.bases) >= Node.SIZE_THRESHOLD:
        #    return False

        if self.in_degree() == 0 or self.out_degree() == 0:
            return True

        total_degree = 0
        predecessors = list(self.predecessors())
        for p in predecessors:
            total_degree += p.out_degree()
        if len(predecessors)>0:
            average_out = float(total_degree) / len(predecessors)
            if average_out < 2:
                return False

        total_degree = 0
        successors = list(self.successors())    
        for s in successors:
            total_degree += s.in_degree()
        if len(successors)>0:
            average_in = float(total_degree) / len(successors)
            if average_in < 2:
                return False

        return True

    @staticmethod
    def collapse_all():
        """Collapses all similar nodes.
        """
        while Node.collapse_some():
            pass

    @staticmethod
    def collapse_some():
        """Collapses some similar nodes. Returns True if any similar nodes
        existed.
        """
        collapsed = False
        for node in Node.nodes:
            if node.collapse_out():
                collapsed = True

        Node.remove_destroyed()
        return collapsed

    def collapse(self, other):
        """Collapse two similar nodes: one is destroyed and the other gains the
        prevalence of the smaller node.
        """
        if self.prevalence < other.prevalence:
            other.collapse(self)
            return

        self.prevalence += other.prevalence
        other.full_destroy()

    def collapse_out(self):
        """If two of this node's successors are similar enough, collapse them.
        Return True if any two nodes were collapsed.
        """
        successors = [n for n in self.successors() if not hasattr(n, 'destroyed')]
        for i in range(len(successors)):
            for j in range(i + 1, len(successors)):
                if successors[i].similar(successors[j]):
                    successors[i].collapse(successors[j])
                    return True
        return False

    def similar(self, other):
        """Returns True iff SELF's bases are of the same length as OTHER and
        the Hamming distance is less than Node.HAMMING_FRACTION, and all
        in- and out-edges are identical.
        """
        if hasattr(self, 'destroyed') or hasattr(other, 'destroyed'):
            return False

        if len(self.bases) != len(other.bases):
            return False
        mismatches = len([i for i in range(len(self.bases)) if self.bases[i] != other.bases[i]])
        if float(mismatches) / max(len(self.bases), 1) >= Node.HAMMING_FRACTION:
            return False
        if set(self.successors()) != set(other.successors()):
            return False
        if set(self.predecessors()) != set(other.predecessors()):
            return False

        return True

    @staticmethod
    def destroy_some_suspicious():
        """Destroy all suspicious nodes. This may create new
        suspicious nodes.


        Return True iff there were some suspicious nodes.
        """
        suspicious_nodes = list(n for n in Node.nodes \
            if n.is_suspicious())
        if len(suspicious_nodes) == 0:
            return False
        suspicious_nodes.sort(key = lambda n: n.average_prevalence())

        for node in (n for n in suspicious_nodes
            if not hasattr(n, 'destroyed')):
            adjacent = [e.in_node for e in node.in_edges]
            adjacent += [e.out_node for e in node.out_edges]
            node.full_destroy()
            for n in adjacent:
                n.local_condense()
        Node.remove_destroyed()
        return True

    def local_condense(self):
        """Try to condense every edge incident to SELF, if valid.
        """
        if len(self.out_edges) == 1:
            self.out_edges[0].local_condense()
        if len(self.in_edges) == 1:
            self.in_edges[0].local_condense()

    @staticmethod
    def disregard_loops():
        """Set norm/copycount of all nodes which have an edge to
        themselves to 0.
        """
        for node in Node.nodes:
            if node in [e.out_node for e in node.out_edges]:
                node.norm = 0
                node.copy_count = 0

    @staticmethod
    def write_all_sequences():
        """Write all sequences in the graph to a text file.
        """
        seqs = []
        for node in Node.nodes:
            if len(node.in_edges) == 0:
                seqs += node.all_sequences_from(0, {})
        with open('output/all_sequences.txt', 'w') as f:
            for seq in seqs:
                f.write(seq + '\n')

    def to_string(self):
        return (str(self.hash) + "\t" + self.bases + "\t"
            + str(self.copy_count) + "\t"
            + str(self.norm) + "\n")

    def __str__(self):
        in_string = "["+", ".join(sorted(e.in_node.bases[:500]+str(e.weight) for e in self.in_edges))+"]"
        out_string = "["+", ".join(sorted(e.out_node.bases[:500]+str(e.weight) for e in self.out_edges))+"]"
        return in_string + " => " + self.bases[:500] + " => " + out_string

def known_paths():
    """Generate all known paths by laying reads onto graph.
    """
    #Record K-mers for all nodes
    Read.known_paths = set()
    kmers = {}
    for node in Node.nodes:
        for i in range(len(node.bases) - Read.K + 1):
            kmer = node.bases[i:i+Read.K]
            if kmer not in kmers: kmers[kmer] = []
            kmers[kmer].append((node, i))

    for read in Read.reads:
        #if len(read) > 200: continue
        start_kmer = read[:Read.K]
        if start_kmer not in kmers: continue
        for start_node, start_i in kmers[start_kmer]:
            for path in search_sequence(read, start_node, start_i):
                Read.reads[read].nodes = path
                if len(path) > 2:
                    Read.known_paths.add(tuple(path))
def search_sequence(seq, node, i):
    """Returns a list of all node-paths starting from NODE at index I
    and matching sequence SEQ.
    """
    node_length = len(node.bases) - i
    if not compare(seq, node.bases[i:]):
        return []
    if len(seq) <= node_length:
        return [[node]]
    paths = []
    seq = seq[node_length:]
    for edge in node.out_edges:
        for path in search_sequence(seq, edge.out_node, edge.weight):
            paths.append([node] + path)
    return paths
def compare(s1, s2):
    """Return True iff s1 and s2 are the same up to their
    minimum length.
    """
    length = min(len(s1), len(s2))
    return s1[:length] == s2[:length]

def read_for_path(path):
    bases = [path[0].bases]
    for i in range(len(path) - 1):
        u, v = path[i], path[i+1]
        for edge in u.out_edges:
            if edge.out_node == v:
                bases.append(v.bases[edge.weight:])
                break
    return (0., ''.join(bases))
def construct_reads():
    for path in Read.mate_paths:
        Read.add_read(read_for_path(path), False)

def clear():
    Read.reads = {}
    Read.known_paths = set()
    Node.kmer_dict = {}
    Node.nodes = []

if __name__ == '__main__':
    doctest.testmod()
