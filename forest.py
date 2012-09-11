#
# Implements a genealogy, and an evolutionary model
#

import random
import bisect
import math
import heapq
import array
import sys

max_float = 1e10
fake_random = 1

# for covariance: need to check if tmrca and pos are reported consistently for ms and own code

#
# - try removing checks on on_local_tree in coalescent loop, and rely just on identity of root,
#   to simplify code
#
# - currently, some invariants are updated halfheartedly, and reset after implementing
#   a recombination.  Investigate which we don't need to reset
#
# - recombheight not remembered
#

class FakeRandom(random.Random):
  def __init__(self):
    self.file = open("./scrm/random.numbers", "r")

  def random(self):
    return float(self.file.readline())

  def seed(self, seed):
    if seed == 0: return 0
    for i in range(seed):
      self.file.readline()

  def sampleInt(self, number):
    return int(self.random() * number)

  def sample(self, population, size):
    sample1 = self.sampleInt(len(population)-2)
    sample2 = self.sampleInt(len(population)-1)
    if sample1 == sample2: sample2 = size - 1
    return [population[sample1], population[sample2]]


##########################################################################################
#
# utilities
#
##########################################################################################


def expovariate( lmbda ):
    """ sample from an exponential distribution with rate lambda """
    return -math.log( random.random() ) / lmbda

def argmin( values ):
    """ returns the index of the minimum value of a list """
    return array.array('d',values).index( min(values) )



##########################################################################################
#
# model settings
#
##########################################################################################


class Model():
    """ geometry of data and model parameters.  also contains emissions """

    def __init__(self, 
                 samples,                   # number of samples
                 Ne=10000,                  # effective sample size; 2 Ne is the number of haploid chromosomes
                 mu=2.5e-8,                 # mutation rate (per nt per generation)
                 rho=1e-8,                  # recombination rate (per nt per generation)
                 clades=None,               # clades (not implemented)
                 maxmigrations=0,           # max number of migrations allowed (not implemented)
                 phased=False,              # whether emission data is phased (0/1 coded) or not (0/1/2 coded)
                 rooted=False,              # whether code 0 is always the ancestor
                 remove_recombination_branch=False,   # set True for SMC model
                 remove_nonancestral_twig=False,      # set True for SMC' model
                 wiufhein=False,                      # set for Wiuf and Hein model - complete all recombination events straight away
                 min_terminal_branch_length=0,        # prune short non-sample terminals
                 prune_interval = 1):      # prune every so often
        
        self.samples = samples
        self._Ne = Ne
        self.mu = mu
        self.rho = rho
        self.clades = clades
        self.maxmigrations = maxmigrations
        self.phased = phased
        self.rooted = rooted
        self.remove_recombination_branch = remove_recombination_branch
        self.remove_nonancestral_twig = remove_nonancestral_twig
        self.wiufhein = wiufhein
        self.min_terminal_branch_length = min_terminal_branch_length
        self.prune_interval = prune_interval

        # derived settings
        self.num_nodes = 2*samples + maxmigrations

    def mutationRate(self, tree):
        return self.mu * tree.treeLength()

    def setMutationData(self, data ):
        self.mutationdata = data

    def getMutationData(self):
        return self.mutationdata

    def setSequenceLength(self, length ):
        self.sequencelength = length

    def getSequenceLength(self):
        return self.sequencelength

    def getNe(self):
        """ Returns effective population size.  Will be time-dependent in future """
        return self._Ne

##########################################################################################
#
# the genealogy, simulation methods, and emission code
#
##########################################################################################


class Genealogy():
    """ A modifyable genealogy with recombination/coalescence event leading to it.
        Includes methods for simulation of recombination/coalescence and emissions, 
        and likelihood calculation of emission.
        Also holds recombination and coalescent heights of events leading to (left of) this genealogy """

    def __init__(self, model, randomstart = False):

        # geometry
        self.samples = model.samples

        # state
        self.parents = [0] * model.num_nodes        # 0 = empty, -1 = root
        self.heights = [0.0] * model.num_nodes      # can be max_float for grand MRCA root (no longer true I think)
        self._recomb_idx = 0             # branch in which recombination occurred
        self._recomb_height = 0.0        # height at which recombination occurred
        self._coalescence_height = 0.0   # height at which coalescence occurred
        
        # cache
        self.on_local_tree = [0] * model.num_nodes  # codes: 0, disjunct from sample graph; 
                                                    #        1, ancestral to samples; 
                                                    #        2, on local tree
        self.onechild = [-1] * model.num_nodes      # one of the node's children, or -1 if none
        self.siblings = [-1] * model.num_nodes      # defined if parent > 0 and parent has 2 children
        self.local_mrca = -1                        # mrca of local tree
        self.local_tree_root = -1                   # root of local tree
        self._emitpatterns = []
        self._treelength = -1
        self._forestlength = -1

        # for iterator implementation
        self._twigiter = 0
        self._previouseventidx = 0
        self._contemporaries = []

        # for tree printing
        self.tree_separation = 3   # separation between trees
        self.legend_separation = 4 # separation between trees and legend
        self.node_spacing = 1      # vertical separation between successive nodes
        self.branch_spacing = 2    # horizontal separation between branches of a node
        self._scaling_Ne = model.getNe()

        # tree pruning
        self.prune_counter = 0

        # debug
        self.do_debug = True

        self.reset(model, randomstart)
        

    def debug(self, *msg):
        if self.do_debug: print " ".join(map(str,msg))

    def debug_list(self):
        self.debug("index:   ","  ".join("%4s" % p for p in range(len(self.parents))))
        self.debug("parents: ","  ".join("%4s" % p for p in self.parents))
        self.debug("onechild:","  ".join("%4s" % p for p in self.onechild))
        self.debug("siblings:","  ".join("%4s" % p for p in self.siblings))
        self.debug("height: "," ".join("%5.f" % h for h in self.heights))

    def reset(self, model, randomstart):
        """ reset to star-like tree, or a random sample from the coalescent """

        if randomstart:
            self._reset_sample(model)
        else:
            self._reset_star(model)

        # set siblings and on_local_tree
        self.resetInvariant()

    def _reset_star(self, model):
        """ reset to a star-like tree """

        for i in range(model.num_nodes):
            self.parents[i] = 0
            self.heights[i] = 0.0

        for i in range(model.samples):
            self.parents[i] = max(i + model.samples - 1, model.samples)
            if i > 0:
                self.parents[i + model.samples - 1] = i + model.samples
                self.heights[i + model.samples - 1] = 2*model.getNe()

        self.parents[2*model.samples - 2] = -1

    def _reset_sample(self, model):
        """ sample a tree from the standard coalescent """

        n = model.samples
        new = n
        Ne = model.getNe()
        coalescing = range( n )
        for c in coalescing:
            self.heights[c] = 0.0
        time = 0
        while len(coalescing)>1:
            # compute rate
            rate = (1.0/(2*Ne))*n*(n-1)/2.0
            # get new time
            time += expovariate( rate )
            # get individuals
            i1, i2 = random.sample( coalescing, 2 )
            # implement
            self.parents[i1] = new
            self.parents[i2] = new
            self.heights[new] = time
            coalescing = [c for c in coalescing if c not in [i1,i2]] + [ new ]
            new += 1
        self.parents[ coalescing[0] ] = -1
        
    def siblings_well_paired(self):
        """ check that siblings are well paired """
        if not self.do_debug: return True
        for node, parent in enumerate(self.parents):
            # skip empties
            if parent == 0: continue
            # check siblings
            sib = self.siblings[node]
            if sib > -1:
                if self.siblings[sib] != node:
                    print "*** error *** : siblings ",node," and ",sib," not paired"
                    self.debug_list()
                    return False
        return True

    def __repr__(self):
        # return self.repr_simple()
        # return repr_ms()
        return self.repr_ascii()


    ##########################################################################################
    #
    # various representations
    #
    ##########################################################################################

    def repr_simple(self):
        s = " recombidx=%s  recombheight=%1.1f\n" % (self._recomb_idx,self._recomb_height)
        for i,p in enumerate(self.parents):
            if p != 0:
                s += "%s: %s [%1.1f]\n" % (i,p,self.heights[i])
        return s


    def repr_ms(self, scaling=0.0):
        """ local tree (only) in  ms  style """
        root = self.getLocalMRCA()
        return self._ms_writer( root, scaling ) + ';'


    def _ms_scale(self, height, scaling):
        if scaling > 0.0:
            return height / (self._scaling_Ne * scaling)
        return height


    def _ms_writer(self, node, scaling):
        # first decide if node leaf or binary internal node
        child1 = self.getChild( node )

        # deal with leaf node
        if child1 == -1:
            height = self.getHeight( node )
            if height == 0.0:
                #  contemporary leaf
                return "%s" % (node+1)
            else:
                # hanging leaf
                return "0+%1.3f" % self._ms_scale(self.getHeight( node ), scaling)

        # get other child, and sort by height
        child2 = self.siblings[child1]
        child1, child2 = max(child1,child2), min(child1,child2)

        height_above1 = self._ms_scale( self.getHeight( node ) - self.getHeight( child1 ), scaling )
        height_above2 = self._ms_scale( self.getHeight( node ) - self.getHeight( child2 ), scaling )

        return "(%s:%1.3f,%s:%1.3f)" % (self._ms_writer(child1,scaling), height_above1, 
                                        self._ms_writer(child2,scaling), height_above2)


    def repr_ascii(self, graph_output=False):
        """ ascii-graphical version """

        alloutput = []
        x_offset = 0
        max_y = 0
        for root in self.getRoots():
           x_offset, y_pos, root_pos, output = self.ascii_writer(x_offset, root, graph_output)
           max_y = max(max_y, y_pos)
           x_offset += self.tree_separation
           alloutput += output
        max_x = x_offset - self.tree_separation

        legendoutput, max_x = self.ascii_legend_writer(max_x + self.legend_separation)
        alloutput = legendoutput + alloutput

        # make canvas
        lines = [ [" " for x in range(max_x+1)] for y in range(max_y+1) ]

        # draw on canvas
        for x,y,string in alloutput:
            for idx,c in enumerate(string):
                lines[y][x+idx] = c
        
        # output
        return '\n'.join( ''.join(lines[y]) for y in range(max_y,-1,-1) )


    def ascii_height_calculator(self, node):

        duplicates = 0
        for n in range(node):
            if self.getHeight(n) == self.getHeight(n+1):
                duplicates += 1
        return self.node_spacing * max(0, node - duplicates)


    def ascii_legend_writer(self, xpos):

        maxx = xpos
        output = []
        for node in range(self.samples-1, len(self.parents)):
            if self.getParent(node) == 0: continue
            
            ypos = self.ascii_height_calculator(node)

            # write height
            heightstr = "%1.1f" % self.getHeight(node)
            maxx = max(maxx, xpos + len(heightstr) - 1)
            output.append( (xpos, ypos, heightstr) )

            # write rule
            if ypos % 5 == 0 and ypos>0:
                for x in range(1, xpos-1, 2):
                    output.append( (x, ypos, ".") )
                    
        return output, maxx
    

    def ascii_writer(self, x_offset, node, graph_output):

        # calculate y position of node
        y_pos = self.ascii_height_calculator(node)

        # first decide if node leaf or binary internal node
        child1 = self.getChild( node )

        # squash numbers if graphical output required
        gofactor = [1,0][graph_output]

        # deal with leaf node
        if child1 == -1:
            output = [(x_offset,y_pos,str(node))]
            top_right_x = x_offset + (len(output[0][2]) - 1)*gofactor
            top_right_y = y_pos
            root_pos = x_offset
            return (top_right_x, top_right_y, root_pos, output)

        # get other child, and sort by height
        child2 = self.siblings[child1]
        child1, child2 = max(child1,child2), min(child1,child2)

        # print left child
        top_right_x1, top_right_y1, root_pos1, output1 = self.ascii_writer(x_offset, child1, graph_output)

        # print vertical part of left branch
        output = []
        for y in range(top_right_y1 + 1, y_pos):
            output.append( (root_pos1, y, "|") )

        # special case: no sibling
        if child2 == -1:

            # position of node
            x_node = root_pos1
            top_right_x2 = max(top_right_x1, x_node + (len(str(node))-1)*gofactor)
            output2 = []
            # to help graph_output:
            child2, root_pos2 = child1, root_pos1

        else:

            # two children.  print right child, at an offset
            top_right_x2, top_right_y2, root_pos2, output2 = \
                self.ascii_writer(top_right_x1 + self.branch_spacing, child2, graph_output)

            # calculate position of node
            x_node = (root_pos1 + root_pos2)//2

            output.append( (root_pos1, y_pos, "/") )
            for x in range(root_pos1 + 1, x_node):
                output.append( (x, y_pos, "-") )

            # print right branch
            for y in range(top_right_y2 + 1, y_pos):
                output.append( (root_pos2, y, "|") )
            output.append( (root_pos2, y_pos, "\\") )
            for x in range(x_node + len(output[0][2]), root_pos2):
                output.append( (x, y_pos, "-") )

        # hack: generate data for graphical output
        if graph_output:
            print "GRAPH\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (root_pos1, x_node, root_pos2, 
                                                             self.getHeight(child1), self.getHeight(node), 
                                                             self.getHeight(child2), 
                                                             self.on_local_tree[child1], self.on_local_tree[child2])

        # print node
        output.append((x_node, y_pos, str(node)))

        # return coordinates and output
        return (top_right_x2, y_pos, x_node, output1 + output2 + output)

    ##########################################################################################
    #
    # simple accessors and iterators
    #
    ##########################################################################################

    def numNodes(self):
        """ number of nodes in graph.  Slow algorithm - used only for reporting """
        nn = 0
        for p in self.parents:
            if p != 0:
                nn += 1
        return nn

    def getParent(self, idx):
        assert 0 <= idx < len(self.parents)
        return self.parents[idx]


    def getChild(self, parent):
        # Fast version, uses cache
        return self.onechild[parent]
        """ simple/slow implementation, used for printing """
        for child, theparent in enumerate(self.parents):
            if theparent <= 0:
                continue
            if theparent == parent:
                return child
        return -1


    def getRoots(self):
        """ returns all roots, in indeterminate order """
        roots = []
        for node, parent in enumerate(self.parents):
            if parent == -1:
                roots.append(node)
        return roots


    def getHeight(self, idx):
        """ Returns height of node, or max_float for root """
        if idx == -1: raise ValueError("height of -1 = <parent of root> requested")
        assert 0 <= idx < len(self.parents)
        return self.heights[idx]


    def getRecombinationHeight(self):
        return self._recombheight


    def getCoalescenceHeight(self):
        return self._coalescenceheight


    def _getMRCA(self, nodes):
        nodes.sort()            # all nodes ancestral to initial ones
        heap = nodes[:]         # nodes whose ancestry needs to be checked
        while len(heap)>1:
            parent = self.getParent( heapq.heappop( heap ) )
            if parent not in nodes:
                nodes.append( parent )
                heapq.heappush( heap, parent )
        return heap[0]


    def getLocalMRCA(self, nodes=None):
        if nodes:
            return self._getMRCA(nodes)
        if self.local_mrca == -1:
            self.populateLocalTreeStatus()
        return self.local_mrca


    def getTMRCA(self, nodes = None):
        return self.getHeight( self.getLocalMRCA( nodes ) )


    def twigIterator(self):
        self._twigiter = 0


    def twigNext(self):
        """ returns indices of twigs, ordered by height """
        # skip empty nodes
        while self._twigiter < len(self.parents) and self.getParent( self._twigiter ) == 0:
            self._twigiter += 1
        if self._twigiter >= len(self.parents):
            return None
        self._twigiter += 1
        return self._twigiter - 1


    def eventIterator(self):
        self.twigIterator()
        self._previouseventidx = self.twigNext()
        self._contemporaries = []


    def eventNext(self):
        """ returns tuples (index, startHeight, endHeight, contemporaries) over all intervals 
            (startHeight,endHeight) without events, with increasing height.  The last event has 
            endHeight == max_float.  'contemporaries' is a list of (height,idx) tuples for twigs 
            that intersect the event interval (startHeight,endHeight) """
        thisidx = self._previouseventidx
        if thisidx == None: return None
        nextidx = self.twigNext()
        startheight = self.getHeight(thisidx)
        parent = self.getParent(thisidx)
        if parent != -1:
            parentheight = self.getHeight(parent)
            # add current twig to list of contemporaries
            heapq.heappush( self._contemporaries, (parentheight, thisidx) )
        # remove twigs from contemporary list when no longer intersecting
        while len(self._contemporaries)>0 and self._contemporaries[0][0] <= startheight:
            heapq.heappop( self._contemporaries )
        self._previouseventidx = nextidx
        if nextidx == None:
            return (thisidx, startheight, max_float, self._contemporaries)
        nextheight = self.getHeight(nextidx)
        # do not return 0-length event intervals, but rather call eventNext again
        if nextheight > startheight:
            return (thisidx, startheight, self.getHeight(nextidx), self._contemporaries)
        else:
            return self.eventNext()


    def treeLength(self, _useforest=False):
        """ Returns length of local tree """
        if _useforest:
            if self._forestlength > 0:
                return self._forestlength
        else:
            if self._treelength > 0:
                return self._treelength
        # update cache
        treelength = 0
        forestlength = 0
        self.twigIterator()
        while True:
            idx = self.twigNext()
            if idx == None: break
            parent = self.getParent( idx )
            if parent != -1:
                branchlength = self.getHeight( parent ) - self.getHeight( idx )
                if self.on_local_tree[ idx ] == 2: 
                    treelength += branchlength
                forestlength += branchlength
        self._treelength = treelength
        self._forestlength = forestlength
        if _useforest: 
            return self._forestlength
        else:
            return self._treelength

    ##########################################################################################
    #
    # simple samplers
    #
    ##########################################################################################

    def sampleTreePoint(self, _useforest=False):
        """ Samples point uniformly on tree.  Returns node index and height above point.
            I believe this routine takes up a fair amount of processing time.  There are
            clever algorithms to do this sampling faster, but I'm not sure they would
            really work in this context """
        length = 0
        lengths = []
        nodes = []
        self.twigIterator()
        while True:
            idx = self.twigNext()
            if idx == None: break
            parent = self.getParent( idx )
            # ignore the infinite branch above the root, and the zero branches above roots
            if parent > -1 and self.getHeight( parent ) < max_float:
                if _useforest or self.on_local_tree[ idx ] == 2:
                    length += self.getHeight( parent ) - self.getHeight( idx )
                    lengths.append( length )
                    nodes.append( idx )
        p = random.random() * length
        idx = bisect.bisect( lengths, p )
        if idx == 0:
            return nodes[idx], p
        else:
            return nodes[idx], p-lengths[idx-1]


    def sampleForestPoint(self):
        return self.sampleTreePoint(_useforest=True)

    ##########################################################################################
    #
    # modifiers
    #
    ##########################################################################################

    def getEmptyNodes(self, number):
        """ Returns a number of empty (unused) nodes.  Destroys invariant if nodes are added """
        nodes = []
        for idx, node in enumerate(self.parents):
            if node == 0:
                nodes.append( idx )
                if len(nodes) == number:
                    return nodes
        while len(nodes) < number:
            nodes.append( len(self.parents) )
            self.parents.append(0)
            self.heights.append(0.0)
            self.on_local_tree.append(0)
            self.onechild.append(-1)
            self.siblings.append(-1)
        return nodes


    def make_unifurcating(self, node, model):
        """ takes a root node, and ensures it is unifurcating.  This is used prior to
            splicing and coalescences """

        # must be a root
        assert self.getParent( node ) == -1

        # must not be a leaf
        child = self.getChild( node )
        assert child != -1

        # already unifurcating?
        if self.siblings[child] == -1: return node

        newroot = self.getEmptyNodes(1)[0]
        self.parents[ node ] = newroot
        self.parents[ newroot ] = -1
        self.heights[ newroot ] = self.getHeight( node )
        self.onechild[ newroot ] = node
        self.siblings[ newroot ] = -1
        return newroot


    def cut(self, node, height, model):
        """ Cuts node at height, and returns index of newly formed leaf and root.  Invariants destroyed """

        nodeheight = self.getHeight( node )
        nodeparent = self.getParent( node )
        leaf, root = self.getEmptyNodes( 2 )

        assert nodeparent != -1                         # can't cut branch above a root node
        assert nodeheight < height                      # must cut branch above node
        assert height < self.getHeight( nodeparent )    # must cut branch below node's parent

        # cut off lower part
        self.parents[ node ] = root

        # the node no longer has a sibling
        node_sibling = self.siblings[ node ]
        self.siblings[ node ] = -1

        # create new root node atop the lower part
        self.parents[ root ] = -1
        self.heights[ root ] = height
        self.onechild[ root ] = node
        self.siblings[ root ] = -1
        self.on_local_tree[ root ] = self.on_local_tree[ node ]

        # create new leaf node
        # (this invalidates on_local_tree up the branch; updating is postponed)
        self.parents[ leaf ] = nodeparent
        self.heights[ leaf ] = height
        self.onechild[ leaf ] = -1
        self.siblings[ leaf ] = node_sibling
        self.on_local_tree[ leaf ] = 0

        # set new sibling of parent's other child
        if node_sibling >= 0:
            self.siblings[ node_sibling ] = leaf

        # set child of node's parent
        self.onechild[ nodeparent ] = leaf
        
        return leaf, root

    def splice(self, node, branch):
        """ splices a top node into a branch at specified height.  Destroys invariant """
        
        height = self.getHeight( node )
        parent = self.getParent(branch)
        parentheight = self.getHeight(parent)

        assert self.getHeight(branch) < height
        assert height < parentheight
        
        self.parents[node] = parent
        self.parents[branch] = node

        old_sibling = self.siblings[branch]
        child = self.onechild[node]
        if child >= 0:
            # node had child; it has one more now
            assert self.siblings[child] == -1   # ensure node had just one child
            self.siblings[child] = branch
            self.siblings[branch] = child
        else:
            # node has no child; it has one now
            self.onechild[node] = branch
            self.siblings[branch] = -1

        if old_sibling >= 0:
            self.siblings[node] = old_sibling
            self.siblings[old_sibling] = node
        else:
            self.siblings[node] = -1

        # set child of parent
        self.onechild[parent] = node

    def remove(self, node):
        """ removes a root node, or (non-sample) leaf node connected to the local tree (i.e. parent is binary),
            or the leaf and root node of a floating branch """
        
        parent = self.getParent(node)
        self.parents[node] = 0

        # handle case of a root node
        if parent < 0:
            child = self.onechild[node]
            # lone node
            if child < 0: return
            sibling = self.siblings[child]
            if sibling > -1:
                # has two children -- we're splitting a tree.  probably won't happen; 
                # implement but throw error for caution
                self.parents[sibling] = -1
                self.parents[child] = -1
                self.siblings[sibling] = -1
                self.siblings[child] = -1
                assert False
            else:
                self.parents[child] = -1
                return
                 
        sibling = self.siblings[node]

        # handle case of floating branch
        if sibling == -1:
            assert self.getParent(parent) == -1
            self.parents[parent] = 0
            return
            
        assert self.getParent(sibling) == parent  # parent is binary

        # the parent now has only one child, so splice it out
        grandparent = self.getParent(parent)
        if grandparent > -1:
            self.parents[sibling] = grandparent
            self.onechild[grandparent] = sibling
            self.parents[parent] = 0
            sibling_parent = self.siblings[parent]
            if sibling_parent > -1:
                self.siblings[sibling] = sibling_parent
                self.siblings[sibling_parent] = sibling
            else:
                self.siblings[sibling] = -1
        else:
            # parent was root
            self.onechild[parent] = sibling
            self.siblings[sibling] = -1

    def pruneTree(self, min_terminal_branch_length):

        pruned = False
        for idx in range(len(self.parents)):

            # skip empty nodes
            if self.parents[idx] == 0: continue

            # skip leaves
            child1 = self.onechild[idx]
            if child1 == -1: continue

            # get other child, if any
            child2 = self.siblings[child1]

            # ensure that samples are not removed
            if child1 < self.samples:
                child1 = -1
            if child2 < self.samples:
                child2 = -1

            # get height of child if a leaf, or -1 if not
            if child1 > -1 and self.onechild[child1] == -1:
                height1 = self.getHeight( child1 )
            else:
                height1 = -1

            if child2 > -1 and self.onechild[child2] == -1:
                height2 = self.getHeight( child2 )
            else:
                height2 = -1

            # get highest child (shortest branch to leaf) in child1
            if height2 > height1:
                child1, height1 = child2, height2
            if height1 == -1: continue

            # test for minimum branch length
            branch_length = self.getHeight( idx ) - height1
            if branch_length > min_terminal_branch_length: continue

            # now remove
            self.debug("* Removing node ",idx," with short branch to child ",child1,"; branch length",branch_length)
            self.remove(child1)
            pruned = True
        
        return pruned

    ##########################################################################################
    #
    # management of invariants
    #
    ##########################################################################################

    def resetInvariant(self, index_to_translate = None):
        """ Ensures that nodes are sorted by height; that the on_local_tree array is
            populated correctly; and that emission patterns and tree length caches are invalidated """
        index_to_translate = self.sortTree( index_to_translate )
        self.resetCache()
        self.debug_list()
        self.populateLocalTreeStatus()
        self.populateSiblings()
        return index_to_translate

    def sortTree(self, index_to_translate ):
        """ ensures that nodes are numbered by increasing height.  NOTE: this means that if leaves are not (all) 
            contemporary, their index will not be invariant.  If we move to modeling time-structured data, leaf
            nodes will need an identifier.  """
        sortdata = sorted( (height, 
                            idx, 
                            self.parents[idx], 
                            self.onechild[idx], 
                            self.siblings[idx], 
                            self.on_local_tree[idx]) 
                           for idx,height in enumerate(self.heights) 
                           if self.parents[idx] != 0)
        translatedict = {-1:-1}
        # build translation dictionary
        for newidx, data in enumerate(sortdata):
            height, oldidx, parent, onechild, sibling, isancestral = data
            translatedict[oldidx] = newidx
        # build new translated arrays
        for newidx, data in enumerate(sortdata):
            height, oldidx, parent, onechild, sibling, onlocaltree = data
            self.heights[newidx] = height
            self.parents[newidx] = translatedict[parent]
            self.onechild[newidx] = translatedict[onechild]
            self.siblings[newidx] = translatedict[sibling]
            self.on_local_tree[newidx] =  onlocaltree
        # clear any trailing entries
        for idx in range(len(sortdata), len(self.parents)):
            self.parents[idx] = 0
        # return any translated index
        if index_to_translate:
            return translatedict[index_to_translate]

    def resetCache(self):
        self._emitpatterns = []
        self._treelength = -1
        self._forestlength = -1

    def populateLocalTreeStatus(self):
        """ annotate each node by whether it is on the local tree.  Worst case linear """
        for idx in range(len(self.parents)):
            self.on_local_tree[idx] = 0

        # first mark a chain to the MRCA and beyond to the root, as on local tree
        node = 0
        while node >= 0:
            self.on_local_tree[ node ] = 2
            node = self.getParent( node )

        # find MRCA
        mrca = 0
        for node in range(1,self.samples):
            while self.on_local_tree[node] == 0:
                self.on_local_tree[ node ] = 2
                node = self.getParent( node )
            mrca = max(mrca, node)

        # store local MRCA
        self.local_mrca = mrca

        # remove MRCA and its ancestors from the local tree, and mark them ancestral to samples
        self.local_tree_root = -1
        while mrca >= 0:
            self.on_local_tree[ mrca ] = 1
            self.local_tree_root = mrca
            mrca = self.getParent( mrca )


    def populateSiblings(self):
        """ populate child and sibling arrays.  Linear """

        self.onechild = [-1] * len(self.parents)
        self.siblings = [-1] * len(self.parents)

        for child in range(len(self.parents)):

            parent = self.getParent(child)
            if parent > 0:
                otherchild = self.onechild[parent]
                if otherchild == -1:
                    self.onechild[parent] = child
                else:
                    self.siblings[otherchild] = child
                    self.siblings[child] = otherchild

    ##########################################################################################
    #
    # core helpers
    #
    ##########################################################################################

    def makeEmissionPatterns(self, model):
        """ builds emission bit patterns for every node.  A 1 in position k on some node means that
            a mutation in the corresponding branch results in a derived character at tip k """

        if len(self._emitpatterns)>0: return  # reset by coalescences

        ep = [0] * len(self.parents)

        for idx in range(len(self.parents)):

            if idx < self.samples: ep[idx] = 1<<idx

            if self.getParent(idx) > 0 and self.getHeight(self.getParent(idx)) < max_float:  
                # not empty, not root, and not the infinite top branch
                ep[ self.getParent(idx) ] |= ep[idx]

        if not model.phased:
            # assume samples are ordered pairs of haploid genomes;
            # encode emissions as 0=hom ref, 1=het, 3=hom alt
            left_phase = int( (math.pow(2,model.samples)-1)/3 + 0.001 )  # bit pattern 0101010101...
            right_phase = 2*left_phase                                   # bit pattern 1010101010...
            for idx in range(len(self.parents)):
                het = (ep[idx] & left_phase) | ((ep[idx] & right_phase) >>1 )
                hom = ((ep[idx] & left_phase) <<1 ) & (ep[idx] & right_phase)
                ep[idx] = het | hom
        self._emitpatterns = ep


    def sampleCoalescences(self, model):
        """ Important subroutine of sampleNextGenealogy; samples a coalescent given a recombination point """

        self.eventIterator()

        # loop over the various coalescent events that may be necessary to join the local tree
        coalescences = []
        coalescence_root_node_1 = -1                   # placeholder for root formed after recombination
        coalescence_root_node_2 = self.local_tree_root # root of local tree
        coalescence_root_height_1 = self._recombheight 
        coalescence_root_height_2 = self.getHeight( coalescence_root_node_2 )
        coalescence_1_active = False
        coalescence_2_active = False

        # location of coalescent process
        current_time = min(coalescence_root_height_1, coalescence_root_height_2)

        # get first event
        (idx, startheight, endheight, contemporaries) = self.eventNext()

        while True:

            # use an independent draw from the standard exponential distribution for each new event
            standard_expo_variate = -math.log( random.random() )

            assert self.siblings_well_paired()
            
            # loop over the existing events to find the new coalescent point
            while True:

                # implement the SMC model: if a contemporary is the branch in which the recombination
                # occurred, then remove it, so that back coalescences will not occur
                if model.remove_recombination_branch:
                    self.debug("* smc model: looking for branch w node nr ",self._recomb_idx)
                    if self._recomb_idx in set(n for h,n in contemporaries):
                        self.debug("* smc model: found; was:",contemporaries)
                        contemporaries = [(h,n) for h,n in contemporaries if n != self._recomb_idx]
                        self.debug("* smc model: found; is:",contemporaries)

                # skip event intervals entirely below recombination
                if endheight <= current_time: 
                    (idx, startheight, endheight, contemporaries) = self.eventNext()
                    continue

                coalescence_1_active = current_time >= coalescence_root_height_1
                coalescence_2_active = current_time >= coalescence_root_height_2

                self.debug("* processing event: node=",idx," interval=",startheight,endheight,
                           "active="," 1"[coalescence_1_active]," 2"[coalescence_2_active],
                           " current time=",current_time)
            
                # calculate the rate of any coalescence occurring
                if coalescence_1_active and coalescence_2_active:
                    # lineage 1 into any contemporary (n, which does not include lineage 1),
                    # lineage 2 into any contemporary, or
                    # lineages 1 and 2 together: total 2n + 1
                    total_rate = (2*len(contemporaries) + 1) / (2.0 * model.getNe())
                else:
                    # single lineage
                    total_rate = len(contemporaries) / (2.0 * model.getNe())
                
                # calculate the interval over which coalescences can occur
                coal_interval_start = max(current_time,startheight)
                interval = endheight - coal_interval_start

                # calculate exponential variate with proper rate
                if total_rate < 1e-10:
                    expo_variate = 1e40
                else:
                    expo_variate = standard_expo_variate / total_rate
                if expo_variate > interval:
                    self.debug("* no coalescence")
                    # no coalescence event.  get a new unbiased sample of a standard exponential
                    # variate by conditioning on no event having happened in the (scaled) interval.
                    # This saves a new draw from the random number generator, makes a big difference
                    standard_expo_variate = (expo_variate - interval) * total_rate
                else:
                    # coalescence event
                    break

                # move to the next interval
                current_time = endheight
                (idx, startheight, endheight, contemporaries) = self.eventNext()

            # choose branch to coalesce into, and compute coalescent height
            if coalescence_1_active and coalescence_2_active:
                # allow either root to coalesce, or roots to coalesce with each other
                root_idx, coal_idx = \
                    random.choice( [(coalescence_root_node_1, coal_idx) for height, coal_idx in contemporaries] +
                                   [(coalescence_root_node_2, coal_idx) for height, coal_idx in contemporaries] +
                                   [(coalescence_root_node_1, -2-coalescence_root_node_2)] )  # flag for pairwise coalescence
            else:
                # choose one of the contemporary branches
                height, coal_idx = random.choice( contemporaries )
                if coalescence_1_active: root_idx = coalescence_root_node_1
                else:                    root_idx = coalescence_root_node_2

            coal_height = coal_interval_start + expo_variate

            self.debug("* coalescence of ",root_idx," into node=",coal_idx," height=",coal_height)
            
            # add coalescence from root_idx at coal_height into coal_idx branch; modifyable
            coalescences.append( [root_idx, coal_height, coal_idx] )

            # next:
            # - if root1 and root2 coalesced together, we're done
            # - if root2 coalesced, follow the branch up, and mark "ancestral",
            #      until we hit a root, and continue the process
            # - if root1 coalesced, follow the branch up until 
            #    - we hit an ancestral branch, and we're done, or
            #    - we hit a root, and continue the process

            if coal_idx <= -2:
                self.debug("*  pairwise coalescence, done")
                break

            if root_idx == coalescence_root_node_2:
                while self.getParent(coal_idx) != -1:
                    self.debug("*  (extending local tree:) node ",coal_idx," is not root; continuing")
                    coal_idx = self.getParent(coal_idx)
                self.debug("*  node ",coal_idx," is root, stopped")
                coalescence_root_node_2 = coal_idx
                coalescence_root_height_2 = self.getHeight( coal_idx )
            else:
                while self.on_local_tree[coal_idx] == 0 and self.getParent(coal_idx) != -1:
                    self.debug("*  node ",coal_idx," is not ancestral and not root; continuing")
                    self.on_local_tree[coal_idx] = 1    # mark ancestral
                    coal_idx = self.getParent(coal_idx)
                self.debug("*  node ",coal_idx," is ancestral or root, stopped")

                if self.on_local_tree[coal_idx] != 0:
                    self.debug("* node is ancestral, done")
                    break
                coalescence_root_node_1 = coal_idx
                coalescence_root_height_1 = self.getHeight( coal_idx )

            if coalescence_root_node_1 == coalescence_root_node_2:
                self.debug("** We coalesced, apparently, so STOP!")
                break

            # set time at which to continue process
            current_time = max(coal_height, min(coalescence_root_height_1, coalescence_root_height_2))
                
            self.debug("* restart coalescence at ",
                       (coalescence_root_node_1,coalescence_root_height_1,coalescence_1_active),
                       (coalescence_root_node_2,coalescence_root_height_2,coalescence_2_active))

        return coalescences


    def pruneSometimes(self, model):
        self.prune_counter += 1
        if self.prune_counter % model.prune_interval == 0:
            self.debug("* pruning")
            if self.pruneTree( model.min_terminal_branch_length ):
                self.resetInvariant()


    ##################################################################################################
    #
    # main API is below
    #
    ##################################################################################################


    def getForestLength(self):
        """ Returns total length of forest """
        return self.treeLength(_useforest=True)


    def sampleNextGenealogy(self, model):
        """ Samples a new genealogy, conditional on a recombination occurring """

        # sample recombination point
        idx, heightAbove = self.sampleForestPoint()
        self._recomb_idx = idx
        self._recombheight = self.getHeight( idx ) + heightAbove

        self.debug("* Recombination at ",self._recomb_idx,self._recombheight)
        
        # if the branch is not on the local tree, postpone the coalescence process.
        # instead, implement the recombination and return
        if (not model.wiufhein) and self.on_local_tree[self._recomb_idx] != 2:
            leafidx, rootidx = self.cut( self._recomb_idx, self._recombheight, model )
            self.resetInvariant()
            self.pruneSometimes(model)
            return

        # For efficiency, postpone the tree surgery to implement the recombination.
        # The reason is that surgery would require re-sorting the tree to allow iteration over events.
        # However, the branch above the recombination will be renamed, and coalescences into the branch
        #  (which is always above the recombination) will be renamed appropriately afterwards

        coalescences = self.sampleCoalescences(model)

        # finally, perform tree surgery
        self.debug_list()
        leafidx, rootidx = self.cut( self._recomb_idx, self._recombheight, model )
        self.debug("* implemented recombination at ",self._recomb_idx,self._recombheight,
                   ": leaf=",leafidx," root=",rootidx)
        self.debug_list()

        for coal_number in range(len(coalescences)):

            topidx, coal_height, coal_idx = coalescences[coal_number]
            
            # the index of the new free branch starting the coalescence is not known
            # before implementing the recombination, and is referred to as -1
            if topidx == -1:
                topidx = rootidx

            # rename any coalescences referring to branch _recomb_idx, which is now cut and
            # stops at the recombination, to the upper branch now referred to as leafidx
            if coal_idx == self._recomb_idx:
                coal_idx = leafidx

            # bring (an) active node up to coalescence height
            topidx = self.make_unifurcating(topidx, model)
            self.heights[topidx] = coal_height

            # deal with the possiblity that coalescences target the same branch twice.
            # this happens when two lineages are active, and both coalesce into the same third branch.
            # the third branch cannot be -1 (the branch carrying the recombination, before splitting) 
            # as the active branch above the local tree is above branch -1.
            # remedy: renumbering any reference to coal_idx to the current topidx
            for n in range(coal_number+1, len(coalescences)):
                if coalescences[n][2] == coal_idx:
                    coalescences[n][2] = topidx

            # treat the special case of two active nodes coalescing
            if coal_idx <= -2:

                coal_idx = -(coal_idx + 2)
                self.debug("* coalescing two active nodes: ",topidx," and ",coal_idx," at height",coal_height)
                self.debug_list()

                # get child of other active node
                self.debug("*  unifurcating",coal_idx)
                coal_idx = self.make_unifurcating(coal_idx, model)
                child = self.getChild( coal_idx )
                self.debug("*  result:",coal_idx,"; its child=",child)

                # remove active node
                self.debug("*  removing ",coal_idx)
                self.remove( coal_idx )

                # coalesce child into remaining active node
                self.debug("*  setting parent of ",child," to be ",topidx,"; was ",self.parents[child])
                self.parents[ child ] = topidx
                # reset siblings
                topidx_child = self.getChild( topidx )   # this is not child; 
                self.siblings[ topidx_child ] = child    # and both topidx and coal_idx are roots
                self.siblings[ child ] = topidx_child

            else:

                self.debug("* coalescing node ",topidx," into node ",coal_idx," at height",coal_height)
                # splice node into the branch we coalesced into
                # topidx is already a lone branch, so splicing will not produce trifurcations
                self.splice( topidx, coal_idx )
        
        # implement SMC and SMC' models
        if (model.remove_nonancestral_twig or model.remove_recombination_branch) and leafidx > -1:
            self.debug_list()
            self.debug("* Implementing SMC/SMC' models by removing leaf=",leafidx)
            self.remove( leafidx )
            self.debug_list()

        # reset all invariants, in particular siblings to allow deletions
        self.resetInvariant()

        # for SMC and SMC' models, remove branch above root, if any
        if (model.remove_nonancestral_twig or model.remove_recombination_branch) and \
                self.local_mrca != self.local_tree_root:
            self.debug("* removing local tree root")
            self.remove( self.local_tree_root )
            self.resetInvariant()

        assert self.siblings_well_paired()
        self.debug("* Done:")
        self.debug_list()

        self.pruneSometimes( model )


    def emissionProbability(self, mu, emission, model):
        """ Returns likelihood of [the equivalence class of] an emission.  
            For phased data, emissions are encoded as 0,1 for each sample
            For unphased data, emissions are encoded as 0, 1, 2 for hom ref, het, and hom alt. """

        self.makeEmissionPatterns( model )

        # build binary emission code; make no assumptions about root state
        anc_emission_code = 0
        der_emission_code = 0
        if model.phased:
            assert len(emission) == model.samples
            for code in reversed(emission):
                assert code in (0,1)
                anc_emission_code = anc_emission_code*2 + code
                der_emission_code = der_emission_code*2 + (1-code)
        else:
            assert model.samples % 2 == 0
            assert len(emission) == model.samples//2
            for code in reversed(emission):
                if code == 0:
                    anc_emission_code = anc_emission_code * 4
                    der_emission_code = der_emission_code * 4 + 3
                elif code == 1:
                    anc_emission_code = anc_emission_code * 4 + 1
                    der_emission_code = der_emission_code * 4 + 1
                elif code == 2:
                    anc_emission_code = anc_emission_code * 4 + 3
                    der_emission_code = der_emission_code * 4
                else:
                    assert False

        # if rooted, use that 0 = ancestral
        if model.rooted: der_emission_code = anc_emission_code

        # calculate likelihood
        likelihood = 0
        for idx in range(len(self.parents)):
            if self.getParent(idx) <= 0:  # empty or root
                continue
            if self._emitpatterns[idx] in [anc_emission_code, der_emission_code]:
                branch_length = self.getHeight( self.getParent(idx) ) - self.getHeight(idx)
                likelihood += branch_length * mu
        return likelihood


    def sampleEmission(self, model):

        """ Samples an emission from the current genealogy """
        self.makeEmissionPatterns( model )

        idx, heightAbove = self.sampleTreePoint()
        empat = self._emitpatterns[idx]

        # convert binary emission pattern into a list
        emission = []
        if model.phased:
            for i in range(model.samples):
                emission.append( empat & 1 )
                empat >>= 1
        else:
            for i in range(model.samples//2):
                code = empat & 3
                if code == 0 or code == 1:
                    emission.append(code)
                else:
                    emission.append(2)
                empat >>= 2

        return emission


##########################################################################################
#
# end: class Genealogy
#
##########################################################################################






##########################################################################################
#
# Test code
#
##########################################################################################


class Accumulator:

    def __init__(self):
        self.s0 = 0.0
        self.s1 = 0.0
        self.s2 = 0.0

    def append(self, intervals):
        for dx, y in intervals:
            self.s0 += dx
            self.s1 += dx*y
            self.s2 += dx*y*y

    def getData(self):
        return self.s0, self.s1, self.s2
        


def test_simulation( iters1, iters2, printGenealogy=False, pars={} ):
    """ many short histories """

    for i in range(iters1):
        random.seed(i)
        model = Model( **pars )
        genealogy = Genealogy( model )
        for iter in range(iters2):
            if printGenealogy:
                print "iteration:",iter
                print genealogy
                print "-"*80
            genealogy.sampleNextGenealogy( model )
        if printGenealogy:
            print "*"*100
            print "*"*100

def test_all( printGenealogy=False, pars={} ):

    for numtrees, iters1, iters2 in [#(5,10,100),
                                     (2,1000,40),
                                     (2,50,400),
                                     (2,1,4000)]:
                                     #(3,1000,40),
                                     #(3,100,400),
                                     #(3,10,4000)]:
        print "Testing %s trees, %s iterations, depth %s" % (numtrees, iters1, iters2)
        pars['samples'] = numtrees
        test_simulation( iters1, iters2, printGenealogy=printGenealogy, pars=pars )

def test_mean_tmrca(iters=100, seqlen=200000, pars={}, seed=0, 
                    mode="default", printmode=None, randomstart=False, excursion_posns=[]):

    rho = 1e-8

    samples = pars.get('samples')

    sampledata = {'parameters':["rho: %s" % rho,
                                "numtrees: %s" % samples,
                                "iters: %s" % iters,
                                "seqlen: %s" % seqlen,
                                "pars: %s" % pars]}

    excursion_num = [0 for r in excursion_posns]
    excursion_den = [0 for r in excursion_posns]
    excursion_posns.append(1e10)  # add sentinel

    for i in range(iters):

        random.seed(50*i+seed)
        model = Model( **pars )
        genealogy = Genealogy( model, randomstart=randomstart )
        excursion_i = 0

        sample   = [(0,genealogy.getTMRCA( [0,1] ))]
        if samples > 3:
            sample02 = [(0,genealogy.getTMRCA( [0,2] ))]
            sample23 = [(0,genealogy.getTMRCA( [2,3] ))]
        else:
            sample02 = sample23 = None

        lasttmrca = sample[0][1]
        iters, tmrcas = 0, 1

        pos = 0

        times = set()

        while pos < seqlen:

            f = genealogy.getForestLength()
            interval = expovariate( f*rho )
            genealogy.sampleNextGenealogy( model )
            pos += interval
            tmrca = genealogy.getTMRCA( [0,1] )
            sample.append( (pos, tmrca) )
            times.add( tmrca )

            if samples > 3:
                sample02.append( (pos, genealogy.getTMRCA( [0,2] ) ) )
                sample23.append( (pos, genealogy.getTMRCA( [2,3] ) ) )

            while excursion_posns[excursion_i] < pos:
                if lasttmrca == sample[0][1]:
                    # condition on being in the starting position
                    excursion_den[excursion_i] += 1
                    if tmrcas>1:
                        # we've made an excursion
                        excursion_num[excursion_i] += 1
                excursion_i += 1

            if tmrca != lasttmrca: 
                tmrcas += 1
                lasttmrca = tmrca

            if printmode and printmode.endswith("all"):
                print genealogy.repr_ascii( graph_output=printmode.endswith("graph") )
                print "*" * 80

            iters += 1

        if printmode and printmode.startswith("ascii"):
            print genealogy.repr_ascii( graph_output=printmode.endswith("graph") )
            print "*" * 80
            print "*" * 80

        if mode != "print":
            addsampledata( sampledata, sample, sample02, sample23, seqlen )
        else:
            printsampledata( sample, seqlen )
        # store total sequence length, and tmrca-changing recombinations
        sampledata['seqlen'] = sampledata.get('seqlen',0) + pos
        sampledata['iterations'] = sampledata.get('iterations',0) + 1
        sampledata['tmrcas'] = sampledata.get('tmrcas',0) + tmrcas
        sampledata['times'] = sampledata.get('times',0) + len(times)
        # store number of nodes at end of simulation
        numnodes = genealogy.numNodes()
        sampledata['numnodes'] = sampledata.get('numnodes',0) + numnodes
        print "Finished iteration ",i+1

    sampledata['excursion_num'] = excursion_num
    sampledata['excursion_den'] = excursion_den
    sampledata['excursion'] = excursion_posns[:-1]

    report( sampledata )


def printsampledata( sample, seqlen ):

    Ne = 10000 ## HARDCODED
    pos0 = 0
    tmrca0 = 0
    lastsegment = -1
    lasttmrca = -1
    print "//"
    for pos,tmrca in sample:
        pos = min(pos, seqlen)
        segment = pos-pos0
        pos0 = pos
        if segment>0:
            # output segment and tmrca0.
            # accumulate if tmrca0 == lasttmrca -- save disk space
            if tmrca0 == lasttmrca:
                lastsegment += segment
            else:
                if lastsegment > -1:
                    print "[%1.3f](1:%1.5f,2:%1.5f);" % (lastsegment, lasttmrca/(2*Ne), lasttmrca/(2*Ne) )
                lastsegment = segment
                lasttmrca = tmrca0
        tmrca0 = tmrca
    if lastsegment > -1:
        print "[%1.3f](1:%1.5f,2:%1.5f);" % (lastsegment, lasttmrca/(2*Ne), lasttmrca/(2*Ne) )
    

def test_mean_tmrca_ms(numtrees=2, iters=100, seqlen=200000, tolerance=1e10, mode="normal", filename=None):

    import ms_simulate
    rho = 1e-8  # also hard=coded in ms_simulate.py

    sampledata = {'parameters':["rho: %s" % rho,
                                "numtrees: %s" % numtrees,
                                "iters: %s" % iters,
                                "seqlen: %s" % seqlen,
                                "pars: %s" % {}]}

    if filename:
        correction_factor = 0.5  # self generated, units of 2Ne
    else:
        correction_factor = 1.0  # ms uses units of 4Ne

    data = ms_simulate.sim_ms(numtrees, iters, seqlen, tolerance, filename=filename, correction_factor = correction_factor)

    for i in range(len(data)):

        lasttmrca = -1
        iters, tmrcas = 0, 0
        pos = 0
        times = set()

        for pos, tmrca in data[i]:

            times.add( tmrca )

            if tmrca != lasttmrca: 
                tmrcas += 1
                lasttmrca = tmrca

            iters += 1

        if mode != "print":
            addsampledata( sampledata, data[i], None, None, seqlen )
        else:
            printsampledata( data[i], seqlen )
        # store total sequence length, and tmrca-changing recombinations
        sampledata['seqlen'] = sampledata.get('seqlen',0) + pos
        sampledata['iterations'] = sampledata.get('iterations',0) + 1
        sampledata['tmrcas'] = sampledata.get('tmrcas',0) + tmrcas
        sampledata['times'] = sampledata.get('times',0) + len(times)
        # store number of nodes at end of simulation
        numnodes = 2
        sampledata['numnodes'] = sampledata.get('numnodes',0) + numnodes
        print "Finished iteration ",i+1

    report( sampledata )


def addsampledata( sampledata, sample, sample02, sample23, seqlen ):
    """ sample is a list of (time, value) tuples.  Prepare data to report mean, variance, and autocorrelation c_ij,ij, and c_ij,ik c_ij,kl correlations """

    s = sampledata.get('intervals',[])
    s.append( makeintervals( sample, seqlen ) )
    sampledata['intervals'] = s

    autocorrelation_distances = range(0,50000,1000)
    sampledata['distances'] = autocorrelation_distances

    for s1, s2, label_offset in [ (sample,sample,0), (sample, sample02,100), (sample, sample23,200) ]:
        if s2 == None: continue
        for distance in autocorrelation_distances:
            intervals1, intervals2, intervals12 = makecorrelation( s1, s2, distance, seqlen )
            for label, data in [ (1, intervals1), (2, intervals2), (12, intervals12) ]:
                label += label_offset
                s = sampledata.get( (label,distance) , Accumulator() )
                s.append( data )
                sampledata[ (label,distance) ] = s


def makeintervals(sample, seqlen):

    intervals = []
    for idx, (x, y) in enumerate(sample[:-1]):
        dx = min(sample[idx+1][0],seqlen) - min(x,seqlen)
        appendinterval( intervals, dx, y )
    return intervals

def makecorrelation(sample1, sample2, distance, seqlen):

    intervals1 = []
    intervals2 = []
    intervals_corr = []
    points1 = [ x for (x,y) in sample1 ]
    points2 = [ x - distance for (x,y) in sample2 ]
    points = sorted( points1 + points2 )
    for idx, x in enumerate(points[:-1]):
        if x > seqlen-distance: continue
        dx = max(0,min(points[idx+1],seqlen-distance)) - max(0,min(x, seqlen-distance))
        idx1 = bisect.bisect( points1, x ) - 1
        idx2 = bisect.bisect( points2, x ) - 1
        y1 = sample1[idx1][1]
        y2 = sample2[idx2][1]
        appendinterval( intervals1, dx, y1)
        appendinterval( intervals2, dx, y2)
        appendinterval( intervals_corr, dx, y1*y2 )
    return intervals1, intervals2, intervals_corr

def appendinterval( intervals, dx, y ):

    if len(intervals)>0 and intervals[-1][1] == y:
        intervals[-1] = (dx + intervals[-1][0], y)
    else:
        intervals.append( (dx, y) )

def report( sampledata ):

    mean, var = meanvar( sampledata['intervals'] )
    print "Parameters: %s" % sampledata['parameters']
    print "Mean tmrca: %1.1f +- %1.1f" % (mean, math.sqrt(var) )
    print " Length of constant tmrca: %1.1f " % (sampledata['seqlen'] / (0.001+sampledata['tmrcas']))
    print " Independent coalescence times: %1.2f" % (sampledata['times'] / float(0.001+sampledata['iterations']))
    print " Avg End Num Nodes: %1.2f" % (sampledata['numnodes'] / float(0.001+sampledata['iterations']))

    if 'excursion' in sampledata:
        posns = sampledata['excursion']
        for idx,p in enumerate(posns):
            print "%s\t%s\t%s\t%1.3f" % (p, sampledata['excursion_num'][idx], sampledata['excursion_den'][idx],
                                            sampledata['excursion_num'][idx]/(0.001+sampledata['excursion_den'][idx]))

    distances = sampledata['distances']
    for label, indexoffset in [ ("Autocorrelation", 0), ("Correlation_ijik", 100), ("Correlation_ijkl",200) ]:
        for d in distances:
            if (indexoffset + 1, d) in sampledata:
                i1, i2, i12 = sampledata[ (indexoffset+1,d) ], sampledata[ (indexoffset+2,d) ], sampledata[ (indexoffset+12,d) ]
                m1, v1 = meanvar( i1 )
                m2, v2 = meanvar( i2 )
                m12, v12 = meanvar( i12 )
                autocorr = m12 - m1*m2
                print "%s %5s: %12.1f (fraction %2.3f) (sd %s)" % (label, d,autocorr,autocorr / (0.001+var), math.sqrt( abs(m12 - m1*m2) ))

def meanvar( intervals ):

    if type(intervals) == type([]):

        s0, s1, s2 = 0, 0, 0
        for interval in intervals:
            for dx, y in interval:
                s0 += dx
                s1 += dx*y
                s2 += dx*y*y

    else:

        s0, s1, s2 = intervals.getData()

    return s1/(s0+1e-10), s2/(s0+1e-10) - (s1/(s0+1e-10))*(s1/(s0+1e-10))


#############################################################################################
#
# Main
#
#############################################################################################

if __name__ == "__main__":

    samples = 4

    if fake_random: random = FakeRandom()

    wiufhein_pars = {'wiufhein': True, 'samples':samples}
    full_pars = {'samples':samples}
    truncated_pars_8000 = {'min_terminal_branch_length':8000, 'prune_interval':10, 'samples':samples}
    truncated_pars_4000 = {'min_terminal_branch_length':4000, 'prune_interval':10, 'samples':samples}
    truncated_pars_2000 = {'min_terminal_branch_length':2000, 'prune_interval':10, 'samples':samples}
    truncated_pars_1000 = {'min_terminal_branch_length':1000, 'prune_interval':10, 'samples':samples}
    smc_pars = {'remove_recombination_branch':True, 'samples':samples}
    smcprime_pars = {'remove_nonancestral_twig':True, 'samples':samples}

    # A min term branch length of 2000 gives results nearly identical to full model

    #test_all( pars=smcprime_pars, printGenealogy=True )

    iters = 1
    seqlen = 150000
    #seqlen = int(16/0.0004)
    #seqlen = 75000
    #seqlen = 80000
    #seqlen = 100000

    filename = None
    seed = 0

    if len(sys.argv)>2:
        try:
            seed = int(sys.argv[2])
        except:
            filename = sys.argv[2]

    mode = "print"
    #mode = "normal"

    printmode = "asciiall"     # None, "asciigraph", "ascii", "asciiall"
    randomstart = True
    excursion_posns = range(0,75001,5000)

    if sys.argv[1] == "0":
        test_mean_tmrca_ms( iters=iters, seqlen=seqlen, tolerance=0.01, mode=mode, filename=filename )
    elif sys.argv[1] == "1":
        test_mean_tmrca( iters = iters, seqlen=seqlen, pars=smc_pars, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )
    elif sys.argv[1] == "2":
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=smcprime_pars, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )
    elif sys.argv[1] == "3":
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=truncated_pars_8000, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )
    elif sys.argv[1] == "4":
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=truncated_pars_4000, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )
    elif sys.argv[1] == "5":
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=truncated_pars_2000, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )
    elif sys.argv[1] == "6":
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=truncated_pars_1000, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )
    elif sys.argv[1] == "7":
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=full_pars, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )            
    else:
        test_mean_tmrca( iters = iters, seqlen=seqlen, seed=seed,  pars=wiufhein_pars, 
                         mode=mode, printmode = printmode, randomstart = randomstart, excursion_posns=excursion_posns )            
