import copy, sys

class NewickTreeError:
    def __init__(self,err):
        self.err = err

    def __repr__(self):
        return "NewickTreeError: " + self.err


class NewickTree:
    def __init__( self, root = None ):
        """
        t = NewickTree( root )

        Creates a new instance of a newick tree with 'root object' (string).
        """
        self.root = root
        self.left = None  #  NewickBinaryTree object
        self.right = None #  NewickBinaryTree object
        self.branch = None
        self.leaves = None # only leaves : [A,B,C] for (( C, B )a, A)b;
        self.inodes = None # only internal nodes: [a,b] for (( C, B )a, A)b;
        self.other = [] #  liste of NewickTree object

    def getIntNodes( self ):
        """
        t.getIntNodes() --> objet

        Return internal nodes of tree as objet
        """
        return self.inodes

    def isSetIntNodes( self ):
        """
        t.isSetIntNodes() --> bool

        Return True if internal nodes of tree is set, False otherwise
        """
        return self.inodes is not None

    def setIntNodes( self ):
        """
        t.setIntNodes()

        Set label's internal nodes of the tree as sorted list ( root is
        considered as an internal node if tree has children ).
        """
        nodes = []
        if ( not self.isLeaf() ) and self.root is not None :
            nodes.append( self.root )
        if self.left is not None:
            if self.left.inodes is None:
                self.left.setIntNodes()
            nodes += self.left.inodes
        if self.right is not None:
            if self.right.inodes is None:
                self.right.setIntNodes()
            nodes += self.right.inodes
        if self.other != []:
            for o in self.other:
                if o.inodes is None:
                    o.setIntNodes()
                nodes += o.inodes

        nodes.sort()
        self.inodes = tuple( nodes )

    def getLeaves( self ):
        """
        t.getLeaves() --> object

        Return leaves of tree as object.
        """
        return self.leaves

    def isSetLeaves( self ):
        """
        t.isSetLeaves() --> bool

        Return True if leaves is set, False otherwise
        """
        return self.leaves is not None

    def setLeaves( self ):
        """
        t.setLeaves()

        Set label's leaves of tree.
        """
        leaves = []

        if self.isLeaf() and self.root is not None:
            leaves.append( self.root )
        if self.left is not None:
            if self.left.leaves is None:
                self.left.setLeaves()
            leaves += self.left.leaves
        if self.right is not None:
            if self.right.leaves is None:
                self.right.setLeaves()
            leaves += self.right.leaves
        if self.other != []:
            for o in self.other:
                if o.leaves is None:
                    o.setLeaves()
                leaves += o.leaves

        leaves.sort()
        self.leaves = tuple( leaves )

    def isEmpty( self ):
        """
        t.isEmpty() --> bool

        Return True if object is not empty, False otherwise
        """
        return self.root is None and self.left is None and self.right is None \
               and self.branch is None and self.other == []

    def isLeaf( self ):
        """
        t.isLeaf() --> bool

        Return True if object has no children, False otherwise
        """
        return ( self.left is None ) and ( self.right is None ) \
               and ( self.other == [] )

    def getRoot( self ):
        """
        t.getRoot() --> root object

        Return the root object (string).
        """
        return self.root

    def setRoot( self, root ):
        """
        t.setRoot( root )

        Set the root with root object (string).
        """
        self.root = root

    def hasRoot( self ):
        """
        t.hasRoot() --> bool

        Return True if the newick tree object has a root, False otherwise
        """
        return self.root is not None

    def getBranch( self ):
        """
        t.getBranch() --> digit

        Return the branch object (digit).
        """
        return self.branch

    def setBranch( self, branch ):
        """
        t.setBranch( branch )

        Set the branch with branch object (digits).
        """
        self.branch = branch 

    def hasBranch( self ):
        """
        t.hasBranch() --> bool

        Return True if the newick tree object has a branch value, False
        otherwise
        """
        return self.branch is not None

    def insertLeft( self, left ):
        """
        t.insertLeft( left )

        Creates a new NewickTree object with left as root object and put it as
        the left child of the current node.
        """
        if self.left == None:
            self.left = NewickTree( left )
        else:
            t = NewickTree( left )
            t.left = self.left
            self.left = t


    def getLeftChild( self ):
        """
        t.getLeftChild() --> object

        Return the newick tree object corresponding to the left child of the
        current object.
        """
        return self.left

    def hasLeftChild( self ):
        """
        t.hasLeftChild() --> bool

        Return True if the newick tree object has a left child, False otherwise
        """
        return self.left is not None


    def insertRight( self, right ):
        """
        t.insertRight( right )

        Creates a new NewickTree object with right as root object and put it as
        the right child of the current node.
        """
        if self.right == None:
            self.right = NewickTree( right )
        else:
            t = NewickTree( right )
            t.right = self.right
            self.right = t

    def getRightChild( self ):
        """
        t.getRightChild() --> object

        Return the newick tree object corresponding to the right child of the
        current object.
        """
        return self.right

    def hasRightChild( self ):
        """
        t.hasRightChild() --> bool

        Return True if the newick tree object has a right child, False otherwise
        """
        return self.right is not None


    def hasChildren( self ):
	"""
	t.hasChildren() --> bool

	Return True if the newick tree object has children, False otherwise
        """
	return self.left is not None and self.right is not None

    def insertOther( self, other, pos ):
        """
        t.insertOther( other, pos )

        Creates a new NewickTree object with other as root object and put it as
        the other child of the current node.
        """
        ## Si right empty, on met a right aulieu de other
        if self.right == None:
            self.right = NewickTree( other )
        elif self.other == []:
            self.other.append( NewickTree( other ) )
        else:
            t = NewickTree( other )
            t.right = ( self.other[ pos ] ) 
            self.other[ pos ] = t

    def appendOther( self, other ):
        """
        t.appendOther( other, pos ) --> int

        Creates a new NewickTree object with other as root object and append it as
        the list of other child of the current node. Return the index of the new
        NewickTree in the list.
        """
        if self.right is None:
            print >>sys.stderr,  NewickTreeError("do not insert other child before righ child")
            sys.exit()
        self.other.append( NewickTree( other ) )
        return len( self.other ) -1

    def getOtherChild( self , pos ):
        """
        t.getOtherChildren() --> list

        Return the  NewickTree object corresponding to the other child at pos of
        the current node.
        """
        return self.other[pos]


    def getOtherChildren( self ):
        """
        t.getOtherChildren() --> list

        Return list of  NewickTree object corresponding to the other children of
        the current node.
        """
        return self.other


    def hasOtherChild( self ):
        """
        t.hasOtherChild() --> bool

        Return True if the NewickBinaryTree object has a other child, False
        otherwise
        """
        return self.other != []

    def getSummary( self ):
        """
        t.getSummary() --> dict

        Return for all subtrees leaves label, branch length and root in
        dictionnary

        ex: (C,((A,B)a:1, D,(E,F)b:2)c:3)d:4;
        --> {(A,B,C,D,E,F): [d,4],(A,B,D,E,F): [c,3],(A,B): [a,1],
             (D):[None, None], (E,F),b,2)}
        """
        summary = {}

        if self.leaves is None:
            self.setLeaves()

        l = self.leaves

        if len( l ) > 1:
            summary[ l ] = [ self.root, self.branch ]

        if self.left is not None:
            tls = self.left.getSummary()
            summary.update( tls )
        if self.right is not None:
            trs = self.right.getSummary()
            summary.update( trs )
        if self.other != []:
            for o in self.other:
                tos = o.getSummary()
                summary.update( tos )
        return summary

    def Set_SelfInternalLeaves_WITH_otherBranchlabels( self, other ):
        """
        t.Set_SelfInternalLeaves_WITH_otherBranchlabels( other )

        Set all internal leave's value of T  with all corresponding branch
        labels of other ( NewickTree object ).
        """
        if self.leaves is None:
            self.setLeaves()
        if other.leaves is None:
            other.setLeaves()

        if self.leaves != other.leaves:
            print >>sys.stderr, NewickTreeError( 'Could not map informations between trees with differents leaves' )
            sys.exit()
        oSum = other.getSummary()
        self._map( oSum )

    def _map( self, oSum ):
        lvs = self.leaves
        if lvs in oSum:
            self.setRoot( oSum[ lvs ][1] )
        if self.left is not None:
            self.left._map( oSum )
        if self.right is not None:
            self.right._map( oSum )
        if self.other !=[]:
            for o in self.other:
                o._map( oSum )

    def Set_SelfInternalLeaves_WITH_bootstrap( self, bootstrap ):
        """
        t.Set_SelfInternalLeaves_WITH_bootstrap( Boot )

        Set all internal leave's value of T  with all corresponding
        bootstrap values from Bootstrap object.
        """
        boot = bootstrap.getBoot()

        s = self.getSummary().keys()
        s.sort()
        b = boot.keys()
        b.sort()
        if s != b:
            print >>sys.stderr, NewickTreeError( 'Could not map informations between trees with differents leaves' )

        self._boot( boot )


    def _boot( self, boot ):
        lvs = self.leaves
        if lvs in boot:
            self.setRoot( boot[ lvs ] )

        if self.left is not None:
            self.left._boot( boot )
        if self.right is not None:
            self.right._boot( boot )
        if self.other != []:
            for o in self.other:
                o._boot( boot )

    def _copy( self, t ):
        nT = NewickTree( t.root )
        if t.left:
            nT.left = self._copy( t.left )
        else:
            nT.left = None
        if t.right:
            nT.right = self._copy( t.right )
        else:
            nT.right = None
        nT.other = []
        if t.other:
            for o in t.other:
                nT.other.append( self._copy( o ) )

        nT.branch = t.branch
        nT.leaves = copy.deepcopy( t.leaves )
        nT.inodes = copy.deepcopy( t.inodes )
        return nT
        

    def _preorder( self ):
        """
        Visit Root node first, then recursively visit left and finally
        recursively visit right.
        """
        print self.root
        if self.left:
            self.left.preorder()
        if self.right:
            self.right.preorder()
        if self.other:
            for o in self.other:
                o.preorder()

    def _postorder( self ):
        """
        Recursively visit left and then recursively visit right and finally
        visit Root node.
        """
        if self.left:
            self.left.postorder()
        if self.right:
            self.right.postorder()
        if self.other:
            for o in self.other:
                o.preorder()
        print self.root

    def _inorder( tree ):
        if self.left:
            self.left.inorder()
        print tree.getRoot()
        if self.right:
            self.right.inorder()
        if self.other:
            for o in self.other:
                o.preorder()

    def equal( self, other, branch = False  ):
        """
        t.equal( other )

        Labels comparison. Branch comparison is optional.
        Note: 0(n!)
        """
        # label
        if self.root != other.root: return False
        # branch
        if branch:
            if self.branch  != other.branch: return False
        # nb of children
        cs, co = 0, 0
        if self.left:
            cs +=1 
        else:
            cs +=0 
        if self.right:
            cs +=1
        else:
            cs +=0
        cs += len( self.other ) 
        if other.left:
            co +=1
        else:
            co +=0
        
        if other.right:
            co +=1
        else:
            co +=0
        co += len( other.other )

        if cs != co : return False

        # test all possible permutations of subtree
        S = [ self.left, self.right ] + self.other
        O = [ other.left, other.right ] + other.other

        s = True
        for Ol in self._all_perms( O ):
            s = True
            for i in range( cs ):
                if not S[i].equal( Ol[i], branch ):
                    s = False
                    break
            if s == True: return True

        return s  == True

    def __eq__( self, other ):
        """
        t.__eq__( other ) <==> t == other

        Label and branch comparison.
        Note: 0(n!)
        """
        if isinstance(other, NewickTree):
            if self.branch != other.branch:  return False
            return self.equal( other , True )

    def _all_perms( self, str ):
        if len(str) <=1:
            yield str
        else:
            for perm in self._all_perms( str[1:] ):
                for i in range( len( perm ) + 1 ):
                    yield perm[:i] + str[0:1] + perm[i:]

    def _2str( self ):
        s = ''
        if self.isEmpty(): return s

        if self.right is not None:
            s += '('
            if self.left is not None:
                s +=  self.left._2str()
            s += ',' + self.right._2str()
            if self.other != []:
                for o in self.other:
                    s += ',' + o._2str()
            s +=  ')'

        else:
            if self.left is not None:
                s += '(' + self.left._2str() + ')'

        if self.root is not None:
            s +=  str( self.root )

        if self.branch is not None:
            s += ':' + str( self.branch )
        
        return s

    def __str__( self ):
        s = self._2str() + ';'
        return s
    
    def _print2taxoptimizer( self ):
        s = ''
        if self.isEmpty(): 
            return s
        if self.root is not None:
            s +=  str( self.root )
        if self.left is not None:
            s += '(' + self.left._print2taxoptimizer() +','
        
            if self.right is not None:
                s += self.right._print2taxoptimizer() 
                if self.other != []:
                    for o in self.other:
                        s += ',' + o._print2taxoptimizer()
            else:
                s = s[:-1]
            
            s +=  ')'
        if self.branch is not None:
            s += ':' + str( self.branch )
        
        return s

    def print2taxoptimizer( self ):
        print self._print2taxoptimizer()
    

class BootstrapError:
    def __init__( self, err ):
        self.err = err

    def __repr__( self ):
        return "BootstrapError: " + self.err


class Bootstrap:
    def __init__( self ):
        """
        new empty dictionnary 
        """
        self.boot = {}

    def setValue ( self, label ):
        if self.boot.has_key( label ):
            self.boot[ label ] += 1
        else:
            self.boot[ label ] = 1

    def setValues ( self , labels ):
        for n in labels:
            self.setValue( n )

    def getValue( self, label ):
        if self.boot.has_key( label ):
            return self.boot[ label ]
        else:
            print >>sys.stderr, BootstrapError( 'not a bootstrap value: ' + label )

    def getBoot( self ):
        return self.boot
    

    def __repr__( self ):
        st = ''
        for k,v in self.boot.iteritems():
            st += k + '\t' + str (v)  + '\n'
        return st

if __name__=='__main__':
    t = NewickTree()
    print "--> t.insertLeft(None)"
    t.insertLeft(None)
    currentTreeL = t.getLeftChild()
    print t
    print "--> t.insertRight(None)"
    t.insertRight(None)
    currentTreeR = t.getRightChild()
    print t
    print "--> currentTreeR.insertLeft('C)"
    currentTreeR.insertLeft('C')
    print "--> currentTreeR.insertRight('D')"
    currentTreeR.insertRight('D')
    print t
    print "--> currentTreeL.insertLeft('E)"
    currentTreeL.insertLeft('E')
    print "--> currentTreeL.insertRight('F')"
    currentTreeL.insertRight('F')
    print t
    
