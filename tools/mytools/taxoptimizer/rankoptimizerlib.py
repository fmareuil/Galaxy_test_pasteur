#! /usr/local/bin/python

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# maufrais@pasteur.fr
#

# version 1.0

import os, sys, getopt
try:
    LIB = os.environ['RANKOPTIMIZERLIB']
except:
    LIB = '/usr/local/bin'
    

if LIB not in sys.path:
    sys.path.append( LIB )

class RankOptimizerError:
    def __init__( self, err ):
        self.err = err

    def __repr__( self ):
        return "[rankoptimizer] " + self.err

class Taxon:
    def __init__( self, name=None, rank='' ):
        """
        t = Taxon( name )

        Creates a new instance of a binary tree with 'name object' (string).
        """
        self.name = name
        self.parent = None  # Taxon object
        self.childs = []    # list of Taxon object
        self.queriesS = []   # list of tuple [ (query, posLine) ]
        self.repr = ''
        self.rank = rank
        self.nbQuerys = 0
        

    def hasRank( self ):
        return self.rank

    def hasChild( self, child ):
        """
        t.hasChild( val ) --> bool
        """
        for chld in self.childs:
            if chld.name == child:
                return True
        return False

    def hasChilds( self ):
        """
        t.hasChilds() --> bool

        Return true if t has a list of childs, false otherwise. A child, it's a
        Taxon object.
        """
        return not self.childs == []

    def hasQueriesS( self ):
        """
        t.hasSpecies() --> bool
        """
        return self.queriesS != []


    def hasQuery( self, name ):
        for qry, pos in self.queriesS:
            if qry == name:
                return True
        return False

    def addQueries( self, queryAll ):
        self.queriesS.append( queryAll ) # queryAll = (query,  posLine )


    def addChild( self, child, rank='' ):
        """
        t.addChild( child )

        Creates a new Taxon object with child as name and append it as
        the childs list of t.
        """
        t = Taxon( child, rank )
        t.parent = self
        if self.childs == []:
            #t.repr = '|\t'
            t.repr = '|'
            self.childs = [ t ]
        else:
            #self.childs[-1].repr = '|\t'
            self.childs[-1].repr = '|'
            #t.repr = '.\t'
            t.repr = '.'
            self.childs.append( t )


    def getChild( self, child ):
        for chld in self.childs:
            if chld.name == child:
                return chld

            
  
##################### Krona
class ElementXML:


    def __init__( self, outfh = None, indent = None ):
        self.outfh = outfh
        self.indent = indent  # xml file indentation

    def startTag( self, element, attr = ''):
        # <element attributes>
        print >>self.outfh, '%s<%s%s>' % ( self.space(), element, attr )

    def endTag( self, element ):
        #</element>
        print >>self.outfh, '%s</%s>' % ( self.space(), element )

    def completeTag( self, element, value, attr=''):
        # <element attributes>value</element>
        print >>self.outfh, '%s<%s%s>%s</%s>' % ( self.space(), element, attr, self.html_str( value ), element )

    def newLine( self ):
        print >>self.outfh, ''

    def increaseIndent( self ):
        self.indent +=1

    def decreaseIndent( self ):
        self.indent -=1

    def getIndent( self ):
        return self.indent

    def space( self ):
        return '  ' * self.indent

    def html_str( self, value ):
        if isinstance( value, str ):
            #value = value.replace( '&', '&amp;' )
            value = value.replace( '>', '&gt;')
            value = value.replace( '<', '&lt;')
            return value
        else:
            return str( value )


class KronaDTD ( ElementXML ):

    def __init__( self, outfh = None, indent = None, collapse = 'true', key='true' ):
        ElementXML.__init__( self, outfh = outfh, indent = indent)
        self.global_attributes = { 'krona_attribute': {'collapse': collapse, 'key': key},
                                'attributes_attribute':{'magnitude': 'reads' },
                                 }

    def headerHtml( self,  kronaURL):
        print >>self.outfh, """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <title>krona</title>
  <meta charset="utf-8"/>
  <link rel="shortcut icon" href="%s/img/favicon.ico"/>
    <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"%s\"."}</script>
    <script src="%s/src/krona-2.0.js"></script>
</head>
<body>
  <img id="hiddenImage" src="%s/img/hidden.png" style="display:none"/>
  <img id="loadingImage" src="%s/img/loading.gif" style="display:none"/>
  <noscript>Javascript must be enabled to view this page.</noscript>
  <div style="display:none">
        """ % ( kronaURL, kronaURL, kronaURL, kronaURL, kronaURL)
        
    def footerHtml(self):
        print >>self.outfh, """
        </div></body></html>
        """
    
    def startKrona( self ):
        # <krona ....>
        #krona_attribute={'collapse':'true', 'key':'true' }
        self.startTag( 'krona', self._attr2str( self.global_attributes['krona_attribute'] ) )
        self.increaseIndent()

    def endKrona( self ):
        # </krona>
        self.decreaseIndent()
        self.endTag( 'krona' )
        
    def color( self, color_values ):
        # <color attribute="..." valueStart="..." valueEnd="..." hueStart="..." hueEnd="..." ...></color>
        # color_values = {attribute:'',  valueStart:'', valueEnd:'', hueStart:'',  hueEnd:''} 

        self.completeTag( 'color', '', self._attr2str( color_values ) )

    def datasets( self, datasets_values ):
        #datasets_values:[set1, set2, ...]
        self.startTag( 'datasets' )
        self.increaseIndent()
        for lds in datasets_values:
            self.dataset( lds )
        self.decreaseIndent()
        self.endTag( 'datasets' )
    
    def dataset( self, value ):
        #<dataset>set1</dataset>
        self.completeTag( 'dataset', value )

    def attributes(self, attributes_values):
        #<attributes magnitude="sample_attr_1">...</attributes>
        #attributes_values={'sample_attr': [(sample_attr_name, attribute_values)], }
        #                   'sample_list': [(sample_list_name)],
        self.startTag( 'attributes', self._attr2str( self.global_attributes['attributes_attribute'] ) )
        self.increaseIndent()
        for san, av in attributes_values['sample_attr']:
            #<attribute ...>sample_attr_1</attribute>
            self.attribute( san,av )
        for sln in attributes_values['sample_list']:
            #<list>sample_list_1</list>
            self.listf( sln )
        self.decreaseIndent()
        self.endTag( 'attributes' )
    
        
    def attribute( self, value, attribute_values ):
        # <attribute display='' mono='' listNode='' listAll='' hrefBase='' target=''>...<attribute>
        # attribute_values ={display:'' mono:'' listNode:'' listAll:'' hrefBase:'' target:''}
        self.completeTag( 'attribute', value, self._attr2str( attribute_values ) )

    
    def listf(self, value ):
        #<list>sample_list_1</list>
        self.completeTag( 'list', value )
    

    def node(self, child_nodes): #taxon
        #<node name='' href=''></node>
        #node_attributes={'name:'','href':''}
        node_attributes, samples_values, list_values, child_nodes = self._toKronaNode(child_nodes)
        self.startTag( 'node', self._attr2str( node_attributes ) )
        self.increaseIndent()
        #<sample_attr_1>
        #samples_values = [(name, [(v1,href),(v,href)...])] plusieurs val, si plusieurs dataset
        for san, sav in samples_values:
            self.sample_attr( san, sav )
        #<sample_list_1>
        
        for c in child_nodes:
            self.node( c )
        for sln, slv in list_values:
            self.sample_list( sln, slv )
        self.decreaseIndent()
        self.endTag( 'node' )
    
    def _toKronaNode( self, taxon ):
        node_attributes = {'name':taxon.name }
        if taxon.rank:
            samples_values_rank = ('rank', [ (taxon.rank,{}) ] ) 
        else:
            samples_values_rank = ('rank', [] ) 
      
        child_nodes = taxon.childs
        #for c in taxon.childs:
        #    self.toKronaNode( c  )
        
        #samples_values_read = ('reads',[ (len(taxon.queriesS), {}) ])
        samples_values_read = ('reads',[ (taxon.nbQuerys, {}) ])
        #samples_values_blast = ('blast',[ (taxon.nbQuerys, {}) ])
        list_values_read = []
        list_values_blast = []
        if len( taxon.childs ) == 0:
            if taxon.hasQueriesS():
                for query in taxon.queriesS:
                    list_values_read.append(("%-40s\t%s"% (query[0],query[1]) ,{}))
                    #list_values_blast.append((query[1],{}))
        #list_values = [('read_members',list_values_read), ('blast_members',list_values_blast)]
        list_values = [('read_members',list_values_read)]
        
        #return node_attributes, [samples_values_read, samples_values_blast, samples_values_rank], list_values, child_nodes
        return node_attributes, [samples_values_read, samples_values_rank], list_values, child_nodes

    
    def sample_attr( self, sample_attr_name, sample_attr_values ):
        
        self.startTag( sample_attr_name )
        self.increaseIndent()
        for vl, av in sample_attr_values: # plusieurs val, si plusieurs dataset
            #<val href=''>7</val>
            self.val( vl, av )
        self.decreaseIndent()
        self.endTag( sample_attr_name )
        
    def val( self, value, val_attribute ):
        #<val href=''>12</val>
        self.completeTag( 'val', value, self._attr2str( val_attribute ) )

    
    def sample_list( self, sple_list_name, sple_list_values ):
        self.startTag( sple_list_name )
        self.increaseIndent()
        self.startTag( 'vals' )
        self.increaseIndent()
        for vl,av in sple_list_values:
            #<val herf=''>7</val>
            self.val( vl,av )
        self.decreaseIndent()
        self.endTag( 'vals' )
        self.decreaseIndent()
        self.endTag( sple_list_name )

    def _attr2str( self, attrs ):
        s = ''
        for k,v in attrs.items():
            if v:
                s += ' %s="%s"' % ( k, str(v) )
        return s

    def _setAttr( self, elem, attr, value ):
        self.attributes[elem][attr] = value

class Krona:
    def __init__( self,  outfh = None, fileName=None, taxoTree = None, kronaURL= 'http://krona.sourceforge.net' ):
        """
        * Object for translating one treeobject into one Krona xml file.
        """
        
        self.outfh = outfh
        self.indent = 0
        self.taxoTree = taxoTree
        self.fileName = fileName
        self.kronaURL = kronaURL
        
    def krona(self):
        self.indent = 0
        xmlKrona = KronaDTD( self.outfh, self.indent)
        xmlKrona.startKrona()
        xmlKrona.datasets( [self.fileName] )
        #color_values = {'attribute':'', 'default':'true'}
        #xmlKrona.color( color_values)
        attributes_values={'sample_attr': [('reads', {'display':'Nb of reads', 'listAll':'read_members'}),
                                           #('blast', {'display':'Blast offset', 'listAll':'blast_members'}),
                                           ('rank', {'display':'Rank', 'mono':'true'})],
                           #'sample_list':[]
                           'sample_list': ['read_members',
                                           'blast_members'] 
                           }
        
        xmlKrona.attributes( attributes_values )

        xmlKrona.node( self.taxoTree )
        
        xmlKrona.endKrona()
        
    def kronaHtml(self):
        self.indent =0
        xmlKrona = KronaDTD( self.outfh, self.indent)
        xmlKrona.headerHtml( self.kronaURL)
        xmlKrona.startKrona()
        xmlKrona.datasets( [self.fileName] )
        #color_values = {'attribute':'', 'default':'true'}
        #xmlKrona.color( color_values)
        attributes_values={'sample_attr': [('reads', {'display':'Nb of reads', 'listAll':'read_members'}),
                                           #('blast', {'display':'Blast offset', 'listAll':'blast_members'}),
                                           ('rank', {'display':'Rank', 'mono':'true'})],
                           #'sample_list':[]
                           'sample_list': ['read_members',
                                           'blast_members'] 
                           }
        
        xmlKrona.attributes( attributes_values )

        xmlKrona.node( self.taxoTree )
        
        xmlKrona.endKrona()
        xmlKrona.footerHtml()
    #####################
  
    
import newickTree_rankopti 
def toNewickTree( taxon, tree, Acc ):
    p = 0
    #print len(taxon.childs)
    for taxa in taxon.childs:
        if not tree.left:
            tree.insertLeft( None )
            #tree.left.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '')+'#' + taxa.rank +':' + str( len(taxa.queriesS)) 
            tree.left.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '')+'#' + taxa.rank +':' + str( taxa.nbQuerys)
            tree.left = toNewickTree( taxa, tree.left, Acc)
        elif not tree.right:
            tree.insertRight( None )
            #tree.right.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '')+'#' + taxa.rank +':' + str( len(taxa.queriesS))
            tree.right.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '')+'#' + taxa.rank +':' + str( taxa.nbQuerys)
            tree.right = toNewickTree( taxa, tree.right, Acc)
        else:
            p = tree.appendOther( None )
            #tree.other[p].root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '')+'#' + taxa.rank +':' + str( len(taxa.queriesS))
            tree.other[p].root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '')+'#' + taxa.rank +':' + str( taxa.nbQuerys)
            tree.other[p] = toNewickTree( taxa, tree.other[p], Acc)
            p+=1
    return tree
    

def toDnd( taxon, Acc=False ):
    tree = newickTree_rankopti.NewickTree( 'root' )
    return  toNewickTree( taxon, tree, Acc )


def _toTree( taxon, name=True, s='', sv='', sd='.', qrName=False ):
    if name:
        if taxon.rank:
            s +=  '+'  + taxon.name + ' (' + taxon.rank + ')'  + ';'
        else:
            s +=  '+'  + taxon.name + ';'
  
    if len( taxon.childs ) > 1 :
        s += '#' + str( taxon.nbQuerys ) +'\n'
        for c in taxon.childs:
            s, sv =  _toTree( c, True, s + sd, sv, sd + c.repr, qrName )
        
    elif len( taxon.childs ) == 1:
        c = taxon.childs[0]
        if c.rank:
            s += c.name + ' (' + c.rank + ')'  +';'
        else:
            s += c.name +   ';'
        s,sv = _toTree(c, False, s, sv , sd, qrName )
    else:
        #s += '#' + str( len(taxon.queriesS)) + '\n'
        s += '#' + str( taxon.nbQuerys ) + '\n'
        if taxon.hasQueriesS():
            if qrName:
                for query,posLine  in taxon.queriesS:
                    s += "%s\n" % ( (sd + ' - ' + query + ' ( ' + str(posLine) + ' ) ') )
            #else:
            #    s += "%s\n" % ( (sd + ' - ' +  str( len(taxon.queriesS)) + ' queries') )
    
    return s,sv

def toTree( taxon, name=True, s='', sv='', sd='.', qrName=False ):
    s, sv = _toTree( taxon, name, s, sv, sd, qrName )
    s += sv + '\n'
    return s





    
