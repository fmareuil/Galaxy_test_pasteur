"""

          ARIA -- Ambiguous Restraints for Iterative Assignment

                 A software for automated NOE assignment

                               Version 2.3


Copyright (C) Benjamin Bardiaux, Michael Habeck, Therese Malliavin,
              Wolfgang Rieping, and Michael Nilges

All rights reserved.


NO WARRANTY. This software package is provided 'as is' without warranty of
any kind, expressed or implied, including, but not limited to the implied
warranties of merchantability and fitness for a particular purpose or
a warranty of non-infringement.

Distribution of substantively modified versions of this module is
prohibited without the explicit permission of the copyright holders.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""



from aria.TypeChecking import *
from aria.ariabase import AriaBaseClass

class XMLReaderError(Exception):
    pass

class XMLTagError(XMLReaderError):

    def __str__(self):
        return 'XML tag "%s" not known or missing.' % \
               Exception.__str__(self)

class XMLDocument:

    def __init__(self, doc_type, dtd):
        self.__doc_type = doc_type
        self.__dtd = dtd

    def get_doc_type(self):
        return self.__doc_type

    def get_dtd(self):
        return self.__dtd

class XMLElement:

    def __init__(self, name = None, attr = None, tag_order = None):

        self.__cdata = None

        if name is not None:
            self.__name = name

        if attr is not None:
            for name, value in attr.items():
                setattr(self, name, value)

        self.set_tag_order(tag_order)

    def set_cdata(self, data):
        if self.__cdata is None:
            self.__cdata = data
        else:
            self.__cdata += data

    def get_cdata(self):
        return self.__cdata

##     def get_name(self):
##         raise
##         attr = '_%s__name' % self.__class__.__name__

##         if attr in self.__dict__:
##             return self.__dict__[attr]
##         else:
##             return None

    def set_tag_order(self, o):
        self.__tag_order = o

    def get_tag_order(self):
        return self.__tag_order

    def __getattr__(self, name):

        attr = '_%s__name' % self.__class__.__name__

        if attr in self.__dict__:
            my_name = self.__dict__[attr]
        else:
            my_name = None
        
        name = str(my_name) + '.' + name
        
        raise XMLTagError(name)

    def __str__(self):

        return 'XMLElement(name=%s, tag_order=%s)' % \
               (self.get_name(), str(self.get_tag_order()))

class ContentConverter:

    """
    class which interprets and converts XMLElements instances
    into Python data stuctures.
    """

    dtd_name = []

    def set_dtd_name(self, name):

        if name is None:
            raise

        ContentConverter.dtd_name.append(name)

    def release(self):

        if ContentConverter.dtd_name:
            ContentConverter.dtd_name = ContentConverter.dtd_name[:-1]
        else:
            raise 'Inconsistency: ContentConverter cannot be released.'

    def get_dtd_name(self):
        """
        returns the value of a class variable 'dtd_name'.
        since definition of element names may differ
        in different DTDs, the main purpose of the method
        is to determine the DTD within the '_xml_state'
        method in order to know, how to map an object
        to an XML element structure.
        """

        return ContentConverter.dtd_name[-1]

        return ContentConverter.dtd_name#self.__class__.dtd_name

    def _xml_state(self, x):
        """
        converts the python data structure, x, into
        a instance of XMLElement.
        """
        raise NotImplementedError

    def load_from_element(self, e):
        """
        converts an element, 'e', into a Python data structure.
        """
        return e

class ContentHandler(ContentConverter):

    """
    Abstract base class for content handler. Event-driven parsers
    must be subclassed from ContentHandler.
    """
    
    def startElementHandler(self, name, attr):
        raise NotImplementedError

    def endElementHandler(self, name):
        raise NotImplementedError

    def charDataHandler(self, data):
        raise NotImplementedError

class BaseReader:

    """
    Abstract reader class
    """

    def __init__(self, handler):
        self.__handler = handler

    def createParser(self):
        raise NotImplementedError

    def getContentHandler(self):
        return self.__handler

class XMLContentHandler(ContentHandler):

    def __init__(self):
        self.elements = []
        self.current = None

    def setDocument(self, d):
        """
        sets the main document of the xml file, that is the
        root of the hierarchy of XMLElement instances.
        """
        self.elements = [d]
        self.current = d

    def startElementHandler(self, name, attr):

        try:
            name = str(name).split(':')[1]
        except:
            pass

        new_element = XMLElement(name, attr)
        self.elements.append(self.current)
        self.current = new_element

    def endElementHandler(self, name):

        try:
            name = str(name).split(':')[1]
        except:
            pass

        if not len(self.elements):
            return

        e = self.load_from_element(name, self.current)
        self.current = self.elements.pop()

        try:
            value = getattr(self.current, name)
            try:
                value.append(e)
            except:
                setattr(self.current, name, [value, e])

        except:
            setattr(self.current, name, e)

    def charDataHandler(self, data):

        data = data.strip()
        
        if not len(data):
            return

        self.current.set_cdata(str(data))

    def load_from_element(self, name, e):
        return e

class XMLReader(BaseReader):
    
    """
    is able to read xml files and converts it into
    a hierarchy of XMLElements.
    """

    def __init__(self, content_handler):
        BaseReader.__init__(self, content_handler)
        
    def createParser(self):
        """
        creates the parser
        """

        content_handler = self.getContentHandler()
        
        if content_handler is None:
            raise 'No content handler set.'

        try:
            import xml.parsers.expat as expat
            parser = expat.ParserCreate()

        except:
            from aria.xmlparser import SelfmadeXMLParser as Parser

            print 'Could not import the Python module xml.parsers.expat. Please check your Python setup. Using Python implementation of XML parser. This could lead to a substancial performance loss.'
            
            parser = Parser()

        parser.StartElementHandler = content_handler.startElementHandler
        parser.EndElementHandler = content_handler.endElementHandler
        parser.CharacterDataHandler = content_handler.charDataHandler

        return parser

    def parse_doc_type(self, file):

        seek = file.tell()

        doctype_found = 0

        seek_previous = seek
        
        while not doctype_found:

            line = file.readline()
            doctype_found = line.find('<!DOCTYPE') > -1

            if seek_previous == file.tell(): break

            seek_previous = file.tell()
            
        file.seek(seek)

        if not doctype_found:
            tag, dtd = None, None

        else:
            tag = line.split()[1]
            dtd = line.split()[3][1:-2]

        try:
            tag = str(tag).split(':')[1]
        except:
            pass

        return tag, dtd

    def load(self, filename, gzip = 0):
        
        import os
        
        filename = os.path.expanduser(filename)

        if gzip:
            from aria.tools import gzip_open as open_func
        else:
            open_func = open

        f = open_func(filename)

        doc = XMLDocument(*self.parse_doc_type(f))

        handler = self.getContentHandler()

        if doc.get_dtd() is not None:
            handler.set_dtd_name(doc.get_dtd())

        handler.setDocument(doc)

        parser = self.createParser()
        parser.ParseFile(f)
        f.close()

        if doc.get_dtd() is not None:
            handler.release()

        failed = 0

        if len(handler.elements) <> 1:
            failed = 1
        elif handler.elements[0] <> doc:
            failed = 1

        if failed:
            raise XMLReaderError, 'XML document misformatted.'

        return doc

class XMLPickler(XMLReader):

    """
    full pickler. reads and writes xml files.
    """

    debug = 0
   
    def _dumps(self, x, name = None, tags = None):
            
        if isinstance(x, XMLElement):
            
##             try:
##                 name = x.get_name()

##             except:
##                 pass

            result = self._dumps(x.__dict__, name, x.get_tag_order())

        elif type(x) == type(self):
            if hasattr(x, '_xml_state'):
                if self.__class__.debug:
                    print x._xml_state.im_class.__name__
                    
                x = x._xml_state(x)
                
            else:
                x = x.__dict__
                
            result = self._dumps(x, name)

        elif type(x) == type({}):

            attribs = ''
            lines = []

            if tags is None:
                tags = x.keys()

            for tag in tags:
                
                try:
                    if tag[:11] == '_XMLElement':
                        continue
                    elif tag[:12] == '_XMLDocument':
                        continue
                except:
                    pass

                block = self._dumps(x[tag], tag)

                if type(block) == type(''):
                    attribs += str(block)
                    
                else:
                    lines += block

            if name is not None:

                if attribs:
                    attribs = ' ' + attribs

                if lines:
                    open_tag = '<' + name + attribs[:-1] + '>'
                    close_tag = '</' + name + '>'
                    lines = ['  ' + s for s in lines]
                    lines.insert(0, open_tag)
                    lines.append(close_tag)
                
                else:
                    open_tag = '<' + name + attribs[:-1] + '/>'
                    lines = [open_tag]
                    
            result = lines

        elif type(x) in (type([]), type(())):

            lines = []

            for i in range(len(x)):

                value = x[i]

                if type(value) in (type(self), type({})):
                    block = self._dumps(value, name)
                    if type(block) == type(''):
                        block = ['<' + name + ' ' + block[:-1] + '/>']
                        
                else:
                    block = ['<' + name + ' value="' + str(value) + '"/>']
                
                lines += block

            result = lines

        else:

            ## represent None as empty string
            
            if x is None:
                x = ''
            result = name + '="' + str(x) + '" '

        return result
                
    def dumps(self, x):
        if type(x) not in (type({}), type(self)):
            s = 'Only instances and dicts can be pickled as XML.'
            raise Exception, s

        if isinstance(x, XMLDocument):
            lines = self._dumps(x)
            doc_type = x.get_doc_type()
            dtd = x.get_dtd()

            header = '<!DOCTYPE %s SYSTEM "%s">' % (doc_type, dtd)
            lines.insert(0, header)
            
        else:
            lines = self._dumps(x)

        return '\n'.join(lines)

    def dump(self, x, filename, gzip = 0):
        import os

        filename = os.path.expanduser(filename)
        
        if gzip:
            from aria.tools import gzip_open as open_func
        else:
            open_func = open

        xml = self.dumps(x)

        f = open_func(filename, 'w')
        f.write(xml)
        f.close()

## TODO: move that one to AriaXML.py, change name to AriaBaseXMLContentHandler

class XMLBasePickler(AriaBaseClass, ContentConverter):
    pass

