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



from UserDict import UserDict
from aria.ariabase import *
from aria.tools import make_block

NOT_INIT = '<not initialized>'

class EntityError(Exception):
    pass

class EntityCastError(EntityError):
    pass
##     def __str__(self):
##         return 'Could not cast value. %s' % EntityError.__str__(self)

class EntityValueError(EntityError):
    pass

class Entity:

    ## for older versions, this variable did not exist. hence,
    ## during unpickling of such objects, self._callback would
    ## cause an error. handling the variable as a class variable
    ## circumvents such problems.

    _callback = None

    def __init__(self, description = None, error_message = None):

        check_type(description, NONE, STRING)
        check_type(error_message, NONE, STRING)
        
        self.__description = description
        self.__error_message = error_message
        self._setName(None)
        self.__descr_wrapped = None
        self.__err_msg_wrapped = None

    def set_default_value(self, v):

        try:
            self.validate(v)
            
        except EntityValueError:
            s = 'Entity "%s" <%s>: cannot be set to default value %s.' % \
                (self.getName(), self.__class__.__name__, str(v))
            raise EntityValueError, s
            
        self.__default = v

    def get_default_value(self):
        return self.__default

    def reset(self):
        """
        sets entity to its default value
        """
        self.__value = self.__default

    def getErrorMessage(self, value = None):

        msg = self.__error_message

        if value is not None:
            try:
                msg = self.__error_message % value
            except:
                pass

        return msg

    def setErrorMessage(self, msg):
        check_string(msg)
        self.__error_message = msg

    def getDescription(self):
        if self.__descr_wrapped is None:
            if self.__description is not None:
                
                s = make_block(self.__description,
                               AriaBaseClass.description_length)
                
                self.__descr_wrapped = '\n'.join(s)
            
        return self.__descr_wrapped

    def setDescription(self, s):
        check_string(s)
        self.__description = s
        self.__descr_wrapped = None

    def _is_valid(self, value):
        return self.is_valid(value)

    def is_valid(self, value):
        return 0

    def validate(self, value):
        if not self._is_valid(value):
            
            err_msg = self.getErrorMessage(value)
            
            if err_msg is None:
                err_msg = '%s <%s>: cannot be set to value %s.' % \
                          (self.getName(),
                           self.__class__.__name__, str(value))

            raise EntityValueError, err_msg

    def set(self, value):
        self.validate(value)
        self.__value = value

        ## try to call callback

        if self._callback is not None:
            self._callback(self)
            
    def get(self):

        try:
            return self.__value
        
        except:
            s = 'Entity %s <%s> has never been initialized.'
            raise ValueError, s % (str(self.getName()), self.__class__.__name__)

    def cast(self, v):
        """
        can be overridden to cast values to the basic data-type
        supported by the entity
        """
        
        return v

    def set_callback(self, f):
        """
        if a callback is set, it is called whenever the value
        of the entity has changed. callback's first argument
        is the entity instance.
        """

        self._callback = f

    def get_callback(self):
        return self._callback

    def is_initialized(self):
        try:
            self.get()
            return 1
        except:
            return 0

    def _setName(self, n):
        self.__name = str(n)

    def getName(self):
        if self.__name is None:
            name = '<no_name>'
        else:
            name = self.__name
            
        return name

    def __doc__(self):
        return str(self.getDescription())

    def __str__(self):
        if not self.is_initialized():
            val = NOT_INIT
        else:
            val = str(self.get())

        return val

    def __repr__(self):
        name = str(self.getName())
        clazz = self.__class__.__name__

        if not self.is_initialized():
            val = NOT_INIT
        else:
            val = str(self.get())

        descr = str(self.getDescription())

        s = '%s(name=%s, value=%s, description=%s)'
        return s % (clazz, name, val, descr)

class TypeEntity(Entity):

    def __init__(self, type, description = None, error_message = None):

        check_string(type)
        self.__type = type

        Entity.__init__(self, description, error_message)

    def is_valid(self, value):
        return is_type(value, self.__type)

class MultiTypeEntity(Entity):

    def __init__(self, types, description = None, error_message = None):

        check_tuple(type)
        self.__types = types

        Entity.__init__(self, description, error_message)

    def is_valid(self, value):
        return len([t for t in self.__types if is_type(value, t)]) > 0

class Float(TypeEntity):

    def __init__(self, description = 'Floating-point number',
                 error_message = 'Floating-point number expected.'):

        TypeEntity.__init__(self, FLOAT, description, error_message)

    def cast(self, v):
        try:
            return float(v)
        except Exception:
            err_msg = self.getErrorMessage(v)
            raise EntityCastError, err_msg

class PositiveFloat(Float):

    def __init__(self, description = 'Positive floating-point number',
                 error_message = 'Positive floating-point number expected.'):

        Float.__init__(self, description, error_message)

    def is_valid(self, value):

        return Float.is_valid(self, value) and value > 0.

class Integer(TypeEntity):

    def __init__(self, description = 'Integer',
                 error_message = 'Integer number expected.'):

        TypeEntity.__init__(self, INT, description, error_message)

    def cast(self, v):
        try:
            return int(v)
        
        except Exception, msg:
            err_msg = self.getErrorMessage(v)
            raise EntityCastError, err_msg

class PositiveInteger(Integer):

    def __init__(self, description = 'Positive Integer',
                 error_message = 'Positive integer number expected.'):

        Integer.__init__(self, description, error_message)

    def is_valid(self, value):

        return Integer.is_valid(self, value) and value > 0

class String(TypeEntity):

    def __init__(self, description = 'String',
                 error_message = 'String excepted.', optional = 0):
        
        TypeEntity.__init__(self, STRING, description, error_message)
        self.__optional = optional

    def is_valid(self, v):
        if not self.is_optional():
            return TypeEntity.is_valid(self, v)
        else:
            return v is None or TypeEntity.is_valid(self, v)

    def is_optional(self):
        try:
            o = self.__optional
        except:
            print 'DEBUG: no optional:',self.__class__.__name__,dir(self)
            
        return self.__optional

    def cast(self, v):
        try:
            return str(v)
        except Exception, msg:
            err_msg = self.getErrorMessage(v)
            raise EntityCastError, err_msg

class NonEmptyString(String):

    def __init__(self, description = 'Non-empty string',
                 error_message = 'Non-empty String excepted.',
                 mandatory = 0):
        
        String.__init__(self, description, error_message)

        self.mandatory(mandatory)

    def is_valid(self, v):
        if self.is_mandatory():
            x = v <> ''
        else:
            x = 1

        return x and String.is_valid(self, v)

    def mandatory(self, m = 1):
        check_int(m)
        self.__mandatory = m

    def is_mandatory(self):
        return self.__mandatory

class FourLetterString(String):

    def is_valid(self, v):
        valid = String.is_valid(self, v)
        
        if v is not None:
            valid = valid and len(v) <= 4
            
        return valid

class NonNegativeInt(Integer):

    def __init__(self, description = 'Non-negative integer',
                 error_message = 'Non-negative integer expected.'):

        Integer.__init__(self, description, error_message)

    def is_valid(self, value):

        return Integer.is_valid(self, value) and value >= 0

class NonNegativeFloat(Float):

    def __init__(self, description = 'Non-negative floating-point number',
                 error_message = 'Non-negative floating-points number expected.'):

        Float.__init__(self, description, error_message)

    def is_valid(self, value):

        return Float.is_valid(self, value) and value >= 0.

class ChoiceEntity(Entity):

    def __init__(self, elements, description = None, error_message = None):

        check_type(elements, LIST, TUPLE)
        self.__elements = elements

        s = str(elements)

        if error_message is None:
            error_message = 'ChoiceEntity: given choice (%s) unknown. ' + \
                            'Valid choices are: %s' % s

        Entity.__init__(self, description, error_message = error_message)

    def is_valid(self, value):

        return value in self.__elements

    def getErrorMessage(self, value=None):

        msg = Entity.getErrorMessage(self, value)

        try:
            return msg % repr(value)
        except:
            return msg

class TrueFalseChoice(ChoiceEntity):

    ## TODO: what about case-sensitivity?

    def __init__(self, description = None,
                 error_message = None):

        check_type(description, STRING, NONE)
        check_type(error_message, STRING, NONE)

        if error_message is None:
            error_message = 'Choice must be "true" or "false"'

        choices = ('true', 'false')

        ChoiceEntity.__init__(self, choices, description, error_message)

class YesNoChoice(ChoiceEntity):

    ## TODO: what about case-sensitivity?

    def __init__(self, description = None,
                 error_message = None):

        check_type(description, STRING, NONE)
        check_type(error_message, STRING, NONE)

        if error_message is None:
            error_message = 'Choice must be "%s" or "%s"' % (str(YES), str(NO))

        choices = (YES, NO)

        ChoiceEntity.__init__(self, choices, description, error_message)

class GZipChoice(ChoiceEntity):

    def __init__(self, description = None,
                 error_message = None):

        check_type(description, STRING, NONE)
        check_type(error_message, STRING, NONE)

        if error_message is None:
            error_message = 'Choice must be "%s", "%s" or "%s"' \
                            % (str(YES), str(NO), str(GZIP))

        choices = (YES, NO, GZIP)

        ChoiceEntity.__init__(self, choices, description, error_message)

class Weight(NonNegativeFloat):

    def __init__(self, description = None, error_message = None):
        if error_message is None:
            error_message = 'Weight w must be 0 < w <= 1.0; %s given.'

        NonNegativeFloat.__init__(self, description, error_message)

    def is_valid(self, value):
        return NonNegativeFloat.is_valid(self, value) and value <= 1.

class PeakType(ChoiceEntity):

    def __init__(self, description = None, error_message = None):

        if error_message is None:
            error_message = 'Value must be volume or intensity; %s given.'

        ChoiceEntity.__init__(self, ('volume', 'intensity'), description,
                              error_message)

class Path(String):

    ## If the path shall not be checked, set this to 0
    global_mandatory = 1

    DESCR = {0: 'Path [string]',
             1: 'Existing Path (enforced) [string].',
             -1: 'Existing Path (if <> '') [string].'}

    def __init__(self, description = None, error_message = None, exists = 1,
                 optional = 0):

        check_int(exists)
        check_int(optional)

        if error_message is None:
            error_message = 'Path "%s" does not exist.'

        String.__init__(self, description, error_message, optional)

        self.mandatory(exists)

    def is_valid(self, value):

        import os

        if self.__class__.global_mandatory and self.__mandatory == 1:
            exists = os.path.exists(value)

        elif self.__mandatory == -1:
            exists = os.path.exists(value) or value.strip() == ''
            
        else:
            exists = 1

        ok = String.is_valid(self, value)
        
        if not ok and self.is_optional():
            ok = value is None

        return ok and exists

    def mandatory(self, m = 1):
        
        check_type(m, INT, BOOL)

        m = {True: 1, False: 0, -1: -1}[m]

        self.__mandatory = m

        ## Depending on m, set new description, if it has
        ## not yet been given by the user.
        ## Note: if global mandatory flag is set to zero,
        ##       the value of m is irrelevant.

        descr = self.getDescription()
        
        if descr in self.DESCR.values() or descr is None:
            m = m or self.__class__.global_mandatory
            descr = self.DESCR[m]
            self.setDescription(descr)

    def is_mandatory(self):
        return self.__mandatory

class AbsolutePath(Path):

    def _abspath(self, path):
        import os
        return os.path.abspath(os.path.expanduser(path))

    def set(self, value):

        if value.strip() <> '':
            value = self._abspath(value)
        
        Path.set(self, value)

    def is_valid(self, value):
        return Path.is_valid(self, self._abspath(value))

class Settings(UserDict, AriaBaseClass):

    def __init__(self, keywords = None, default_settings = None):
        """
        if default_settings is given (must be of class Settings
        too) its settings are just copied.

        problem: if keywords have to specified in the constructor,
                 a Settings class cannot be ealisy sub-classed.
                 therefore it is better to introduce a method
                 'create' which returns the 'keywords'
                 dictionary. a sub-class could then simply override
                 that method.

                 the keywords argument is still supported, but a
                 DEPRECATED message is displayed. in future, keywords
                 should be specified by overriding the 'create'
                 method.
        """

        AriaBaseClass.__init__(self)
        UserDict.__init__(self)

        if keywords is not None:

            if not hasattr(self.__class__, 'deprecated_displayed'):
            
                ## keywords have been specified -> DEPRECATED message

                class_name = self.__class__.__name__

                s = '"keywords" argument in constructor ' + \
                    'of %s(Settings) class ' % class_name + \
                    'should not be used any more. Override "create" ' + \
                    'method instead.'

                self.deprecated(s)

                setattr(self.__class__, 'deprecated_displayed', 0)

        else:
            keywords = self.create()
            if keywords == {}:
                class_name = self.__class__.__name__
                self.warning('%s: no entities defined.' % class_name)

        check_dict(keywords)
        check_elements(keywords.values(), 'Entity')

        self.data = keywords

        for name, entity in keywords.items():
            if entity.getName() is None:
                entity._setName(name)

        if default_settings is not None:
            check_type(default_settings, 'Settings')
            self.update(default_settings)

        for key, value in self.create_default_values().items():

            entity = self.getEntity(key)
            entity.set_default_value(value)

        ## loop through settings and set entities name

        for name, entity in self.data.items():
            entity._setName(name)

    def create_default_values(self):
        """
        Override this method to implement default values.
        """
        return {}

    def create(self):
        """
        must return dict. key: settings name, value: entity
        """
        return {}

    def __getitem__(self, key):

        check_string(key)

        if not key in self.data:
            name = self.__class__.__name__
            s = '%s: entity "%s" not known.'
            raise ValueError, s % (name, key)

        entity = self.data[key]

        try:
            return entity.get()

        except ValueError:
            name = self.__class__.__name__
            s = '%s: setting "%s" exists, but has never been initialized.'
            self.error(ValueError, s % (name, key))

    def __setitem__(self, key, value):

        check_string(key)

        if not key in self.data:
            name = self.__class__.__name__
            s = '%s: entity "%s" not known.'
            raise ValueError, s % (name, key)

        entity = self.data[key]

        try:
            entity.set(value)
            
        except EntityValueError, err_msg:

            entity_type = entity.__class__.__name__
            
            s = '%s: Value "%s" for entity "%s" (%s) is invalid.\n%s'
            name = self.__class__.__name__
            self.error(EntityValueError, s % (name, str(value), key,
                                              entity_type, err_msg))

    def getEntity(self, key):

        check_string(key)

        if not self.has_key(key):
            name = self.__class__.__name__
            s = '%s: entity "%s" not known.'
            raise ValueError, s % (name, key)

        return UserDict.__getitem__(self, key)

    def as_dict(self):

        d = {}

        for key, entity in UserDict.items(self):
            try:
                d[key] = entity.get()
            except ValueError:
                d[key] = NOT_INIT

        return d

    def update(self, other):

        for k in other.keys():
            try:
                self[k] = other[k]
            except:
                pass

    def values(self):
        return map(self.__getitem__, self.keys())

    def items(self):
        return zip(self.keys(), self.values())

    def reset(self):
        """
        attempts to reset all settings with assigned
        default values.
        """

        for entity in UserDict.values(self):
            
            try:
                entity.reset()
            except:
                pass

    def str(self, name):
        """
        returns str(settings.getEntity(name))
        """

        return str(self.getEntity(name))

    def __repr__(self):

        class_name = self.__class__.__name__

        l = []

        for name in self.keys():
            entity = self.getEntity(name)
            name = "'%s' <%s>" % (name, entity.__class__.__name__)

            if not entity.is_initialized():
                val = '<not_initialized>'
            else:
                val = repr(entity.get())
            
            l.append('%s: %s' % (name, val))

        return '%s({%s})' % (class_name, ', '.join(l))
        
    __str__ = __repr__

if __name__ == '__main__':
    pass
