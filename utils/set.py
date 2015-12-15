#!/usr/bin/env python

class SetItem(object):

    """ Convenience class test the pSet code.
    """

    def __init__(self, name):
        """ Constructor.
        """
        self.name = name

    def __str__(self):
        """ String formatting.
        """
        return self.name


class Set(object):

    """ Class describing an item collection providing both sequential and
    random access. The list is the fundamental object.
    """

    def __init__(self, *items):
        """ Constructor.
        """
        self.list = []
        self.dict = {}
        self.index = 0
        #if name is None or isinstance(name,str):
        #    self.name = name
        #elif isinstance(name,SetItem):
        #    # Forgot to give name
        #    self.name = None
        #    self.add(name)
        #else: 
        #    raise Exception("Unexpected argument of type %s"%type(name))
        self.add(*items)


    def add(self, *items):
        """ Add one or more items to the collection. 
        """
        for item in items:
            self.list.append(item)
            if item.name is not None:
                self.dict[item.name] = item
        self.names = self.keys()


    def pop(self, key):
        """ Pop an item from the collection. 
        """
        item = self.dict.pop(key)
        index = self.list.index(item)
        self.names = self.keys()
        return self.list.pop(index)

    def __copy__(self):
        return self.__class__(*self.list[:])

    def __len__(self):
        """ Return the length of the collection.
        """
        return len(self.list)

    def __add__(self, other):
        """ Add (i.e concatenate) two Sets or a Set and a SetItem.
        """
        if isinstance(other,SetItem):
            items = self.list + [other]
        else:
            items = self.list + other.list
        return self.__class__(*items)


    def __getitem__(self, key):
        """ Support indexing, slicing and direct access.
        """
        if isinstance(key, str):
            return self.dict[key]
        elif isinstance(key, int):
            return self.list[key]
        elif isinstance(key, slice):
            return self.__class__(*self.list[key])

    def __call__(self, *keys):
        """ Support selection of an arbitrary subset of aliases.
        """
        items = []
        for key in keys:
            item = self[key]
            if item is not None:
                items.append(item)
        return self.__class__(*items)

    def __iter__(self):
        """ Provide iterator behavior.
        """
        return self

    def next(self):
        """ Provide iterator behavior.
        """
        if self.index >= len(self):
            self.index = 0
            raise StopIteration
        item = self[self.index]
        self.index += 1
        return item

    def keys(self):
        """ Ordered list of keys """
        return [ val.name for val in self.list ]
        
    def vals(self):
        """ Ordered list of values """
        return self.list

    def items(self):
        """ Ordered list of items """
        return [ (val.name, val) for val in self.list ]



    
if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()

    s = Set(SetItem('hello'), SetItem('world'))
    for item in s:
        print item
    print s[0]
    print s['world']
