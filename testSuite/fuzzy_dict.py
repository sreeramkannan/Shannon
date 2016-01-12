import doctest

class FuzzyDict(object):
    """A dictionary that matches strings with a certain
    tolerance. That is, string keys map to objects, but
    when a certain key is queried, a list of all objects
    with "similar" keys is returned.

    A key is similar if it is identical under at most
    TOLERANCE substitutions.

    This is done by using splitting the key into multiple
    portions. If any of the portions match a previously
    inserted portion, the full key is looked up and
    compared for substitution distance.

    Test normal dictionary functionality.

    >>> f = FuzzyDict()
    >>> f['ABCD'] = 1
    >>> f['ABCD']
    [1]
    >>> f['LMNO'] = 10
    >>> f['LMNO']
    [10]
    >>> f['ABCD'] = 2
    >>> f['ABCD']
    [2]

    Test fuzziness.
    >>> f['ABCX'] = 3
    >>> f['AXCX'] = 4
    >>> f['XXCX'] = 5
    >>> sorted(f['ABCD'])
    [2, 3, 4]
    """
    TOLERANCE = 2 #Tolerance in key searching
    PARTS = 3 #Number of parts to split each key into

    def __init__(self):
        self.part_dicts = []
        for _ in range(FuzzyDict.PARTS):
            self.part_dicts.append({})
        self.dict = {}

    def __setitem__(self, key, value):
        for p in range(FuzzyDict.PARTS):
            part = self.part(p, key)
            if part not in self.part_dicts[p]:
                self.part_dicts[p][part] = set()
            self.part_dicts[p][part].add(key)
        self.dict[key] = value

    def __getitem__(self, key):
        partial_keys = set()
        for p in range(FuzzyDict.PARTS):
            part = self.part(p, key)
            if part in self.part_dicts[p]:
                partial_keys.update(self.part_dicts[p][part])
        matched_keys = (k for k in partial_keys if self.matches(k, key))
        return [self.dict[k] for k in matched_keys if k in self.dict]

    def __contains__(self, key):
        return len(self.__getitem__(key)) > 0
        
    def matches(self, key1, key2):
        if len(key1) != len(key2):
            return False
        differences = 0
        for char in range(len(key1)):
            if key1[char] != key2[char]:
                differences += 1
        if differences > FuzzyDict.TOLERANCE:
            return False
        return True

    def part(self, p, key):
        return key[(p*len(key)) / FuzzyDict.PARTS:
            ((p+1)*len(key)) / FuzzyDict.PARTS]

#Note: This dictionary takes space proportional to the number of keys
#ever inserted, not the number of keys actually existing in the dictionary.

doctest.testmod()