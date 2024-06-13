# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 09:46:38 2019.

@author: Kev
"""
import numbers
import inspect
from frac import frac


class gseries():
    """Geometric series 1,k,k^2,k^3,..."""

    def __init__(self, r, c0=1):
        self.r = r
        self.c0 = c0

    def __iter__(self):
        """Return iterator."""
        self.term = self.c0  # should be new variable ???
        return self

    def __next__(self):
        """Iterate term."""
        result = self.term
        self.term *= self.r
        return result


class tseries():
    """
    Define series terms based on taylor funcions fn/n! .

    common series include
       exp  = tseries([1])
       cosh = tseries([1,0])    [f(0),f'(0),f''(0),...] where these repeat
       sinh = tseries([0,1])
       cos  = tseries([1,0,-1,0])
       sin  = tseries([0,1,0,-1])
    """

    def __init__(self, terms, cycle=None):
        self.terms = terms
        if cycle is None:
            cycle = len(terms)
        self.cycle = cycle

    def __iter__(self):
        """Iterate function."""
        idx = 0
        k = 0
        kfact = 1
        while True:
            yield frac(self.terms[idx], kfact)
            k += 1
            kfact *= k
            idx += 1
            if idx == len(self.terms):
                idx -= self.cycle


class ifun():
    """Indexed variable base on function."""

    def __init__(self, fun):
        """Allow [] notation rather than ()."""
        self.fun = fun

    def __getitem__(self, idx):
        """Allow [] notation rather than ()."""
        return self.fun(idx)

    def __call__(self, idx=None):
        """Return a iteratable generator function."""
        if idx is None:
            return self.gen()
        return self.fun(idx)

    def gen(self):
        """Generate terms from function."""
        idx = 0
        while True:
            yield self.fun(idx)
            idx += 1


class isit():
    """Istream iterator."""

    def __init__(self, stream):
        """Instance with own idx."""
        self.st = stream
        self.idx = 0

    def __iter__(self):
        """Start iterator."""
        return self

    def __next__(self):
        """Iterate Istream."""
        try:
            result = self.st[self.idx]
            self.idx += 1
        except IndexError:
            raise StopIteration
        if result is None:
            raise StopIteration
        return result


class istream():
    """Indexed stream (oxymoron) logs stream for indexing."""

    def __init__(self, p_s):
        coef = []

        p_end = None
        if isinstance(p_s, istream):
            p_stream = p_s.p_stream
            coef = p_s.p_coef
            if p_s.end is not None:
                p_end = p_s.end
        elif isinstance(p_s, numbers.Number):
            p_stream = iter([p_s])
        elif isinstance(p_s, frac):
            p_stream = iter(p_s)
        elif isinstance(p_s, list):
            p_stream = iter(p_s)
            p_end = len(p_s)
        elif inspect.isgenerator(p_s):
            p_stream = p_s
        elif inspect.isgeneratorfunction(p_s):
            p_stream = p_s()
#        elif isinstance(p_s, poly):
#            p_stream = iter(p_s)
#        elif isinstance(p_s, kpoly):
#            p_stream = p_s()
        elif isinstance(p_s, tseries):
            p_stream = iter(p_s)
        elif isinstance(p_s, gseries):
            p_stream = iter(p_s)
        elif isinstance(p_s, isit):
            p_stream = p_s
#        elif isinstance(p_s, class_or_tuple)
        else:
            print('error in istream init w/', type(p_s))
            raise (TypeError)
        self.p_stream = p_stream
        self.p_coef = coef
        self.idx = 0  # allow each istream instance to iterate
        self.end = p_end
        self.it = isit(self)

    def __call__(self):
        """Reset index for new iteration."""
        return istream(self)

    def __iter__(self):
        """Return fresh iterator."""
        return isit(self)
    #        my_idx = 0
    #        while True :
    #            try :
    #                yield self[my_idx]
    #                my_idx += 1
    #            except StopIteration :
    #                raise( StopIteration)

    def __next__(self):
        """Return Iterator method."""
        try:
            current_item = self[self.idx]
            self.idx += 1
            return current_item
        except IndexError:
            #           return
            raise (StopIteration)

    def __getitem__(self, idx=None):
        """Index into stream."""
        if idx is None:
            return isit(self)
        elif isinstance(idx, slice):
            r = []
            step = idx.step
            if step is None:
                step = 1
            for i in range(idx.start, idx.stop, step):
                r.append(self[i])
            return r
#        if idx < self.start :
#            if self.end == None :
#                raise(IndexError)
#            else :
#                return self.p_coef[self.end+1+idx]
        if idx < 0 or (self.end is not None and idx > self.end):
            return None
        while len(self.p_coef) < idx + 1:
            try:
                self.p_coef.append(next(self.p_stream))
            except StopIteration:
                self.end = len(self.p_coef) - 1
                return None
        if len(self.p_coef) < idx + 1:
            return None
        return self.p_coef[idx]

    def __len__(self):
        """Return length - return None if unknown or infinite."""
        if self.end is None:
            return 0
        else:
            return len(self.p_coef)

    def drop(self, cnt=1):
        """Remove element(s) from the beginning of istream-NULL for stream0."""
        if len(self.p_coef) < cnt or cnt < 0:
            raise (IndexError)
        self.p_coef = self.p_coef[cnt:]
        if self.end is not None:
            self.end -= cnt

    def pop(self, cnt=1):
        """Remove element(s) from the end of the istream."""
        if self.end is None or len(self.p_coef) < cnt or cnt < 0:
            raise (IndexError)
        self.p_coef = self.p_coef[:-cnt]
        self.end -= cnt
