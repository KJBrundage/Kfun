# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:12:33 2020

@author: Kev
"""

"""
kite's define streams of K terms as constants + recursive offset values

pi = kite([(1,0),(4,1),(1,3),(ite(0,[1],[1]),ite(2,[],[1]))],1)
   = 0 + 4/(1+1/(3+i^2/(2i+1+...

e = kite([1,0,1,1,(1,ite(2,[],[1,0,0]))],3)
  = 1+1/(1/(1+1/(1+1/(2+1/(1+1/(1+1/(4+1/(1+1/(1+1/(6+1/(1+1/(1+1/(8+ ...
                                                                   
phi = kite([1, 1], 1)

sqrt(2) = kite([1,2], 1)

sqrt(n) = kite([m, (n-m**2,2*m)],1)

any square root will be a simple repeating palindrome
solution to quadratic surd <=> simple repeating K

tan = kite([0, (x,1), (-x**2,ite(2,[],[1]))],1) {where x = poly([0,1])}
"""

class ite():
    "iterator element"
    def __init__(self, const, r=[], s=[] ):
        "create iterator element"
        self.const = const
        if len(r) < len(s):
            r += [0]*(len(s)-len(r))
        elif len(r) > len(s):
            s += [0]*(len(r)-len(s))
        self.r     = r
        self.s     = s
        self.len   = len(r)
    def __repr__(self):
        "reference to ite"
        return "ite( "+str(self.const)+","+str(self.r)+","+str(self.s)+")"
class kite():
    def __init__(self, items=[], cycle=1):
        "2-d iterator class. a's and b's are ite's"
        self.items = items
        self.cycle = cycle
        self.longx = 0
        for ix in items:
            if isinstance(ix,tuple):
                h,k = ix
            else:
                h = 1
                k = ix
            if isinstance(h,ite) and h.len > self.longx :
                self.longx = h.len
            if isinstance(k,ite) and k.len > self.longx :
                self.longx = k.len
    def __iter__(self):
        "iterator"
        ix  = 0
        h,k = 0,0
        if isinstance(self.items[ix],tuple):
            h,k = self.items[ix]
        else:
            h,k = 1,self.items[ix]
        out = [(h,k)]
        yield (h,k)
        while True:
            ix += 1
            if ix >= len(self.items):
                ix -= self.cycle
                out = out[-self.longx:]
            if isinstance(self.items[ix],tuple):                
                h,k = self.items[ix]
            else:
                h,k = 1,self.items[ix]
            if isinstance(h,ite):
                h0 = h.const
                lr = len(h.r)
                for icx in range(0,lr) :
                    h0 += out[-lr+icx][0]*h.r[icx]
                    h0 += out[-lr+icx][1]*h.s[icx]
            else :
                h0 = h
            if isinstance(k,ite):
                k0 = k.const
                lr = len(k.r)
                for icx in range(0,lr) :
                    k0 += out[-lr+icx][0]*k.r[icx]
                    k0 += out[-lr+icx][1]*k.s[icx]
            else:
                k0 = k
            out.append((h0,k0))
            yield (h0,k0)
    def __repr__(self):
        "representation of kite"
        rstr = "kite("+repr(self.items)+" , "+repr(self.cycle)+")"
        return rstr