# -*- coding: utf-8 -*-
"""

Rx - Recursive xform.

Created on Sun Aug  4 08:49:30 2019.
@author: Kev

"""
from frac import frac
from poly import poly
from series0 import ifun


class rx():
    """Transform vector part of rxform."""

    def __init__(self, xform=None):
        if isinstance(xform, rx):
            self.pt = xform.pt
            self.qt = xform.qt
            self.off = xform.off
        else:
            self.pt = [[1], []]
            self.qt = [[], [1]]
            self.off = 0

    def div_p(self, m):
        """Divide p by x^m."""
        mx = m
        while mx > 0:
            if (self.qt[0] == [] or self.qt[0][0] == 0) and \
                    (self.qt[1] == [] or self.qt[1][0] == 0):
                self.qt[0] = self.qt[0][1:]
                self.qt[1] = self.qt[1][1:]
                self.off += 1
            else:
                self.pt[0] = [0] + self.pt[0]
                self.pt[1] = [0] + self.pt[1]
            mx -= 1

    def div_q(self, m):
        """Divide q by x^m."""
        mx = m
        while mx > 0:
            if (self.pt[0] == [] or self.pt[0][0] == 0) and \
                    (self.pt[1] == [] or self.pt[1][0] == 0):
                self.pt[0] = self.pt[0][1:]
                self.pt[1] = self.pt[1][1:]
                self.off += 1
            else:
                self.qt[0] = [0] + self.qt[0]
                self.qt[1] = [0] + self.qt[1]
            mx -= 1

    def sub_p(self, r):
        """Subtract q*r from p."""
        i = 0
        for pterm in self.qt[0]:
            while len(self.pt[0]) <= i:
                self.pt[0].append(0)
            self.pt[0][i] -= r * pterm
            i += 1
        i = 0
        for qterm in self.qt[1]:
            while len(self.pt[1]) <= i:
                self.pt[1].append(0)
            self.pt[1][i] -= r * qterm
            i += 1

    def sub_q(self, r):
        """Subtract p*r from q."""
        i = 0
        for pterm in self.pt[0]:
            while len(self.qt[0]) <= i:
                self.qt[0].append(0)
            self.qt[0][i] -= r * pterm
            i += 1
        i = 0
        for qterm in self.pt[1]:
            while len(self.qt[1]) <= i:
                self.qt[1].append(0)
            self.qt[1][i] -= r * qterm
            i += 1

    def __repr__(self):
        """Representation."""
        rstr = 'px off = ' + repr(self.off) + ' vec =' + repr(self.pt) + \
            '\nqx off = ' + repr(self.off) + ' vec =' + repr(self.qt)
        return rstr


class rxform():
    """Recursive xform for infinite streams."""

    def __init__(self, p, q=None, rx0=None):
        """Initialize rxform with single function."""
        self.p = p
        if q is None:
            self.q = poly([1])
            self.anx = True
        else:
            self.q = q
            self.anx = False  # q<>1?
        if rx0 is None or not isinstance(rx0, rx):
            self.rx = rx()
        else:
            self.rx = rx0

    def reduce_p(self):
        """Reduce p term."""
        px = poly(ifun(self.eval_p).gen())
        qx = poly(ifun(self.eval_q).gen())
        m = qx.first - px.first
        p0 = px[px.first]
        q0 = qx[qx.first]
        if q0 == 0:
            raise (StopIteration)
        r = frac(p0, q0, True)
#        print px.first,qx.first,m,r
        self.rx.div_q(m)
#        self.rx.norm()
        self.rx.sub_p(r)
        return (m, r)

    def reduce_q(self):
        """Reduce q term."""
        px = poly(ifun(self.eval_p).gen())
        qx = poly(ifun(self.eval_q).gen())
        m = px.first - qx.first
        p0 = px[px.first]
        q0 = qx[qx.first]
        if p0 == 0:
            raise (StopIteration)
        r = frac(q0, p0, True)
#        print px.first,qx.first,m,r
        self.rx.div_p(m)
#        self.rx.norm()
        self.rx.sub_q(r)
        return (m, r)

    def eval_p(self, idx):  # Note: powers (pt) are inverted (based on eval)
        """Evaluate transformed p0."""
        new_p = 0
        no_p = True
        i = idx + self.rx.off
        for r in self.rx.pt[0]:  # p portion from p terms
            if self.p.last is None or i <= self.p.last:
                try:
                    new_p += r * self.p[i]
                    no_p = False
                except IndexError:
                    pass
            i += 1
#        if self.anx:
#            if idx == 0 and len(self.rx.pt[1]) <> 0:
#                new_p += self.rx.pt[1][0]
#                no_p = False
#        else:
        i = idx + self.rx.off
        for r in self.rx.pt[1]:  # p portion from q terms
            if self.q.last is None or i <= self.q.last:
                try:
                    new_p += r * self.q[i]
                    no_p = False
                except IndexError:
                    pass
            i += 1
        if no_p:
            raise (StopIteration)
        else:
            return new_p

    def eval_q(self, idx):  # Note powers (qt) are inverted (based on eval)
        """Evaluate transformed q0."""
        new_q = 0
        no_q = True
        i = idx + self.rx.off
        for r in self.rx.qt[0]:  # q portion from p terms
            if self.p.last is None or i <= self.p.last:
                try:
                    new_q += r * self.p[i]
                    no_q = False
                except IndexError:
                    pass
            i += 1
#        if self.anx:
#            if idx == 0 and len(self.rx.qt[1]) <> 0:
#                new_q += self.rx.qt[1][0]
#                no_q = False
#        else:
        i = idx + self.rx.off
        for r in self.rx.qt[1]:  # q portion from q terms
            if self.q.last is None or i <= self.q.last:
                try:
                    new_q += r * self.q[i]
                    no_q = False
                except IndexError:
                    pass
            i += 1
        if no_q:
            raise (StopIteration)
        else:
            return new_q

    def __repr__(self):
        """Representation of rxform."""
        return repr(self.rx)
