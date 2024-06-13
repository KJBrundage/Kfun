# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:10:32 2020.

Continued fractions

  K is the base for all the continued fraction classes.
  Each K provides a stream of tuples of the form (a x^m, b)
                  where a, b are typically integer (complex integers)
                  and for constants m=0 (i.e these return a series of (a,b) )
        sub classes of K have names of the form
        K I,C,P, S,G S {integer,complex,poly} {simple, generalized}
        KISS are the most common constant (simple continued fractions)
        KIGS are constants represented as generalized continued fractions.
        KPSS are the most common function type where a=1 for all terms.

    K(float) returns KISS
    K()
@author: Kev
"""

from series0 import istream
from cx import cxr, cxi, cxi_r
from digits import digits
from frac import frac
from poly import poly, polyq
from rx import rxform
from ops import mxform, txform, tensor
from approx import approx

from numbers import Number
from inspect import isgenerator
from types import FunctionType

from functools import total_ordering
"""
kite's define streams of K terms as constants + recursive offset values

pi = kite([(4,1),(1,3),(ite(0,[1],[1]),ite(2,[],[1]))],1)

e = kite([1,0,1,1,(1,ite(2,[],[1,0,0]))],3)
"""

global max_chars


def set_max_chars(new_max=60):
    """Set the char displayed for power."""
    global max_chars
    max_chars = new_max


set_max_chars(60)


@total_ordering
class K():
    """Base class for continued fractions."""

    def __init__(self, *argv):
        """Continued fraction instance."""
        narg = len(argv)
        self.n = narg
        self.a = argv[:]
        if narg == 0:
            self.kid = None
        elif narg == 1:
            a0 = argv[0]
            if isinstance(a0, int):           # Ki     integer
                self.kid = Ki(a0)
            elif isinstance(a0, float):       # Kf     float (approximation)
                self.kid = Kf(a0)
            elif isinstance(a0, complex):     # Kc     complex (approximation)
                self.kid = Kc(a0)
            elif isinstance(a0, cxi) \
                    or isinstance(a0, cxr):   # Kcx    complex rational
                self.kid = Kcx(a0)
                # print('CXR type')
            elif isinstance(a0, poly):        # Kp     Polynomial
                self.kid = Kp(a0)
            elif isinstance(a0, polyq):       # .      polyq
                self.kid = Kp(a0.num, a0.den)
            elif isinstance(a0, frac):        # -      frac
                self.kid = Kr(a0.num, a0.den)
            elif isinstance(a0, FunctionType):  # Kg   function
                self.kid = Kg(a0)
            elif isgenerator(a0):             # -      generator
                self.kid = Kg(a0)
            elif isinstance(a0, approx):      # -      approximator
                self.kid = Kg(iter(a0))
            elif isinstance(a0, list):        # Kl     list
                self.kid = Kl(a0)
            else:
                raise TypeError
            self.kiss = self.kid.kiss
        elif narg == 2:
            a0 = argv[0]
            a1 = argv[1]
            # print("2 arg ", type(a0), type(a1))
            if isinstance(a0, int) and isinstance(a1, int):
                self.kid = Kr(a0, a1)       # kr    Rational
            elif isinstance(a0, list) and not isinstance(a0[0], list):
                self.kid = Kl(a0, a1)
            elif isinstance(a0, list) and isinstance(a1, K):
                # Km     Mobius transform (1 var)
                if isinstance(a1.kid, Kt):
                    new_tx = a0 * a1.kid.xform
                    self.kid = Kt(new_tx, a1.kid.x_var, a1.kid.y_var)
                else:
                    self.kid = Km(a0, a1)
            else:
                raise TypeError
            self.kiss = self.kid.kiss
        elif narg == 3:
            tx0 = txform(argv[0])
            x0 = argv[1]
            y0 = argv[2]
            #  Compose Mobius transforms into tesor
            if isinstance(x0.kid, Km):
                mx_x = x0.kid.mx
                my_x = x0.kid.vx
            else:
                mx_x = 1
                my_x = x0

            if isinstance(y0.kid, Km):
                mx_y = y0.kid.mx
                my_y = y0.kid.vx
            else:
                mx_y = 1
                my_y = y0

            if mx_x != 1 or mx_y != 1:
                tx0 = tx0.comp(mx_x, mx_y)

            self.kid = Kt(tx0, my_x, my_y)   # Kt     Tensor transform (2 var)
            self.kiss = self.kid.kiss

    def txmod(self, xvar, yvar):
        """Ubsorb mobius transforms into tensor."""
        print("txmod", self, xvar, yvar)

    def __iter__(self):
        """Return iterator of K yields (m,b) terms."""
        return iter(self.kiss)

    def __pos__(self):
        """Unary + converts to a simple continued fraction."""
        return K([[1, 0], [0, 1]], self)

    def __neg__(self):
        """Negate."""
        return -1*self

    def __add__(self, other):
        """Addition with cf."""
        if isinstance(other, int) or isinstance(other, cxi):
            return K([[1, other], [0, 1]], self)
        elif isinstance(other, frac):
            return K([[other.den, other.num], [0, other.den]], self)
        elif isinstance(other, cxr):
            # print("K add with ", other)
            return K([[1, cxi(other.re, other.im)], [0, other.den]], self)
        elif isinstance(other, K):
            return K([[[0, 1], [1, 0]], [[0, 0], [0, 1]]], self, other)

    def __radd__(self, other):
        """Right addition with cf."""
        return self + other

    def __sub__(self, other):
        """Subtract with cf."""
        if isinstance(other, int):
            return K([[1, -other], [0, 1]], self)
        elif isinstance(other, frac):
            return K([[other.den, -other.num], [0, other.den]], self)
        elif isinstance(other, K):
            return K([[[0, -1], [1, 0]], [[0, 0], [0, 1]]], self, other)

    def __rsub__(self, other):
        """Subtract cf from other."""
        return -1*self + other

    def __mul__(self, other):
        """Multiply Continued fraction."""
        # print("in cf", type(other))
        if isinstance(other, int) or isinstance(other, cxi):
            return K([[other, 0], [0, 1]], self)
        elif isinstance(other, frac):
            return K([[other.num, 0], [0, other.den]], self)
        elif isinstance(other, cxr):
            return K([[cxi(other.re, other.im), 0], [0, other.den]], self)
        elif isinstance(other, K):
            return K([[[1, 0], [0, 0]], [[0, 0], [0, 1]]], self, other)

    def __rmul__(self, other):
        """Right multiplication with cf."""
        return self*other

    def __truediv__(self, other):
        """Division with cf."""
        if isinstance(other, int):
            return K([[1, 0], [0, other]], self)
        elif isinstance(other, frac):
            return K([[other.den, 0], [0, other.num]], self)
        elif isinstance(other, K):
            return K([[[0, 0], [1, 0]], [[0, 1], [0, 0]]], self, other)

    def __rtruediv__(self, other):
        """Right division with cf."""
        if isinstance(other, int):
            return K([[0, other], [1, 0]], self)
        elif isinstance(other, frac):
            return K([[0, other.num], [other.den, 0]], self)
        elif isinstance(other, K):
            return K([[[0, 1], [0, 0]], [[0, 0], [1, 0]]], self, other)

    def __pow__(self, other):
        """Continued fraction to power."""
        if isinstance(other, int):
            if other == 0:
                return K([1])
            elif other == 1:
                return self
            elif other > 1:
                return self*(self**(other-1))
            else:
                return (1/self)**(-other)

    def __repr__(self):
        """Represent K form."""
        return repr(self.kid)

    def __str__(self):
        """Generate string of K's iterator."""
        k_string = ""
        pc = 0
        kit = iter(self.kid)
        try:
            xm, b = next(kit)
        except StopIteration:
            return '*0'
        fini = False
        if b != 0:
            k_string += repr(b)
        while len(k_string) < max_chars and not fini:
            try:
                xm, b = next(kit)
            except StopIteration:
                if k_string == '':
                    k_string = '0'
                fini = True
                # print( 'terminate output')
                break
                # raise (StopIteration)

            if (k_string != '' and
                ((isinstance(xm, int) and xm >= 0)
                    or (isinstance(xm, frac) and xm.num >= 0)
                    or (isinstance(xm, cxi) and xm.re >= 0)
                    or (isinstance(xm, poly) and xm.c0 >= 0))):
                k_string += '+'
            if (isinstance(b, int) or isinstance(b, frac)) and b == 0:
                k_string += repr(xm) + '/('
                pc += 1
            else:
                k_string += repr(xm) + '/(' + repr(b)
                pc += 1
        if not fini:
            k_string += "+ ..."
        else:
            if k_string[-3:] == '/(1':
                k_string = k_string[:-3]
                pc -= 1
            for _ in range(pc):
                k_string += ')'
        return k_string

    def __call__(self, var=None):
        """Approximator for K variable.

         When the continued fraction is a function
        this returns the Pade approximations.
        """
        if var is None:
            # print("K without arg")
            return approx(self)
        else:
            # print("K with arg")
            # return K(iter(approx(self, var)))
            return approx(self, var)

    def __getitem__(self, idx):
        """Use the kiss member to index into it's istream."""
        return self.kiss[idx]

    def __lt__(self, other):
        """Is variable less than other."""
        sflag = True
        maxit = 1000
        scnt = 0
        if isinstance(other, K):
            otherk = other
        else:
            otherk = K(other)
        sit = iter(self)
        oit = iter(otherk)
        while scnt < maxit:
            try:
                sterm = next(sit)[1]
            except StopIteration:
                sterm = None
                sit = maxit
            try:
                oterm = next(oit)[1]
            except StopIteration:
                oterm = None
                sit = maxit
            if sterm is None and oterm is None:
                return False
            if oterm is None or sterm < oterm:
                if sflag:
                    return True
                else:
                    return False
            if sterm is None or sterm > oterm:
                if sflag:
                    return False
                else:
                    return True
            sflag = not sflag
            scnt += 1
        return False

    def __eq__(self, other):
        """Is variable equal to other."""
        maxit = 1000
        scnt = 0
        if isinstance(other, K):
            otherk = other
        else:
            otherk = K(other)
        sit = iter(self)
        oit = iter(otherk)
        while scnt < maxit:
            try:
                sterm = next(sit)[1]
            except StopIteration:
                sterm = 0
                sit = maxit
            try:
                oterm = next(oit)[1]
            except StopIteration:
                oterm = 0
                sit = maxit
            if sterm != oterm:
                return False
            scnt += 1
        return True

    # def __gt__(self, other):
    #     """Is first variable greater than second."""
    #     if isinstance(other, K):
    #         return not self < other
    #     else:
    #         return not self < K(other)

        # sgen = iter(self)
        # h0, h1 = 0, 1
        # k0, k1 = 1, 0
        # done = False
        # if var is None:
        #     while not done:
        #         try:
        #             a, b = next(sgen)
        #         except StopIteration:
        #             a, b = 0, 1
        #             done = True
        #         h0, h1 = h1, h0*a + h1*b
        #         k0, k1 = k1, k0*a + k1*b
        #         yield frac(h1, k1)
        # elif isinstance(var, int):
        #     while not done:
        #         for _ in range(2*var):
        #             try:
        #                 a, b = next(sgen)
        #             except StopIteration:
        #                 a, b = 0, 1
        #                 done = True
        #             h0, h1 = h1, h0*a + h1*b
        #             k0, k1 = k1, k0*a + k1*b
        #         yield frac(h1, k1)
        #         var = 1
        # elif isinstance(var, frac):
        #     return approx(pade())
        # my_kiss = []
        # my_app = approx(self(0))
        # f0 = +my_app  # skip first pade
        # print("skip first pade", f0)
        # while not done:
        #     f0 = +my_app  # next pade approximation
        #     print("pade approximation", f0)
        #     yield approx(f0(var))
        # my_ix = 0

        # itx0 = f0(var)
        # itx1 = f1(var)
        # t0 = next(itx0)
        # t1 = next(itx1)
        # while t0 == t1:
        #     if my_ix == len(my_kiss):
        #         yield((1, t0))
        #         my_kiss.append((1, t0))
        #     elif t0 != my_kiss[my_ix]:
        #         raise ArithmeticError

        #     my_ix += 1
        #     t0 = next(itx0)
        #     t1 = next(itx1)

        # else:
        #     print('comming soon')
    def sig(self, terms):
        """Return the first terms of the signature of a continued fraction."""
        myit = iter(+self)
        mysig = []
        for _ in range(terms):
            try:
                term = next(myit)
                if term[0] == 1:
                    mysig.append(term[1])
                else:
                    mysig.append(term)
            except StopIteration:
                break
        return mysig

    def out(self, ndigits, base=10):
        """Return string of ndigits."""
        dgen = digits(self, base)
        out_str = next(dgen)
        if isinstance(out_str, tuple):
            out_real, out_img = out_str
            for _ in range(ndigits-1):
                out_str = dgen.send(False)
                out_real += out_str[0]
                out_img += out_str[1]
            out_str = dgen.send(True)
            out_real += out_str[0]
            out_img += out_str[1]
            return out_real + "+" + out_img + cxi_r
        else:
            for _ in range(ndigits-1):
                out_str += dgen.send(False)
            out_str += dgen.send(True)
        return rup_str(out_str, base)


def ap_k(kvar):
    """Approximator function for k variable."""
    sgen = iter(kvar)
    h0, h1 = 0, 1
    k0, k1 = 1, 0
    done = False
    while not done:
        try:
            a, b = next(sgen)
        except StopIteration:
            a, b = 0, 1
            done = True
#        print("in ap_k", h0, h1, a, b, type(a), type(b))
        h0, h1 = h1, h0*a + h1*b
        k0, k1 = k1, k0*a + k1*b
        yield frac(h1, k1)


def pade(fun, var):
    """Pade approximation generator for a cf function."""
    my_app = approx(fun)
    done = False
    while not done:
        yield approx(+my_app(var))


def rup_str(out_str, base):
    """Correct string for round up."""
    if out_str[-1] == '*':  # rounded up 9 (base -1), must carry
        if base <= 10:
            hi_c = chr(47 + base)  # chr for highest digit in base (base-1)
        else:
            hi_c = chr(54 + base)
        cidx = -2
        out_str = out_str[:-1] + '0'
        while out_str[cidx] == hi_c or out_str[cidx] == '.':
            if out_str[cidx] == '.':
                cidx -= 1
            else:
                out_str = out_str[:cidx] + '0' + out_str[cidx + 1:]
                cidx -= 1
            if -cidx > len(out_str):
                out_str = '0' + out_str
        c_c = chr(ord(out_str[cidx]) + 1)
        out_str = out_str[:cidx] + c_c + out_str[cidx + 1:]
    return out_str


"""
  these are K's kid's
  each kid has an istream that yields (a,b) values (kiss)
"""


class Ki(K):
    """Continued fraction integer."""

    def __init__(self, ival):
        """Initialize KI."""
        self.ival = ival
        self.kiss = istream(i2k(ival))

    def __iter__(self):
        """Iterate."""
        return iter(self.kiss)

    def __repr__(self):
        """Representation of KI."""
        return "K(" + str(self.ival) + ")"


def i2k(i):
    """Integer to k iterator."""
    yield (1, i)


class Kr(K):
    """Continued fraction fraction."""

    def __init__(self, num, den):
        """Initialize KI."""
        f0 = frac(num, den)
        self.num = f0.num
        self.den = f0.den
        self.kiss = istream(r2k(self.num, self.den))

    def __iter__(self):
        """Iterate."""
        return iter(self.kiss)

    def __repr__(self):
        """Representation of KI."""
        return "K(" + str(self.num) + "," + str(self.den) + ")"


def r2k(num, den):
    """Rational to k iterator."""
    while (den != 0):
        t0 = num // den
        yield (1, t0)
        num, den = den, num - t0*den


class Kg(K):
    """Generator type."""

    def __init__(self, igen):
        """Initialize Kg."""
        self.gen = igen
        self.kiss = istream(igen)

    def __iter__(self):
        """Iterate Kg."""
        return iter(self.kiss)

    def __repr__(self):
        """Represent Kg."""
        return "K("+str(self.gen)+")"


class Kf(K):
    """Continued fraction float."""

    def __init__(self, fval):
        """Initialize KI."""
        self.fval = fval
        self.kiss = istream(f2k(fval))

    def __iter__(self):
        """Return iterator."""
        return iter(self.kiss)

    def __repr__(self):
        """Return representation of KI."""
        return "K("+str(self.fval) + ")"


def f2k(f):
    """Float to k iterator."""
    if f < 0:
        n = int(f) - 1
    else:
        n = int(f)
    yield (1, n)
    r = f - n
    while r > 0:
        r = 1/r
        if r > 10000000000:  # ???
            break
        yield (1, int(r))
        r -= int(r)


class Kc(K):
    """Continued fraction complex type."""

    def __init__(self, cval):
        """Initialize KI."""
        self.cval = cval
        self.kiss = istream(c2k(cval))

    def __iter__(self):
        """Iterate."""
        return iter(self.kiss)

    def __repr__(self):
        """Represent KI."""
        return "K(" + str(self.cval) + ")"


def c2k(c):
    """Complex to K."""
    rc = int((c.real+.5) // 1)
    ic = int((c.imag+.5) // 1)
    yield (1, cxi(rc, ic))
    c -= rc + (0+1j)*ic
    while c != 0:
        c = 1/c
#        print(c)
        rc = int((c.real+.5) // 1)
        ic = int((c.imag+.5) // 1)
        yield (1, cxi(rc, ic))
        c -= rc + (0+1j)*ic


class Kcx(K):
    """Continued fraction (CX) complex."""

    def __init__(self, cval):
        """Initialize KI."""
        self.cval = cval
        self.kiss = istream(cx2k(cval))

    def __iter__(self):
        """Iterate."""
        return iter(self.kiss)

    def __repr__(self):
        """Representation of Kcx."""
        return "K(" + str(self.cval) + ")"


def cx2k(c):
    """Complex to k interator."""
    yield (1, cxi(c))
    r = c - cxi(c)
    while r != 0:
        r = 1/r
        # print("cxck", r, cxi(r))
        yield (1, cxi(r))
        r -= cxi(r)


class Kp(K):
    """Continued fraction of poly or polyq."""

    def __init__(self, p, q=poly([1])):
        self.p = p
        self.q = q
        self.kiss = istream(p2k(p, q))

    def __iter__(self):
        """Iterate."""
        return iter(self.kiss)

    def __repr__(self):
        """Representation for Kp."""
        r_x = "K("+repr(self.p)
        if self.q != poly([1]):
            r_x += ","+repr(self.q)
        r_x += ")"
        return r_x


def p2k(p, q=poly([1])):
    """Polynomial to K terms."""
    my_rx = rxform(p, q)
    done = False
    if p.first == q.first:
        xm, b = my_rx.reduce_p()
        term = (poly([1], xm), b)
        yield term
    else:
        yield (poly([1]), 0)
    if p.first < q.first:
        print("skipq")
        yield (poly([1]), 0)
        skipq = True
    else:
        skipq = False
    while not done:
        if not skipq:
            xm, b = my_rx.reduce_q()
            term = (poly([1], xm), b)
            yield term
        else:
            skipq = False
        xm, b = my_rx.reduce_p()
        term = (poly([1], xm), b)
        yield term


class Kl(K):
    """Continued fraction from list (+offset)."""

    def __init__(self, terms, offset=[]):
        """Instance of K list."""
        if isinstance(offset, int):
            offset = [0]*offset
        my_terms = []

        for termx in terms:
            if isinstance(termx, int) or isinstance(termx, cxi):
                termx = (1, termx)
            my_terms.append(termx)

        len_pre = len(terms) - len(offset)
        self.pre = my_terms[0:len_pre]
        self.cycle = my_terms[len_pre:]
        self.offset = offset
        self.kiss = istream(iter(self))

    def __iter__(self):
        """Iterate Kl."""
        for term in self.pre + self.cycle:
            yield term
        if self.cycle == []:
            return   # ???
        my_c = self.cycle[:]
        my_x = self.offset[:]
        while True:
            term = my_c[0]
            xvec = my_x[0]
            term = (term[0], term[1]+xvec)
            yield term

            my_c = my_c[1:]
            my_c.append(term)
            my_x = my_x[1:]
            my_x.append(xvec)

    def __repr__(self):
        """Representation for Kl."""
        len_c = len(self.offset)
        terms = self.pre + self.cycle
        t2 = []
        for termx in terms:
            a, b = termx
            if a == 1:
                t2.append(b)
            else:
                t2.append((a, b))
        if self.offset == [0]*len_c:
            return "K("+repr(t2) + "," + repr(len_c)+")"
        else:
            return "K("+repr(t2) + "," + repr(self.offset)+")"


class Km(K):
    """Continued fraction mobius."""

    def __init__(self, xform, var):
        # print("init",var,xform)
        if isinstance(var.kid, Km):
            # print("compose xforms")
            xform = mxform(xform) * var.kid.mx
            var = var.kid.vx
        # elif isinstance(var.kid, Kt):
        #     if xform[1][0] == 0:  #  affine transform

        self.vx = var
        self.mx = mxform(xform)
        self.kiss = istream(iter(self))

    def __iter__(self):
        """Iterate k terms."""
        kiss = iter(self. vx.kid)
        my_rx = mxform(self.mx)
        first = True
        while True:
            cx, c0 = my_rx.p
            dx, d0 = my_rx.q
            if d0 == 0:
                r0 = None
            else:
                r0 = c0 // d0
            if dx == 0:
                r1 = None
            else:
                r1 = cx // dx
            if r0 is None and r1 is None:
                # print("Km:",my_rx)
                return
            while r0 is None or r0 != r1:
                try:
                    xmb = next(kiss)
#                    print('imgen in_x', xmb)
#                    print('my_rx in =',my_rx,xmb)
                    my_rx.in_x(xmb)
#                    print('my_rx out =',my_rx)
                except StopIteration:
                    my_rx.in_inf()
#                    print('my_rx inf =',my_rx)
                cx, c0 = my_rx.p
                dx, d0 = my_rx.q
                if dx == 0:
                    r0, r1 = r1, None
                else:
                    r0, r1 = r1, cx // dx
#            print("mobius iter",my_rx,type(r0),r0,r1)
            if isinstance(r0, Number) or isinstance(r0, cxi):
                yield (1, r0)
                # print("ops number")
                my_rx.sub_p(r0)
                my_rx.swap()
            elif isinstance(r0, poly):
                # print("mobius yield",r0)
                if first:
                    mx = r0.first
                    b = r0[mx]
#                    print("first",b,mx)
                    if mx != 0:
                        yield (poly([1]), 0)
                        if mx < 0:
                            yield (poly([1]), 0)
                            my_rx.swap()
#                        my_rx.swap()
#                        print("after swap 0",my_rx)
                    else:
                        yield (poly([1]), b)
                        my_rx.sub_p(b)
                    my_rx.swap()
                    first = False
                else:
                    mx = r0.first
                    b = r0[mx]
#                    print("b,mx=",b,mx)
#                    print("before yield",my_rx)
                    yield (poly([1], -mx), b)  # K form
                    my_rx.mul_p(poly([1], -mx))
                    my_rx.sub_p(b)
                    my_rx.swap()
#                    print("after swap",my_rx)
        return

    def __repr__(self):
        """Representation of mobius."""
        mstr = "K("+repr(self.mx)+", "+repr(self.vx)+")"
        return mstr


class Kt():
    """Continued fraction tensor (2 variable operator)."""

    def __init__(self, tensor, x_var, y_var):
        """Initialize tensor."""
        self.x_var = x_var  # ??? +
        self.y_var = y_var  # ??? +
        self.xform = txform(tensor)
        self.kiss = istream(iter(self))

    def __iter__(self):
        """Iterate tensor."""
        myrx = txform(self.xform)
        mytx = tensor(myrx, self.x_var, self.y_var)
        done = False
        ecnt = 0
        close = False
        first = True
        while not done:
            #            print('myrx:',myrx,close)
            r0 = mytx.rat(close)
            # rl = mytx.rat2(close)
            # if isinstance(rl, list):
            #     r0 = rl[0]
            # else:
            #     r0 = rl
            # print("Kt_iter", close, r0)
            if isinstance(r0, poly):
                mx = r0.first
                b = r0[mx]
                r0 = (poly([1], -mx), b)
            if r0 == 0:
                return
            if close and not isinstance(r0, tuple):
                if mytx.inf():
                    return
            if isinstance(r0, tuple):
                ecnt = 0
                if first and r0[0] != 1:
                    r0 = (poly([1]), 0)
                    yield r0
#                    mytx.out((poly([1]),0))
                    first = False
                    # print('before out0',mytx)
                    mytx.out(r0)
                    # print('after out0',mytx)
                else:
                    yield r0
                    # print('before out',mytx)
                    mytx.out(r0)
                    # print('after out',mytx)
            else:
                if r0 == 1:
                    mytx.in_r()
                elif r0 == 2 or r0 == 3:
                    mytx.in_y()
                elif r0 == 4 or r0 == 5:
                    mytx.in_x()
                elif r0 == 6 or r0 == 7:
                    # print("close", ecnt, repr(mytx))
                    mytx.in_x()
                    mytx.in_y()
                    mytx.in_x()
                    mytx.in_y()
                    ecnt += 1
                    # mytx.in_r()
                    # mytx.in_r()
                    # mytx.in_r()
                    if ecnt > 4:
                        close = True
                elif r0 == 15:
                    done = True
                else:
                    # print("err=", r0, repr(mytx))  # zero dxy in denominator
                    mytx.in_x()
                    mytx.in_y()
                    mytx.in_x()
                    mytx.in_y()
            if mytx.inf():
                done = True

    def __iter2_(self):
        """Iterate tensor."""
        myrx = txform(self.xform)
        mytx = tensor(myrx, self.x_var, self.y_var)
        done = False
        ecnt = 0
        close = False
        first = True
        while not done:
            #            print('myrx:',myrx,close)
            # r0 = mytx.rat(close)
            rl = mytx.rat2(close)
            if isinstance(rl, list):
                r0 = rl[0]
                r1 = rl[1]
            else:
                r0 = rl
                r1 = None
            # print("Kt_iter", close, r0, r1)
            if isinstance(r0, poly):
                mx = r0.first
                b = r0[mx]
                r0 = (poly([1], -mx), b)
            if r0 == 0:
                return
            if close and not isinstance(r0, tuple):
                if mytx.inf():
                    return
            if isinstance(r0, tuple):
                ecnt = 0
                if first and r0[0] != 1:
                    r0 = (poly([1]), 0)
                    yield r0
#                    mytx.out((poly([1]),0))
                    first = False
                    # print('before out0',mytx)
                    mytx.out(r0)
                    # print('after out0',mytx)
                else:
                    yield r0
                    # print('before out',mytx)
                    mytx.out(r0)
                    # print('after out',mytx)
            else:
                if r0 == 1:
                    mytx.in_r()
                elif r0 == 2 or r0 == 3:
                    if r1 == 1:
                        mytx.in_y()
                    else:
                        mytx.in_y()
                        mytx.in_y()
                elif r0 == 4 or r0 == 5:
                    if r1 == 1:
                        mytx.in_x()
                    else:
                        mytx.in_x()
                        mytx.in_x()
                elif r0 == 6 or r0 == 7:
                    # print("close", ecnt, repr(mytx))
                    mytx.in_x()
                    mytx.in_y()
                    mytx.in_x()
                    mytx.in_y()
                    ecnt += 1
                    # mytx.in_r()
                    # mytx.in_r()
                    # mytx.in_r()
                    if ecnt > 4:
                        close = True
                elif r0 == 15:
                    done = True
                else:
                    # print("err=", r0, repr(mytx))  # zero dxy in denominator
                    mytx.in_x()
                    mytx.in_y()
                    mytx.in_x()
                    mytx.in_y()
            if mytx.inf():
                done = True
#         xg   = iter(self.x_var)
#         yg   = iter(self.y_var)
#         xytog= True
#         done = False
#         while( not done):
#             r0 = myrx.rat()
#             if isinstance(r0, tuple):
#                 yield r0
#                 myrx.outUsage
#             else :
#                 if r0 == 1:
#                     if xytog :
#                         myrx.in_x(next(xg))
#                     else :
#                         myrx.in_y(next(yg))
#                     xytog = not xytog
#                 elif r0 == 2 or r0 == 3 :
#                     myrx.in_y(next(yg))
#                 elif r0 == 4 or r0 == 5 :
#                     myrx.in_x(next(xg))
#                 elif r0 == 6 or r0 == 7:
# #                    print("close")
#                     if xytog :
#                         myrx.in_x(next(xg))
#                     else :
#                         myrx.in_y(next(yg))
#                     xytog = not xytog
#                 elif r0 == 15 :
#                     return
#                 else:
#                     if xytog :
#                         myrx.in_x(next(xg))
#                         myrx.in_x(next(xg))
#                     else :
#                         myrx.in_y(next(yg))
#                         myrx.in_y(next(yg))
#                     xytog = not xytog

    def __repr__(self):
        """Representation for Kt."""
        rstr = "K(" + repr(self.xform) + ", " + \
            repr(self.x_var) + ", " + repr(self.y_var) + ")"
        return rstr


def nfgen(kvar):
    """Iterate normal form (b0+a1/(1+a2/(1+...)."""
    kit = iter(kvar)
    a0, b0 = next(kit)
    yield (a0, b0)
    b0 = 1
    while True:
        try:
            a1, b1 = next(kit)
            if isinstance(a1, poly):
                deg = a1.first
                ac = a1[deg]
                yield (poly([frac(ac, b0*b1, True)], deg), 1)
            else:
                yield (frac(a1, b1*b0, True), 1)
            b0 = b1
        except StopIteration:
            return


def nform(kvar):
    """N-form of continued fraction."""
    return K(nfgen(kvar))


def kfgen(kvar, bp=False, xz=False):
    """Iterate k form (b0 + 1/(b1+1/(b2+...).

    bp = True => bx is positive.
    xz = True => x^n with bx terms, o'wise with 1.
    """
    kit = iter(kvar)
    a0, b0 = next(kit)
    yield (a0, b0)
    cx = 1
    if isinstance(a0, poly):
        dx = a0.first
        a0 = a0[dx]
    while True:
        try:
            a1, b1 = next(kit)
            if isinstance(a1, poly):
                deg = a1.first
                ac = a1[deg]
                if xz:
                    dx = deg - dx
                    deg = 0

                if a0*b1 < 0 and bp:
                    if xz:
                        yield (cx*poly([-1]),
                               frac(-a0*b1, poly([ac], dx)))
                    else:
                        yield (cx*poly([-1], deg), frac(-a0*b1, ac, True))
                    cx = -1
                else:
                    if xz:
                        yield (cx*poly([1]),
                               frac(a0*b1, poly([ac], dx)))
                    else:
                        yield (cx*poly([1], deg), frac(a0*b1, ac, True))
                    cx = 1
                a0 = ac
            else:
                yield (1, frac(a0*b1, a1))
                a0 = a1
        except StopIteration:
            return


def kform(kvar, bp=False, xz=False):
    """K-form of continued fraction.

    bp = True => bx is positive.
    xz = True => x^n with bx terms, o'wise with 1.
    """
    return K(kfgen(kvar, bp, xz))
