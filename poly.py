# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:56:23 2019.

@author: Kev
"""
import numbers
# import matplotlib.pyplot as plt
# import numpy as np

# def PolyCoefficients(x,coeffs,start=0):
#     """ Returns a polynomial for ``x`` values for the ``coeffs`` provided.

#     The coefficients must be in ascending order (``x**0`` to ``x**o``).
#     """
#     o = len(coeffs)
#     print(f'# This is a polynomial of order {ord}.',o,coeffs,start)
#     y = 0
#     for i in range(o):
#         try:
#             y += coeffs[i]*x**(i-start)
#         except:
#             pass
#     return y

from frac import frac
from cx import cxi, cxr
from approx import approx

from series0 import istream

global pow_r


def set_pow_r(char='^'):
    """Set the char displayed for power."""
    global pow_r
    pow_r = char


pow_r = '^'
set_pow_r()


def p_add(p1, p2):
    """Add two poly's."""
    idx = min(p1.first, p2.first)  # min(p1.first,p2.first)
#    print("p_add:", p1[idx], p2[idx], type(p1[idx]), type(p2[idx]))
    psum = p1[idx] + p2[idx]
    while p1.last is None or p2.last is None or idx <= max(p1.last, p2.last):
        yield psum
        idx += 1
        psum = p1[idx] + p2[idx]


def p_mul(p1, p2):
    """Multiply two polynomials."""
    idx = p1.first+p2.first  # p1.first + p2.first
    while p1.last is None or p2.last is None or idx <= (p1.last + p2.last):
        tsum = 0
        last_idx = idx - p2.first
        if p1.last is not None:
            last_idx = min(last_idx, p1.last)
        for p_idx in range(p1.first, last_idx + 1):
            #            print p_idx, idx-p_idx, p1[p_idx], p2[idx-p_idx]
            tsum += p1[p_idx] * p2[idx - p_idx]
#            print (idx,tsum, p_idx, idx-p_idx)
        yield tsum
        idx += 1


# def p_divx(p1, ix):
#     """Divide polynomial by x^px."""
#     idx = 0   # ??? What about Laurent functions?
#     while idx < px:
#         if p1[idx] != 0:
#             raise TypeError
#         idx += 1
#     while p1.last is None or idx <= p1.last:
#         yield p1[idx]
#         idx += 1


def p_mulc(p1, c):
    """Multiply a polynomial by a constant (note:can cast all coef as frac)."""
    idx = p1.first  # p1.first
    pmul = p1[idx] * c
    while p1.last is None or idx <= p1.last:
        yield pmul
        idx += 1
        pmul = p1[idx] * c


def p_differentiate(p):
    """Differentiate polynomial stream."""
    idx = p.first  # ignore constant term
    while p.last is None or idx <= p.last:
        yield idx * p[idx]
        idx += 1


def p_integrate(p, c0):
    """Integrate polynomial stream."""
    yield c0
    idx = p.first
    while p.last is None or idx <= p.last:
        if idx == -1:
            if p[idx] == 0:
                yield 0
            else:
                raise TypeError
        else:
            yield p[idx]/(idx+1)
        idx += 1


class p_iter():
    """Iterator class for poly (account for offset)."""

    def __init__(self, p):
        """Create instance of p_iter."""
        self.poly = p
        self.idx = p.first

    def __next__(self):
        """Next function."""
        while self.poly.last is None or self.idx <= self.poly.last:
            result = self.poly[self.idx]
            if result != 0:
                self.idx += 1
                return pe(result, self.idx - 1)
            else:
                self.idx += 1


class poly():
    """Polynomial class."""

    def __init__(self, p_stream, my_first=0):
        """Poly instance."""
        if isinstance(p_stream, poly):
            self.poff = p_stream.poff + my_first
            self.first = p_stream.first + my_first
            if p_stream.last is not None:
                self.last = p_stream.last + my_first
            else:
                self.last = None
            self.coef = p_stream.coef
            self.c0 = p_stream.c0
#            print ('poly type>?')
        else:
            self.first = my_first  # first is first non-zero term > poff
            self.poff = my_first    # poff is offset into istreams start
            self.coef = istream(p_stream)
            if isinstance(p_stream, list):
                self.last = self.first + len(p_stream) - 1
                while self[self.last] == 0 and self.last > self.first:
                    self.last -= 1
            else:
                self.last = None
#            print('istream')
        self.c0 = self[self.first]
#        print('start',self.first,first_c,self.last)
        while self.c0 == 0 and (self.last is None or self.first < self.last):
            self.first += 1
            self.c0 = self[self.first]
        if self.c0 == 0:
            self.first, self.last = 0, 0
#        self.sig = self[0:10]  # ???
#            print('looking for first',self.first,first_c)
#            try:
#                self.coef.drop(1)
#            except IndexError:
#                self.first, self.last = 0,0
#            first_c = self[self.first]
#        if first_c == 0:
#            self.first = self.last
#        c_cnt = self.first + 1
#        print('still',self.first,self.last)
#        check = 0
#        while self.last is None and check == 0: # to infinity and beond?
#            check = self[c_cnt]
#            c_cnt += 1
#        print('still',self.first,self.last)
#        while self.last is not None and self.first < self.last \
#            and self[self.last] == 0:
#                self.last -= 1
#        print('end',self.first,self.last)

    def zero(self):
        """Is polynomial zero."""
        if self[self.first] == 0:
            return True
        else:
            return False

    def __eq__(self, other):
        """Compare."""
        if isinstance(other, poly):
            if self.first != other.first:
                return False
            cnt = self.first
            while self.last is None or other.last is None:
                if self[cnt] == other[cnt]:
                    cnt += 1
                else:
                    return False
            if self.last is not None and other.last is not None and \
               self.last == other.last:
                for i in range(self.first, self.last+1):
                    #                    print i, self[i],other[i]
                    if self[i] != other[i]:
                        return False
            else:
                return False
            return True
        else:
            return self == poly([other])

    def __abs__(self):
        """Absolute value - lowest order term >= 0."""
        if self[self.first] < 0:
            return -self
        else:
            return self

    def __add__(self, other):
        """Add two polys."""
        if isinstance(other, poly):
            return poly(p_add(self, other), min(self.first, other.first))
        elif isinstance(other, frac):
            #           print("poly add frac #####")
            return self + poly([other])
        else:
            return self + poly(other)

    def __radd__(self, other):
        """Add from right."""
        #        print 'radd'
        return self + other

    def __neg__(self):
        """negate."""
        return -1*self

    def __sub__(self, other):
        """subtract."""
        return self + -1 * poly(other)

    def __rsub__(self, other):
        """Sub with non poly."""
        return poly(other) - self

    def __mul__(self, other):
        """multiply."""
        if isinstance(other, poly):
            return poly(p_mul(self, other), self.first + other.first)
        else:
            if other == 0:
                return poly([0])
            elif isinstance(other, numbers.Number) or isinstance(other, frac):
                return poly(p_mulc(self, other), self.first)
            else:
                return self * poly([other])

    def __rmul__(self, other):
        """mult."""
        return self * other

    def __truediv__(self, other):
        """Divide poly by constant => poly or poly => kpoly."""
        if isinstance(other, poly):
            if other.first == other.last:
                return poly(self, -other.first) / other.c0
            return polyq(self, other).px
        elif isinstance(other, numbers.Integral) or isinstance(self, frac):
            return poly(p_mulc(self, frac(1, other)), self.first)
        else:
            return self/poly(other)

    def __floordiv__(self, other):
        """Floor div for poly."""
        if isinstance(other, numbers.Number):
            my_den = poly(other)
        else:
            my_den = other
        px = self.first - my_den.first
        if my_den == 0:
            raise ZeroDivisionError
#        rx = self[self.first]/other[other.first]
        return poly([frac(self.c0, my_den.c0)], px)
#        return (px,frac(self[self.first],other[other.first]))

    def __rfloordiv__(self, other):
        """Floor div for poly with first arg non-poly."""
        return poly(other) // self

    def __mod__(self, other):
        """Ppoly modulo."""
        return self - (self//other)*other
#    def __truediv__(self,other) :
#        return __div__(self,other)

    def __rdiv__(self, other):
        """Divide with non-poly by poly."""
        return poly(other) / self

    def __pow__(self, other):
        """Polynomial to a power."""
        if isinstance(other, int):
            result = 1
            for _ in range(other):
                result *= self
            return result

    def __getitem__(self, idx):
        """Index into stream."""
        #            print idx, self.first
        if isinstance(idx, slice):
            result = []
            my_start = 0
            if idx.start is not None:
                my_start = idx.start
            for i in range(my_start, idx.stop):
                try:
                    if i - self.poff >= 0:
                        result.append(self.coef[i-self.poff])
                    else:
                        result.append(0)
                except IndexError:  # ??? is exception handled by getitem below
                    self.last = self.first + self.coef.end
                    while self[self.last] == 0 and self.last > self.first:
                        self.last -= 1
                    result.append(0)
            return result
        else:
            coef_idx = idx - self.poff
            if coef_idx >= 0 and \
                    (self.coef.end is None or coef_idx <= self.coef.end):
                try:
                    cv = self.coef[coef_idx]
                    if cv is None:
                        self.last = self.first + self.coef.end
                        cv = 0
                    return cv
                except IndexError:
                    self.last = self.first + self.coef.end
                    while self[self.last] == 0 and self.last > self.first:
                        self.last -= 1
                    return 0
            else:
                return 0

    def __call__(self, x_var=None):
        """Evaluate power series, if x_var == None return generator."""
        #        print 'call poly ',self,' with ',x_var
        if x_var is None:
            # print("poly call changed")
            return approx(self)
        # else:
            # print("comming soon")
        #        print y_out, y_old
#        if isinstance(x_var, numbers.Integral) or isinstance(x_var, frac):
#            p_sum = 0
#            x_pow = 1
#            for _ in range(0, self.first):
#                x_pow *= x_var
#            pidx = self.first
#            while self.last is None or pidx <= self.last:
#                #            for pidx in range(self.first,self.last+1) :
#                p_sum += self[pidx] * x_pow
#                pidx += 1
#                x_pow *= x_var
#            return p_sum
        if isinstance(x_var, cxi) or isinstance(x_var, cxr):
            p_sum = 0
            x_pow = cxr(1, 0, 1)
            if self.first < 0:
                for _ in range(self.first, 0):
                    x_pow /= x_var
            elif self.first > 0:
                for _ in range(0, self.first):
                    x_pow *= x_var
            pidx = self.first
            while self.last is None or pidx <= self.last:
                if self[pidx] is not None and self[pidx] != 0:  # ??? infinite
                    p_sum += self[pidx] * x_pow
                    if self.last is None:
                        return p_sum
                pidx += 1
                x_pow *= x_var
            #            print 'poly_out = ',p_sum
            return p_sum
#        elif isinstance(x_var, k) and isinstance(x_var.kid, kcon):
#            print 'i am k', x_var
#            p_sum = 0
#            x_pow = 1
#            for _ in range(0, self.first):
#                x_pow *= x_var
#            pidx = self.first
#            while self.last is None or pidx <= self.last:
#                #            for pidx in range(self.first,self.last+1) :
#                p_sum += self[pidx] * x_pow
#                pidx += 1
#                x_pow *= x_var
#            return p_sum
        elif isinstance(x_var, complex):
            y_gen = self.conv(x_var)
            y_out = next(y_gen)
            y_old = y_out - 1
            while y_out.real != y_old.real or y_out.imag != y_old.imag:
                y_old = y_out
                try:
                    y_out = next(y_gen)  #
                except StopIteration:
                    pass
            return y_out
        elif isinstance(x_var, frac):
            p_sum = 0
            x_pow = frac(1, 1)
            if self.first < 0:
                for _ in range(self.first, 0):
                    x_pow /= x_var
            elif self.first > 0:
                for _ in range(0, self.first):
                    x_pow *= x_var
            pidx = self.first
            while self.last is None or pidx <= self.last:
                p_sum += self[pidx] * x_pow
                pidx += 1
                x_pow *= x_var
            return p_sum
        else:
            # print('estimate')  # get convergences until at machine precision
            y_gen = self.conv(x_var)
            y_out = next(y_gen)
            y_old = y_out - 1
            while (float(y_out) != float(y_old)):  # ??? **** change?
                y_old = y_out
                try:
                    y_out = next(y_gen)  #
                except StopIteration:
                    pass
                # except OverflowError :
                #     # y_out = y_old
                #     print('overflowError',x_var)
            if isinstance(y_out, frac) and isinstance(y_out.num, float):
                # print('floater')
                return y_out.num/y_out.den
            else:
                return y_out

    def __iter__(self):
        """Iterate."""
        return p_iter(self)

    def __repr__(self):
        """Display poly type."""
        global pow_r
        p_string = ""
        i_ptr = self.first
        while len(p_string) < 70 and (self.last is None or i_ptr <= self.last):
            #            x_val = self[i_ptr+1] #read ahead to find last
            x_val = self[i_ptr]
            while x_val == 0 and (self.last is None or i_ptr <= self.last):
                i_ptr += 1
                x_val = self[i_ptr]
            #            if self.last is not None and i_ptr > self.last :
            #                pass
            #            else :
            if i_ptr == 0:
                pow_string = ''
            elif i_ptr == 1:
                pow_string = 'x'
            else:
                if i_ptr < 0:
                    pow_string = 'x' + pow_r + '(' + str(i_ptr) + ')'
                else:
                    pow_string = 'x' + pow_r + str(i_ptr)
#            print x_val,pow_string
            if x_val == 0:
                pass
            elif x_val == 1:
                if p_string != '':
                    p_string += ' +'
                if pow_string == '':
                    p_string += '1'
                else:
                    p_string += pow_string
            elif x_val == -1:
                p_string += ' -'
                if pow_string == '':
                    p_string += '1'
                else:
                    p_string += pow_string
            elif x_val > 0 and p_string != '':
                if pow_string == '':
                    p_string += ' +' + str(self[i_ptr])
                else:
                    p_string += ' +' + str(self[i_ptr]) + '*' + pow_string
            else:
                if pow_string == '':
                    p_string += ' ' + str(self[i_ptr])
                else:
                    p_string += ' ' + str(self[i_ptr]) + '*' + pow_string

            i_ptr += 1
        if len(p_string) > 40 and self.last is None:
            p_string += '+ ...'
        elif len(p_string) > 40 and i_ptr <= self.last:
            p_string += '+...+'+str(self[self.last])+'x'+pow_r+str(self.last)
        elif p_string == '':
            p_string = '0'
#        elif self.last is None or self.last > 1 or self.first <> self.last:
#            p_string = '('+p_string+')'
        return p_string

    def __len__(self):
        """Length of poly, None if undefined or infinite."""
        return len(self.coef)

    def sig(self, sig_len):
        """Signature of poly."""
        return frac(self, 1).sig(sig_len)

    def deg(self):
        """Degree of poly, None if unknown."""
        return self.last

    # def div_xm(self, m):
    #     """Divide polynomial by x^m - shift start by m."""
    #     self.first -= m
    #     self.poff -= m
    #     if self.last is not None:
    #         self.last -= m

    def mul_xm(self, m):
        """Mmultiply polynomial by x^m - shift start by m."""
        self.first += m
        self.poff += m
        if self.last is not None:
            self.last += m

    def drop(self, cnt=1):
        """Drop coefficients."""
        if len(self.coef) < cnt or cnt < 0:
            raise (IndexError)
        self.coef.drop(cnt)

    def dx(self):
        """Return Derivative."""
        return poly(p_differentiate(self), self.first-1)

    def ix(self, c0=0):
        """Return Integral."""
        return poly(p_integrate(self, c0), self.first)

    def conv(self, x_var=None):
        """Convergences for poly evaluated at x_var (x_var should be small)."""
        if x_var is None:
            idx = self.first
            poly_ap = poly([])
            while self.last is None or\
                    (self[idx] is not None and idx <= self.last):
                if self[idx] != 0:
                    poly_ap += poly([self[idx]], idx)
                    yield poly_ap
                idx += 1
        else:
            idx = self.first
            x_fact = 1
            if idx > 0:
                for _ in range(idx):
                    x_fact *= x_var
            elif idx < 0:
                for _ in range(-idx):
                    x_fact /= x_var
                #        x_fact = x_var ** idx
            p_sum = 0
            while self.last is None or\
                    (self[idx] is not None and idx <= self.last):
                if self[idx] != 0:
                    p_sum += self[idx] * x_fact
                    yield p_sum
                idx += 1
                x_fact *= x_var

    # def plot(self, xrange=(-2, 2), color='r'):
    #     """Plot a polynomial."""
    #     if isinstance(xrange, tuple):
    #         xrange = np.linspace(xrange[0], xrange[1], 100)
    #     yvals = []
    #     for xval in xrange:
    #         print('x=', xval)
    #         try:
    #             yvals.append(self(xval))
    #         except OverflowError:
    #             yvals.append(0)
    #             # pass
    #             print('overflow')
    #             # xrange = xrange[xrange != xval]
    #     plt.plot(xrange, yvals, color)
    #     print(xrange, yvals)
#     """Iterate for polynomial evaluation of frac."""
#     x_v = frac(f0)
#     x_v.norm()
#     x_pow = frac(1, 1)
#     p_sum = frac(0, 1)
#     p_small = False
#     if px.first < 0:
#         for _ in range(px.first, 0):
#             x_pow /= x_v
#     elif px.first > 0:
#         for _ in range(0, px.first):
#             x_pow *= x_v
#     pidx = px.first
#     while px.last is None or pidx <= px.last:
#         if px[pidx] is not None and px[pidx] != 0:
#             p_term = px[pidx] * x_pow
#             p_sum += p_term
#             p_sum.norm()
# #                    print(p_term, x_var)
#             if px.last is None:
#                 if p_small or abs(p_term.num * f0.den) < p_term.den:
#                     yield p_sum
#                     p_small = True
#         pidx += 1
#         x_pow *= x_v
#     #            print 'poly_out = ',p_sum
#     yield p_sum


# class poly_approx():
#     """Polynomial approximation."""
#     def __init__(self, ap_gen):


class pe(poly):
    """Poly element."""

    def __init__(self, c0=0, deg=0):
        """Create pe= c0*x^deg."""
        self.first = deg
        self.last = deg
        self.poff = deg
        self.coef = istream([c0])
        self.c0 = c0
        self.pow = deg


def pqt(vec0, vec1):
    """PQ transform iterator."""
    icnt = 1
    done = False
    while not done:
        if vec0.last is None or vec0.last >= vec0.first+icnt or \
           vec1.last is None or vec1.last >= vec1.first+icnt:
            term = vec1.c0 * vec0[vec0.first+icnt] - \
                   vec0.c0 * vec1[vec1.first+icnt]
            yield term
            icnt += 1
        else:
            done = True


class PQTx():
    """Rational Polynomial Transform to CF terms."""

    def __init__(self, num, den):
        """Instance of transform."""
        self.pvec = [den, num]

    def __getitem__(self, idx):
        """Index into vectors 0=den,1=num,..."""
        if idx < len(self.pvec):
            return self.pvec[idx]
        else:
            while idx > len(self.pvec) - 1:
                if self.pvec[-1] is None:
                    return None
                else:
                    self.pvec.append(poly(pqt(self.pvec[-2], self.pvec[-1])))
                    if self.pvec[-1] == 0:
                        self.pvec[-1] = None
            return self.pvec[idx]

    def __call__(self):
        """Return cf form."""
        deg0 = self.pvec[0].first
        deg1 = self.pvec[1].first
        nt = poly([self.pvec[1].c0], deg1-deg0)
        dt = self.pvec[0].c0
        icnt = 1
        done = False
        yield(1, 0)
        while not done:
            yield(nt, dt)
            icnt += 1
            if self[icnt] is None:
                done = True
            else:
                nt = poly([self.pvec[icnt].c0], self.pvec[icnt].first+1)
                dt = self[icnt-1].c0


def poly_pow(r, n):
    """Poly for r*x^n."""
    new_p = poly([r], n)
#    new_p.first = n
    new_p.last = n
    return new_p


class pq_iter():
    """Iterate pq."""

    def __init__(self, pq):
        """Create pq iterator."""
        self.p0 = pq.num
        self.q0 = pq.den
        self.pi = self.p0.first
        self.qi = self.q0.first
        self.qc = self.q0[self.qi]
        self.dx = self.pi - self.qi
        self.ic = self.pi
        self.cx = []

    def __next__(self):
        """Iterate."""
        pt = self.p0[self.ic]
#        print(pt)
        qic = self.qi + self.ic - self.pi
        for rx in self.cx:
            pt -= rx * self.q0[qic]
#            print("-", rx, "*", self.q0[qic])
            qic -= 1
        r = frac(pt, self.qc)
        self.cx.append(r)
        self.ic += 1
        return r


class polyq():
    """Quotient of two polynomials."""

    def __init__(self, p_num, p_den):
        """Instance of polyq."""
        self.frac = frac(p_num, p_den)
        self.num = self.frac.num
        self.den = self.frac.den
        self.dx = self.num.first - self.den.first
        self.px = poly(iter(self), self.dx)

    def __call__(self, xval):
        """Evaluate polyQ."""
#        print ('calling polyQ eval')
        return approx(pq_iter(self, xval))

    def __iter__(self):
        """Iterate through the coefficents of PolyQ."""
        pqi = pq_iter(self)
        while True:
            yield next(pqi)
    # def __iter__(self) :
    #     """iterator for K functions"""

        # p = self.num
        # q = self.den
        # my_rx = rxform(p,q)
        # done = False
        # if p.first == q.first :
        #     term =my_rx.reduce_p()
        #     yield term
        # else:
        #     yield (poly([1]),0)
        # while not done :
        #     term = my_rx.reduce_q()
        #     yield term
        #     term = my_rx.reduce_p()
        #     yield term
    def __repr__(self):
        """Representation."""
        return '( ' + repr(self.num) + ') / (' + repr(self.den) + ' )'

    # def plot(self, xrange=(-2, 2), color='r'):
    #     """Plot a polynomial."""
    #     if isinstance(xrange, tuple):
    #         xrange = np.linspace(xrange[0], xrange[1], 100)
    #     yvals = []
    #     for xval in xrange:
    #         try:
    #             yvals.append(self(xval))
    #         except(OverflowError):
    #             print('OverflowError', xval)
    #             yvals.append(yvals[-1])
    #     plt.plot(xrange, yvals, color)
    #     print(xrange, yvals)


class epoly(poly):
    """Exponential Generator polynomial."""

    def __init__(self, terms, cycle=None):
        self.terms = terms
        if cycle is None:
            cycle = len(terms)
        self.cycle = cycle
        poly.__init__(self, self.nexis)
#        print (self.first, self.last)

    def nexis(self):
        """Convert epoly terms to poly."""
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

    def dx(self):
        """Return Derivative wrt x."""
        if len(self.terms) == self.cycle:
            return epoly(self.terms[1:] + [self.terms[0]], self.cycle)
        else:
            return epoly(self.terms[1:], self.cycle)

    def ix(self, c0=0):
        """Integrate wrt x."""
        return epoly([c0]+self.terms, self.cycle)
