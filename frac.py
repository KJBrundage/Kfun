# -*- coding: utf-8 -*-
"""Fraction Class."""

import numbers
# from digits import i2s


def i2s(ival, base):
    """Convert int to str of digits in base."""
    n_str = ''
    while ival > 0:  # convert to 'base' digits
        qrem = ival % base
        ival = ival // base
        if qrem < 10:  # 0,1,...,9
            n_str = str(qrem) + n_str
        else:  # A,B,C,...
            n_str = chr(55 + qrem) + n_str
    if n_str == '':
        n_str = '0'
    return n_str

# import random
# from approx import approx


def gcd(a, b):
    """Greatest Common Denominator."""
    while b != 0:
        a, b = b, a % b
    return abs(a)


def gcd2(a, b):
    """Binary GCD."""
    if a < b:
        a, b = b, a
    while b != 0:
        if a % 2 == 0:
            if b % 2 == 0:  # both even
                return gcd2(a >> 1, b >> 1) * 2
            else:          # even,odd
                return gcd2(b, a >> 1)
        else:
            if b % 2 == 0:
                return gcd2(a, b >> 1)
            else:
                return gcd2(b, (a - b) >> 1)
    return a


def gcd3(a, b):
    """GCD based on shift."""
    pt = 1
    while b != 0:
        if a % 2 == 0:
            if b % 2 == 0:
                pt = pt << 1
                a = a >> 1
                b = b >> 1
            else:
                a = a >> 1
        else:
            if b % 2 == 0:
                b = b >> 1
            else:
                if a > b:
                    a, b = b, (a - b) >> 1
                else:
                    a, b = a, (b - a) >> 1
    return pt*a


# def tg(n, m):
#     for i in range(n):
#         g = gcd(random.randint(1, m), random.randint(1, m))

# def tg2(n,m):
#     for i in range(n):
#         g = gcd2(random.randint(1,m),random.randint(1,m))

# def tg3(n,m):
#     for i in range(n):
#         g = gcd3(random.randint(1,m),random.randint(1,m))


class frac():
    """Fraction class for poly coefficients."""

    def __init__(self, num, den=1, norm=True):
        """Fraction instance."""
        if isinstance(num, tuple):
            n0 = num[0]
            d0 = num[1]
        elif isinstance(num, frac):
            n0 = num.num
            d0 = num.den
        elif isinstance(num, numbers.Number):
            n0 = num
            d0 = 1
        else:
            n0 = num
            d0 = 1
#            print 'frac', type(num), num, type(den), den
#            raise (TypeError)
        if isinstance(den, frac):
            n0 *= frac(den).den
            d0 *= frac(den).num
        else:
            # elif isinstance(den, numbers.Number):  # ???
            d0 *= den

        if isinstance(d0, numbers.Number) and d0 < 0:
            n0 *= -1
            d0 *= -1
#                   Don't normalize by default, support accuracy
#                       use frac(n,d,True) or call norm() to normalize.
        if norm:
            if (isinstance(n0, numbers.Integral) and
               isinstance(d0, numbers.Integral)):
                g = gcd(n0, d0)
                if d0 < 0:
                    g *= -1
                if g != 0:
                    n0 //= g
                    d0 //= g
            # else:
            #     n0 = n0 / d0
            #     d0 = 1
#            print type(den), isinstance(den, poly), isinstance(den, kpoly)
#            raise (TypeError)
        self.num = n0
        self.den = d0

    def cmp(self, other):
        """Old __cmp__ returns -1 if <, 0 if == and +1 if >."""
#        print ('cmp ',type(self),' with ',type(other))
        if other is None:
            return 1
        if isinstance(other, frac):
            vs, vo = self.num * other.den, self.den * other.num
        else:
            vs, vo = self.num, self.den * other
        if vs == vo:
            return 0
        elif (not isinstance(vs, numbers.Integral) or
              not isinstance(vo, numbers.Integral) or vs < vo):
            return -1
        else:
            return 1

    def __lt__(self, other):
        """Less than."""
        return self.cmp(other) == -1

    def __le__(self, other):
        """Less than or equal to."""
        return self.cmp(other) <= 0

    def __eq__(self, other):
        """Equal to."""
        return self.cmp(other) == 0

    def __ne__(self, other):
        """Not equal."""
        return self.cmp(other) != 0

    def __ge__(self, other):
        """Greater than or equal to."""
        return self.cmp(other) >= 0

    def __gt__(self, other):
        """Greater than."""
        return self.cmp(other) == 1

    def __neg__(self):
        """Negate."""
        return frac(-self.num, self.den)

    def __add__(self, other):
        """Addition."""
        if isinstance(other, frac):
            num = self.num * other.den + self.den * other.num
            den = self.den * other.den
            return frac(num, den, True)
        elif isinstance(other, int):
            num = self.num + other * self.den
            den = self.den
            return frac(num, den, True)
        else:
            return (other * self.den + self.num) / self.den

    def __radd__(self, other):
        """Right add."""
        return self + other

    def __sub__(self, other):
        """Subtraction."""
        return self + other * (-1)

    def __rsub__(self, other):
        """Swap elements."""
        return self * (-1) + other

    def __mul__(self, other):
        """Multiplication."""
        if isinstance(other, frac):
            num = self.num * other.num
            den = self.den * other.den
            return frac(num, den, True)
        elif isinstance(other, numbers.Integral):
            num = self.num * other
            den = self.den
            return frac(num, den, True)
        elif isinstance(other, float):
            return frac(other)*self
        else:
            # print("in frac mul", type(self.num), type(self.den),
            #       type(other), other)
            q = other * self
            return q

    def __rmul__(self, other):
        """Right multiplication."""
        return self * frac(other)

    def __div__(self, other):
        """Divide."""
        if isinstance(other, int):
            return frac(self.num, self.den * other)
        elif isinstance(other, frac):
            return self * ~other
            # elif not isinstance(self.den, int) and isinstance(other, int):
            #            print 'div',self,other  # ???
        else:
            return self.num / (self.den * other)

    def __rdiv__(self, other):
        """Swaped elemants."""
        return ~self * other
    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __invert__(self):
        """Invert."""
        return frac(self.den, self.num)

    def __round__(self):
        """Round."""
        return (2*self.num+self.den) // (2*self.den)

    def __float__(self):
        """Floating point."""
        try:
            result = float(self.num / self.den)
        except TypeError:
            result = 0
            #            return self.num/self.den
        except FloatingPointError:
            result = 0.
        return result

    def __trunc__(self):
        """Truncate."""
        if self.num is int and self.den is int:
            return self.num // self.den
        else:
            return int(self.num / self.den)

    def __abs__(self):
        """Absolute value."""
        return frac(abs(self.num), abs(self.den))

    def __iter__(self):
        """Iterate."""
        a, b = self.num, self.den
        cterm = a // b
        if isinstance(cterm, float):
            cterm = int(cterm)
        # elif not isinstance(cterm, int): # ???
#        print a, b, cterm, 'in __iter__'
        a, b = b, a - b*cterm
        yield (1, cterm)
        while b != 0:
            cterm = a // b
            if isinstance(cterm, float):
                cterm = int(cterm)
            a, b = b, a - b*cterm
#            print a, b, cterm
            yield (1, cterm)

    def __pow__(self, other):
        """Power."""
        if isinstance(other, int) or isinstance(other, float):
            if other == 0:
                return 1
            elif other > 0:
                return self*(self**(other-1))
            else:
                return ~self**other
        else:
            print('not implimented', type(other))
            raise TypeError

    def __xor__(self, other):
        """^ form."""
        return self**other

    def __getitem__(self, idx=None):
        """Index into stream."""
        if idx is None:
            return self
        elif isinstance(idx, slice):
            sigl = idx.stop
            my_sig = self.sig(sigl+1)
            return my_sig[idx.start:idx.stop]
#        if idx < self.start :
#            if self.end == None :
#                raise(IndexError)
#            else :
#                return self.p_coef[self.end+1+idx]
        if idx < 0:
            return None
        sigl = idx + 2
        my_sig = self.sig(sigl)
#        print("sig",my_sig)
        return my_sig[idx]

    def __call__(self, other=None):
        """Evaluate returns an approximator."""
        if other is None:
            return self
        else:
            # print("calling with ", other, type(self.num), type(self.den))
            if isinstance(self.num, int):
                return self.num/self.den
            else:
                return self.num(other)/self.den(other)
#        print('eval frac')
#         if isinstance(self.num, numbers.Number) \
#                 and isinstance(self.den, numbers.Number):
#             return self.num/self.den
#         else:
#             done = False
#             n_val, d_val = 0, 0
#             print("frac eval", self.num, self.den)
#             n_app = self.num(other)
#             d_app = self.den(other)
#             my_app = n_app/d_app
#             while not done:
#                 try:
#                     print("eval", n_val, d_val, done)
#                     yield +my_app
#                 except StopIteration:
#                     print("eval done", n_val, d_val)
# #                    yield n_val/d_val
#                     done = True

    def __repr__(self):
        """Representation."""
        if self.den == 0:
            if self.num == 0:
                return "\u2426"
            elif isinstance(self.num, numbers.Number) and self.num < 0:
                return "-\u221E"
            else:
                return "\u221E"
        if self.den == 1:
            return str(self.num)
        else:
            if isinstance(self.num, numbers.Integral):
                if self.num < 0:
                    return "-(" + str(-self.num) + " / " + str(self.den) + ")"
                else:
                    return "(" + str(self.num) + " / " + str(self.den) + ")"
            elif hasattr(self.num, "first")\
                    and self.num.first == self.num.last:
                return "(" + str(self.num) + " / " + str(self.den) + ")"
            else:
                return "(" + str(self.num) + ") / (" + str(self.den) + ")"

    def norm(self):
        """Normalize for Pade expressions."""
        if hasattr(self.den, 'first'):
            mx = frac(1, self.den[self.den.first])
            self.num *= mx
            self.den *= mx
        else:
            if isinstance(self.num, numbers.Integral) and \
                    isinstance(self.den, numbers.Integral):
                mx = gcd(self.num, self.den)
                self.num = int(self.num/mx)
                self.den = int(self.den/mx)

    def real(self):
        """Return floationg point representation."""
        return float(self.num) / self.den

    def out(self, dout, base=10):
        """Return string representation of frac."""
        gen = digits(self, base)
        out_str = ""
        for _ in range(dout):
            try:
                out_str += next(gen)
            except StopIteration:
                break
        return out_str

    def sig(self, sigl=None, big_n=None):
        """Signature of a fraction."""
        sig = []
        tnum, tden = self.num, self.den
        first_term = tnum // tden
        if isinstance(first_term, numbers.Number):
            nt = first_term
        elif hasattr(first_term, "first") and first_term.first == 0:
            nt = first_term
        else:
            nt = 0
        while tden != 0 and (sigl is None or len(sig) <= sigl):
            if big_n is not None and nt > big_n:
                print("Infinity and beyond")
                tden = 0
                if len(sig) == 0:
                    sig.append(0)
            else:
                sig.append(nt)
                tnum, tden = tden, tnum - nt*tden
            if tden == 0:
                nt = None
            else:
                nt = tnum//tden
        while len(sig) > 1 and sig[-1] == 1:
            sig = sig[:-2]+[sig[-2]+1]
        return sig

    def close(self, inum):
        """Is the fraction within 1 of inum?."""
        if abs(self.num - self.den * inum) < self.den:
            return True
        else:
            return False


def digits(var, base=10):
    """Generate a stream of digits representing a fraction."""
    num = var.num
    den = var.den
    if num*den < 0:
        neg = True
    else:
        neg = False
    num, den = abs(num), abs(den)
    t0 = int(num/den)
    s0 = i2s(t0, base)
    num -= t0*den
    num *= base
    if neg:
        yield "-"+s0+"."
    else:
        yield s0+"."
    while num > 0:
        t0 = int(num/den)
        s0 = i2s(t0, base)
        num -= t0*den
        num *= base
        yield s0


def apgen(var, mx=[[1, 0], [0, 1]]):
    """Generate terms for approximation."""
    terms = iter(var)
    cx, c0 = mx[0]
    dx, d0 = mx[1]
    done = False
    while not done:
        try:
            a, b = next(terms)
            cx, c0 = a*c0 + b*cx, cx
            dx, d0 = a*d0 + b*dx, dx
            yield frac(cx, dx)
        except StopIteration:
            done = True


def same_sig(sig_1, sig_2):
    """Compare signatures of two fractions."""
    idx = 0
    while len(sig_1) > idx and len(sig_2) > idx and (sig_1[idx] == sig_2[idx]):
        idx += 1
    return idx
