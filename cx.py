# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:34:51 2019.

@author: Kev
"""
import frac
import numbers
global cxi_r


def cint(val):
    """Integer of complex rational."""
    return val.cint()


def set_cxi_r(char):
    """Set the power string."""
    global cxi_r
    cxi_r = char


cxi_r = 'i'   # default value for complex charactor


class cx():
    """Complex type."""

    def conj(self):
        """Conjunct complex."""
        if isinstance(self, cxi):
            return cxi(self.re, -self.im)
        else:
            return cxr(self.re, -self.im, self.den)


class cxi(cx):
    """Complex integers."""

    def __init__(self, real, imag=0):
        """Complex integer constructor."""
        if isinstance(real, cxi):
            self.re = real.re
            self.im = real.im
        elif isinstance(real, cxr):  # round re and im
            self.im = (2*real.im + real.den) // (2*real.den)
            self.re = (2*real.re + real.den) // (2*real.den)
        else:
            self.re = real
            self.im = imag

#            return 1 #?cxr

    def __eq__(self, other):
        """Equal."""
        if isinstance(other, cxi):
            if self.re == other.re and self.im == other.im:
                return True
            else:
                return False
        else:
            return self == cxi(other)

    def __ne__(self, other):
        """Not equal."""
        return not self == other

    def __gt__(self, other):
        """Greater than."""
        if isinstance(other, cxi):
            if self.re**2 + self.im**2 > other.re**2 + other.im**2:
                return True
            else:
                return False
        else:
            return self > cxi(other)

    def __ge__(self, other):
        """Greater than or equal to."""
        if isinstance(other, cxi):
            if self.re**2 + self.im**2 >= other.re**2 + other.im**2:
                return True
            else:
                return False
        else:
            return self >= cxi(other)

    def __lt__(self, other):
        """Less than."""
        if isinstance(other, cxi):
            if self.re**2 + self.im**2 < other.re**2 + other.im**2:
                return True
            else:
                return False
        else:
            return self < cxi(other)

    def __neg__(self):
        """Cnegate."""
        return cxi(-self.re, -self.im)

    def __add__(self, other):
        """Addition."""
        if isinstance(other, cxi):
            re = self.re + other.re
            im = self.im + other.im
            return cxi(re, im)
        elif isinstance(other, int):
            re = self.re + other
            return cxi(re, self.im)
        elif isinstance(other, frac.frac):  # ???
            re = self.re + int(other + frac.frac(1, 2))
            return cxi(re, self.im)
        else:
            print('not found _add_ type=', type(other))
            return other + self

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
        # print("cxi mult", type(other))
        if isinstance(other, cxi):
            re = self.re * other.re - self.im * other.im
            im = self.re * other.im + self.im * other.re
            return cxi(re, im)
        elif isinstance(other, cxr):
            re = self.re * other.re - self.im * other.im
            im = self.re * other.im + self.im * other.re
            return cxr(re, im, other.den)
        elif isinstance(other, numbers.Integral):
            return cxi(self.re * other, self.im * other)
        elif isinstance(other, float):
            return cxi(self.re * other, self.im * other)
        else:
            return other*self
        #     raise TypeError

    def __rmul__(self, other):
        """Right multiply."""
        return self * other

    def __pow__(self, other):
        """Power."""
        if isinstance(other, int) or isinstance(other, float):
            if other == 0:
                return 1
            elif other > 0:
                return self*(self**(other-1))
            else:
                return (self**(other+1))/self

    def __xor__(self, other):
        """Exclusive or."""  # power???
        return self**other

    def __div__(self, other):
        """Divide - returns cxr."""
        if isinstance(other, cxi):
            re = self.re * other.re + self.im * other.im
            im = self.im * other.re - self.re * other.im
            mx = other.re * other.re + other.im * other.im
            return cxr(re, im, mx)
        elif isinstance(other, numbers.Integral):
            return cxr(self) / cxr(other)
        elif isinstance(other, cxr):
            return cxr(self) / other
        # print ('cxi div by',type(other),'undefined')

    def __rtruediv__(self, other):
        """Swaped elemants."""
        return cxi(other) / self

    def __truediv__(self, other):
        """Truediv vectored to div."""
        return self.__div__(other)

    def __floordiv__(self, other):
        """Floor returns nearest cxi."""
        if isinstance(other, cxi):
            re = self.re * other.re + self.im * other.im
            im = self.im * other.re - self.re * other.im
            mx = other.re * other.re + other.im * other.im
            re = (2*re+mx) // (2*mx)
            im = (2*im+mx) // (2*mx)
            return cxi(re, im)
        elif isinstance(other, cxr):
            return cxr(self)//other
        elif isinstance(other, numbers.Integral):
            return self // cxi(other)
        #        elif isinstance(other,frac.frac) :
        #            return self*frac.frac(other.den,other.num)
        # else:
        #     raise TypeError

    def __rfloordiv__(self, other):
        """Floor div with cxi on right."""
        # print("cxi r div", type(other))
        return cxi(other) // self
#    def __round__(self, other):
#        """round to cxi"""
#        if isinstance(other, cxi):
#            re = self.re * other.re + self.im * other.im
#            im = self.im * other.re - self.re * other.im
#            mx = other.re * other.re + other.im * other.im
#            re = (2 * re + mx) // (2 * mx)
#            im = (2 * im + mx) // (2 * mx)
#            return cxi(re, im)
#        elif isinstance(other, numbers.Integral):
#            return self / cxi(other)
#        elif isinstance(other, cxr) :
#            return cxr(self) / other
#        print ('cxi div by',type(other),'undefined')
#        #        elif isinstance(other,frac.frac) :
#        #            return self*frac.frac(other.den,other.num)
#        #        else :
#        #            return
#        #            raise TypeError
#

    def __mod__(self, other):
        """Modulo."""
        if isinstance(other, cx):
            return self - round(self / other) * other
        elif isinstance(other, numbers.Integral):
            return self % cxi(other, 0)
        #    def __rmod__(self,other) :
        #        return cmplx(other)%self

    def __invert__(self):
        """Invert."""
        return cxi(self.re, -self.im)

    def __float__(self):
        """Floating point."""
        return float(self.re)  # + self.im * (0 + 1J)

    def __trunc__(self):
        """Truncate."""
        return cxi(int(self.re+frac.frac(1, 2)), int(self.im+frac.frac(1, 2)))

    def __abs__(self):
        """Absolute value."""
        return cxi(abs(self.re), abs(self.im))

    def __len__(self):
        """Lenght (magnitude) of a complex integer."""  # ???
        return int(self.re**2 + self.im**2)

    def __round__(self):
        """Round cxi returns itself."""
        return self

    def __repr__(self):
        """Representation."""
        global cxi_r
        if cxi_r == ',':
            return 'cxi(' + str(self.re) + ',' + str(self.im) + ')'
        else:
            if self.im == 0:
                return str(self.re)
            elif self.re == 0:
                if self.re == 1:
                    return cxi_r
                elif self.re == -1:
                    return "-"+cxi_r
                else:
                    return str(self.im) + cxi_r
            elif self.im < 0:
                return '(' + str(self.re) + str(self.im) + cxi_r + ')'
            else:
                return '(' + str(self.re) + '+' + str(self.im) + cxi_r + ')'

#    def d2(self):
#        """return len^2"""
#        return self.re * self.re + self.im * self.im

    def real(self):
        """Real part of cxi."""
        return self.re  # +float(self.im)*(0+1J)

    def imag(self):
        """Imaginary part of cxi."""
        return self.im

    def cx(self):
        """Machine complex."""
        return self.re + (0+1J)*self.im

    def d2(self):
        """Len squared."""
        return (self.re*self.re + self.im*self.im)


class cxr(cx):
    """Complex rational."""

    def __init__(self, real, imag=0, denom=1):
        """Complex rational constructor."""
        if isinstance(real, cxi) and imag == 0:
            self.re = real.re
            self.im = real.im
            self.den = denom
        else:
            if denom == 0:
                # return +/- infinity
                if real < 0:
                    real = -1
                elif real > 0:
                    real = 1
                if imag < 0:
                    imag = -1
                elif imag > 0:
                    imag = 1
            elif denom < 0:
                real *= -1
                imag *= -1
                denom *= -1
            # gcd0 = frac.gcd(real, imag)
            # gcd1 = frac.gcd(denom, gcd0)
            self.re = real
            self.im = imag
            self.den = denom
            # if gcd1 != 0:
            #     self.re //= gcd1
            #     self.im //= gcd1
            #     self.den //= gcd1

    def __eq__(self, other):
        """Equal."""
        if isinstance(other, cxr):
            if self.re * other.den == other.re * self.den and \
               self.im * other.den == other.im * self.den:
                return True
            else:
                return False
        else:
            return self == cxr(other)

    def __ne__(self, other):
        """Not equal."""
        return not self == other

    def __neg__(self):
        """Cnegate."""
        return cxr(-self.re, -self.im, self.den)

    def __add__(self, other):
        """Addition."""
        if isinstance(other, cxi):
            re = self.re + other.re * self.den
            im = self.im + other.im * self.den
            return cxr(re, im, self.den)
        elif isinstance(other, cxr):
            re = self.re*other.den + other.re*self.den
            im = self.im*other.den + other.im*self.den
            den = self.den*other.den
            return cxr(re, im, den)
        elif isinstance(other, numbers.Integral):
            re = self.re + other * self.den
            return cxr(re, self.im, self.den)
        elif isinstance(other, frac.frac):
            re = self.re * other.den + other.num * self.den
            im = self.im * other.den
            den = self.den * other.den
            return cxr(re, im, den)
        elif isinstance(other, cxr):
            re = self.re * other.den + other.re * self.den
            im = self.im * other.den + other.im * self.den
            den = self.den * other.den
            return cxr(re, im, den)
        else:
            raise TypeError  # type(other)

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
        if isinstance(other, cxi):
            re = self.re * other.re - self.im * other.im
            im = self.re * other.im + self.im * other.re
            den = self.den
            return cxr(re, im, den)
        elif isinstance(other, cxr):
            re = self.re * other.re - self.im * other.im
            im = self.re * other.im + self.im * other.re
            den = self.den * other.den
            return cxr(re, im, den)
        elif isinstance(other, numbers.Integral):
            return cxr(self.re * other, self.im * other, self.den)
        elif isinstance(other, frac.frac):
            return cxr(self.re*other.num, self.im*other.num,
                       self.den*other.den)
        else:
            raise TypeError

    def __rmul__(self, other):
        """Right divide."""
        return self * other

    def __div__(self, other):
        """Divide."""
        # print('cxr div')
        if isinstance(other, cxi):
            if other == 0:
                return cxr(self.re, self.im, 0)
            else:
                re = self.re * other.re + self.im * other.im
                im = self.im * other.re - self.re * other.im
                mx = other.re * other.re + other.im * other.im
                den = self.den * mx
                return cxr(re, im, den)
        elif isinstance(other, cxr):
            if other == 0:
                return cxr(self.re, self.im, 0)
            else:
                re = self.re * other.re + self.im * other.im
                im = self.im * other.re - self.re * other.im
                mx = other.re * other.re + other.im * other.im
                den = mx * self.den
                re *= other.den
                im *= other.den
                return cxr(re, im, den)
        elif isinstance(other, numbers.Integral):
            return cxr(self.re, self.im, self.den*other)
        elif isinstance(other, frac.frac):
            return self * ~ other
        else:
            raise TypeError

    def cint(self):
        """Complex version of int, returns cxi."""
        re = int(self.re/self.den)
        im = int(self.im/self.den)
        return cxi(re, im)
#    def __rdiv__(self,other) :
#        """right div"""
#        print('call rdiv',self,other)
#        if isinstance(other,numbers.Integral):
#            return cxr(other,0,1) / self
#        else :
#            raise( TypeError)
#

    def __truediv__(self, other):
        """Vector to __div__."""
        return self.__div__(other)

    def __floordiv__(self, other):
        """Divide."""
        if isinstance(other, cxi):
            re = self.re * other.re + self.im * other.im
            im = self.im * other.re - self.re * other.im
            mx = other.re * other.re + other.im * other.im
            mx *= self.den
            if re < 0:
                re = re // mx + 1
            else:
                re = re // mx
            if im < 0:
                im = im // mx + 1
            else:
                im = im // mx
            # re = re // mx
            # im = im // mx
            return cxi(re, im)
        elif isinstance(other, cxr):
            re = self.re * other.re + self.im * other.im
            im = self.im * other.re - self.re * other.im
            mx = other.re * other.re + other.im * other.im
            den = mx * self.den
            re *= other.den
            im *= other.den
            re = re // den
            im = im // den
            return cxi(re, im)

        elif isinstance(other, int):
            return self // cxi(other)
        #        elif isinstance(other,frac.frac) :
        #            return self*frac.frac(other.den,other.num)
        else:
            print('type of ', type(other), 'not recognized in //')
            raise TypeError

    def __rtruediv__(self, other):
        """Swaped elemants."""
        return cxr(other) / self

    def __round__(self):
        """Round."""
        return cxi((2*self.re + self.den) // (2*self.den),
                   (2*self.im + self.den) // (2*self.den))
#        q = self.den // 2
#        r = (self.re+q) // self.den
#        i = (self.im+q) // self.den
#        return cxi(r,i)

    def __pow__(self, other):
        """Power - complex ^ int ."""
        if isinstance(other, int):
            if other > 0:
                return self ** (other - 1) * self
            elif other < 0:
                return self ** (other + 1) / self
            else:  # other == 0
                return 1
            #        elif isinstance(other,frac.frac) :
            #            return self*frac.frac(other.den,other.num)
        else:
            raise TypeError

    def __mod__(self, other):
        """Modulo of cxr."""
        if isinstance(other, cxi):
            return self - round(self / other) * other
        elif isinstance(other, cxr):
            return self - cxi(self / other) * other
        elif isinstance(other, numbers.Integral):
            return self % cxi(other, 0)
        #    def __rmod__(self,other) :
        #        return cmplx(other)%self

    def __invert__(self):
        """Invert."""
        return cxr(self.den*self.re, -self.den*self.im, self.re**2+self.im**2)

    def __float__(self):
        """Floating point."""
        return float(self.re)/float(self.den)  # + self.im * (0 + 1J)

    def __trunc__(self):
        """Truncate. only accepts integer results. Bummer."""
        #   cxi(int(self.re + frac.frac(1, 2)), int(self.im + frac.frac(1, 2)))
        # return cxi(int(self.re / self.den), int(self.im / self.den))
        return int(self.re/self.den)

    def __iter__(self):
        """Iterate terms of continued fraction."""
        mval = self
        rval = round(mval)
        mval -= rval
        yield rval
        while mval != 0:
            mval = ~mval
            rval = round(mval)
            mval -= rval
            yield rval

    def __abs__(self):
        """Absolute value."""
        return cxi(abs(self.re), abs(self.im))

    def __len__(self):
        """Lenght (magnitude) of a complex integer."""
        return int(abs(self.re) + abs(self.im))

    def __repr__(self):
        """Representation."""
        global cxi_r
        if cxi_r == ',':
            return 'cxr(' + str(self.re) + ',' + str(self.im) + ',' + \
                str(self.den) + ')'
        else:
            if self.im == 0:
                return str(frac.frac(self.re, self.den))
            elif self.re == 0:
                return str(frac.frac(self.im, self.den)) + cxi_r
            elif self.im < 0:
                return '(' + str(frac.frac(self.re, self.den)) \
                       + str(frac.frac(self.im, self.den)) + cxi_r + ')'
            else:
                return '(' + str(frac.frac(self.re, self.den)) \
                       + '+'+str(frac.frac(self.im, self.den)) + cxi_r + ')'

    def d2(self):
        """Length squared."""
        return frac.frac(self.re*self.re+self.im*self.im, self.den*self.den)

    def real(self):
        """Real part of cxr."""
        return frac.frac(self.re, self.den)  # +float(self.im)*(0+1J)

    def imag(self):
        """Imaginary part of cxr."""
        return frac.frac(self.im, self.den)

    def cx(self):
        """Cast as machine complex."""
        return self.re/self.den + (0+1J)*self.im/self.den

    def close(self, other):
        """Is it close to other cxi."""
        if not isinstance(other, cxi):  # cost?
            raise TypeError
        rx = abs(self.re - self.den * other.re)
        ix = abs(self.im - self.den * other.im)
        if 3*rx < 2*self.den and 3*ix < 2*self.den:
            return True
        else:
            return False
