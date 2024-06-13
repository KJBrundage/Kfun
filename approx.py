# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 13:14:48 2022.

@author: kjbru_000
"""
from inspect import isgenerator
from frac import frac, same_sig
from series0 import istream
global fstep, vstep


def set_steps(new_fstep=None, new_vstep=None):
    """Set number of iterations with each update for iteration calculations."""
    global fstep, vstep
    if new_fstep is not None:
        fstep = new_fstep
    if new_vstep is not None:
        vstep = new_vstep


fstep, vstep = 1, 1


class approx():
    """Convergence value."""

    def __init__(self, apv, xvar=None):
        """Create instance of approximation."""
        self.val = None
        if isgenerator(apv):
            self.gen = istream(ita(apv))
            self.done = False
        elif isinstance(apv, int):
            self.val = apv
            self.done = True
        elif isinstance(apv, frac):
            self.val = apv
            self.done = True
        else:
            self.gen = istream(apgen(apv))
            self.done = False
        if xvar is not None:
            self.gen = istream(ita(itfx(apv, xvar)))
        self.var = apv
        self.xvar = xvar
        self.val = +self
        self.psig = []
        self.sig = self.val.sig()

    def __pos__(self):
        """Iterate approximation."""
        self.pval = self.val
        self.psig = self.sig
        if self.done:
            return self.val
            # if self.mom is not None:
            #     self.gen = +self.mom
            #     self.done = False
            #     self.val = +self
        else:
            try:
                self.val = next(self.gen)
                if isinstance(self.val, frac):
                    self.cval = self.val.den/self.val.num
                    self.sig = self.val.sig()
                else:
                    self.cval = frac(self.val, 1)
                    self.sig = self.cval.sig()
                # print("in approx, val=", self.val)
                # if isinstance(self.val, frac):
                #     self.sig = self.val.sig()
                # else:
                #     self.sig = self.val[0:10]  # ???fr
            except StopIteration:
                self.done = True
        return self.val

    def __invert__(self):
        """Return current approximation."""
        return self.val

    def __repr__(self):
        """Representation."""
        if isinstance(self.val, approx):
            return '~'+str(~self.val)
        else:
            return '~('+str(self.val)+')'

    def __add__(self, other):
        """Add approximation to other."""
        if isinstance(other, approx):
            new_v = approx(app_add(self, other))
            if self.done and other.done:
                new_v.done = True
        else:
            new_v = approx(app_add(self, other))
            if self.done:
                new_v.done = True
        return new_v

    def __radd__(self, other):
        """Right add approximation."""
        return self + other

    def __sub__(self, other):
        """Subtract other from approximation."""
        if isinstance(other, approx):
            new_v = approx(app_sub(self, other))
            if self.done and other.done:
                new_v.done = True
        else:
            new_v = approx(app_sub(self, other))
            if self.done:
                new_v.done = True
        return new_v

    def __rsub__(self, other):
        """Sight sub of approx."""
        return (-1 * self) + other

    def __mul__(self, other):
        """Multiply approximation and other."""
        if isinstance(other, approx):
            new_v = approx(app_mul(self, other))
            if self.done and other.done:
                new_v.done = True
        else:
            new_v = approx(app_mul(self, other))
            if self.done:
                new_v.done = True
        return new_v

    def __rmul__(self, other):
        """Right mult of approx."""
        return self * other

    def __div__(self, other):
        """Divide approximation by other."""
        if isinstance(other, approx):
            new_v = approx(app_div(self, other))
            if self.done and other.done:
                new_v.done = True
        else:
            new_v = approx(app_div(self, other))
            if self.done:
                new_v.done = True
        return new_v

    def __rdiv__(self, other):
        """Rdiv."""
        print('wtf')   # ??? wrong

    def __pow__(self, other):
        """Approximation value to integer power."""
        if isinstance(other, int):
            if other == 0:
                return approx(1)
            elif other == 1:
                return self
            elif other > 1:
                return self*(self**(other-1))
            else:
                return (1/self)**(-other)

    def __float__(self):
        """Floating point of current approximation."""
        return float(self.val)

    def __call__(self, val):
        """Return approximation with about that many digits accuracy."""
        f0 = ~self
        if isinstance(f0.den, int):
            digits = val
            while f0.den < 10**digits:
                f0 = +self
            return self.val
        else:
            return f0.num(val)/f0.den(val)
    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __getitem__(self, idx):
        """Return the cf terms from an approximation."""
        if isinstance(idx, slice):
            result = []
            my_start = 0
            if idx.start is not None:
                my_start = idx.start
            for i in range(my_start, idx.stop):
                try:
                    result.append(self[i])
                except IndexError:  # ??? is exception handled by getitem below
                    pass
            return result
        else:

            while (len(self.sig) <= idx+2 or len(self.psig) <= idx+2 or
                   self.sig[idx] != self.psig[idx] or
                   self.sig[idx+1] != self.psig[idx+1] or
                   self.sig[idx+2] != self.psig[idx+2]) and not self.done:
                +self
            if len(self.sig) >= idx:
                return (1, self.sig[idx])

    def __iter__(self):
        """Create an iterator for an approximations. Looks like a cf."""
        my_idx = 0
        done = False
        while not done:
            yield self[my_idx]
            my_idx += 1

    def out(self, ndig):
        """Output ndig of current value."""
        return self.val.out(ndig)

    def norm(self):
        """Normalize Pade expressions (set denominater = 1 + ...)."""
        return self.val.norm()

    def err(self):
        """Return the maximum error value for well behaved CFs and Pade."""
        if self.pval is None:
            return None
        else:
            return abs(self.val - self.pval)

    def sig(self, sigl=None):
        """Signature of an approximation."""
        return self.val.sig(sigl)

    def good_sig(self):
        """Return the trusted part of the approximations signature."""
        s0, s1 = self.pval.sig(), self.val.sig()
        ss = same_sig(s0, s1)
        return s0[0:ss]


def itfun(fun, var):
    """Approximator for function evaluated at var."""
    fs, vs, cb, fb = fstep, vstep, 3, 10
    big_n = 10000000000
    mysi = 0
    done = False
    # uv = True
    apfun = approx(fun)
    apvar = approx(var)
    for _ in range(fs):
        +apfun
    f0 = +apfun
    f1 = +apfun
    for _ in range(vs):
        +apvar
    v0 = +apvar
    v1 = +apvar
    # print("first est", f0, v0, f0(v0))
    s00 = f0(v0).sig(mysi+fb, big_n)
    s01 = f0(v1).sig(mysi+fb, big_n)
    s10 = f1(v0).sig(mysi+fb, big_n)
    s11 = f1(v1).sig(mysi+fb, big_n)
    # print(s00, '\n', s01, '\n', s10, '\n', s11)
    while not done:
        tsi = mysi + cb
        if tsi > len(s11):
            s11 = f1(v1).sig(mysi+fb, big_n)
        if tsi > len(s11):
            tsi = len(s11)
        ss00 = same_sig(s11, s00)
        ss10 = same_sig(s11, s10)
        ss01 = same_sig(s11, s01)
        ck00, ck10, ck01 = ss00 >= tsi, ss10 >= tsi, ss01 >= tsi
        # print(mysi,  ss01 >= tsi, ss10 >= tsi)
        if ck00 and ck01 and ck10:
            yield (1, s11[mysi])
            mysi += 1
            if mysi == len(s11):
                done = True
                # print("I'm outa here ", s11)
        else:
            # print('\n', s01, '\n', s10, '\n', s11)
            if not ck10:
                # print("update var", ss10 <= tsi)
                for _ in range(vs):
                    +apvar
                v0 = +apvar
                v1 = +apvar
                s10 = f1(v0).sig(mysi+fb, big_n)
                s01 = f0(v1).sig(mysi+fb, big_n)
                s11 = f1(v1).sig(mysi+fb, big_n)
            if not ck01:
                # print("update fun", ss01 <= tsi)
                for _ in range(fs):
                    +apfun
                f0 = +apfun
                f1 = +apfun
                s01 = f0(v1).sig(mysi+fb, big_n)
                s10 = f1(v0).sig(mysi+fb, big_n)
                s11 = f1(v1).sig(mysi+fb, big_n)
            # if not ck00 and ck01 and ck10:
            #     if uv:
            #         print("update var*", ss00 <= tsi, uv)
            #         v0, s00, s10 = v1, s01, s11
            #         v1 = +apvar
            #         s01 = f0(v1).sig(mysi+fb, big_n)
            #         s11 = f1(v1).sig(mysi+fb, big_n)
            #     else:
            #         print("update fun*", ss00 <= tsi, uv)
            #         f0, s00, s01 = f1, s10, s11
            #         f1 = +apfun
            #         s10 = f1(v0).sig(mysi+fb, big_n)
            #         s11 = f1(v1).sig(mysi+fb, big_n)
            #     uv = not uv


def itfx(fun, var):
    """Iterate over a function at a variable."""
    apfun = fun()
    apvar = var()
    for _ in range(vstep):
        v0 = +apvar
    for _ in range(fstep):
        f0 = +apfun
    v1 = +apvar
    f1 = +apfun
    idx = 0
    try_v = True
    check = False
    locked = False
    while not locked:
        print(f1, "(", v1, ")")
        a11 = f1(v1)
        a10 = f1(v0)
        while a11[idx] != a10[idx]:
            v0, v1 = v1, +apvar
            a10, a11 = a11, f1(v1)

        a01 = f0(v1)
        if a11[idx] != a01[idx]:
            f0, f1 = f1, +apfun
        else:
            a00 = f0(v0)
            if a11[idx] != a00[idx]:
                if try_v:
                    v0, v1 = v1, +apvar
                    try_v = False
                elif check is False:
                    f0, f1 = f1, +apfun
                    check = True
                else:
                    locked = True
            else:
                try_v = True
                check = False
                yield (1, a11[idx])
                idx += 1


def itf(fun, var):
    """Iterate over function at variable."""
    apfun = fun()
    apvar = var()
    v0 = +apvar
    v1 = +apvar
    f0 = +apfun
    f1 = +apfun
    a00, a01, a10, a11 = f0(v0), f0(v1), f1(v0), f1(v1)
    idx = 0
    try_v = True
    check = False
    locked = False
    while not locked:
        # print(idx, a00[idx], a01[idx], a10[idx], a11[idx])
        if a01[idx] != a11[idx]:
            f0, f1 = f1, +apfun
            a00, a01, a10, a11 = a10, a11, f1(v0), f1(v1)
            # print("new fun")
        elif a10[idx] != a11[idx]:
            v0, v1 = v1, +apvar
            a00, a01, a10, a11 = a01, f0(v1), a11, f1(v1)
            # print("new var")
        elif a00[idx] != a11[idx]:
            # print("cycle?")
            if try_v:
                v0, v1 = v1, +apvar
                try_v = False
            elif check is False:
                f0, f1 = f1, +apfun
                try_v = True
                check = True
            else:
                locked = True
                print("locked")
            a00, a01, a10, a11 = a11, f0(v1), f1(v0), f1(v1)
        else:
            print("sig =", idx, a11[0:idx+1])
            check = False
            yield a11[idx]
            idx += 1


def ita(ig_gen):
    """Return an approximator to a cf or poly generator."""
    first = True
    t0 = next(ig_gen)
    if isinstance(t0, tuple):
        h0, h1 = 0, 1
        k0, k1 = 1, 0
        a, b = t0
        h0, h1 = h1, a*h0 + b*h1
        k0, k1 = k1, a*k0 + b*k1
        done = False
        while not done:
            try:
                a, b = next(ig_gen)
                h0, h1 = h1, a*h0 + b*h1
                k0, k1 = k1, a*k0 + b*k1
                yield frac(h1, k1)
                first = False
            except StopIteration:
                if first:
                    yield frac(h1, k1)
                done = True
    elif isinstance(t0, approx):
        done = False
        while not done:
            # print("approx term")
            try:
                next_term = next(ig_gen)
                yield next_term
            except StopIteration:
                done = True
    else:
        sum = t0
        done = False
        while not done:
            try:
                next_term = next(ig_gen)
                sum += next_term
                yield sum
            except StopIteration:
                done = True


def apgen(var):
    """Create Generator function for approximation."""
    terms = iter(var)
    return ita(terms)
    # t0 = next(terms)
    # if isinstance(t0, tuple):
    #     a, b = t0
    #     cx, c0 = mx[0]
    #     dx, d0 = mx[1]
    #     cx, c0 = a*c0 + b*cx, cx
    #     dx, d0 = a*d0 + b*dx, dx
    #     done = False
    #     while not done:
    #         try:
    #             a, b = next(terms)
    #             cx, c0 = a*c0 + b*cx, cx
    #             dx, d0 = a*d0 + b*dx, dx
    #             yield frac(cx, dx)
    #         except StopIteration:
    #             done = True
    # else:
    #     sum = t0
    #     done = False
    #     while not done:
    #         try:
    #             next_term = next(terms)
    #             sum += next_term
    #             yield sum
    #         except StopIteration:
    #             done = True


def is_done(val0, val1):
    """Is the iteration done?."""
    if isinstance(val0, approx) and not val0.done:
        return False
    elif isinstance(val1, approx) and not val1.done:
        return False
    else:
        return True


def app_add(val0, val1):
    """Create Generator of approximator add."""
    done = False
    while not done:
        new_v = (+val0) + (+val1)
        if is_done(val0, val1):
            done = True
        yield approx(new_v)


def app_sub(val0, val1):
    """Create Generator of approximator sub."""
    done = False
    while not done:
        new_v = (+val0) - (+val1)
        if is_done(val0, val1):
            done = True
        yield approx(new_v)


def app_mul(val0, val1):
    """Create Generator for approxiator mult."""
    done = False
    while not done:
        new_v = (+val0) * (+val1)
        if is_done(val0, val1):
            done = True
        yield approx(new_v)


def app_div(val0, val1):
    """Create Generator for approxiator div."""
    done = False
    while not done:
        new_v = (+val0) / (+val1)
        if is_done(val0, val1):
            done = True
        yield approx(new_v)
