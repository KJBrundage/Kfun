# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 19:57:39 2023.

@author: kjbru_000
"""

# from cf import K
from kmath import sp, sqrt, spat
from frac import gcd, frac
from cf import K


def surd_rst(m_val):
    """Determine rst values from quadratic surd."""
    """ x = (r + sqrt(s))/t : start by finding t"""
    t = 1
    kim = sp(m_val, 2, 20, 50, 100)
    if not isinstance(kim, tuple):
        print("not a surd")
        return None
    s_val = spat(m_val)
    kcoef, koff = kim
    # print(kim)
    # max_len = len(koff)
    done = False
    while(not done):
        t_val = spat(s_val)
        print("surd", type(kcoef), type(koff))
        while len(kcoef) > koff+1 or (kcoef[-1][1] % 2) != 0:
            if t > 0:
                t = -t
            else:
                t = 1-t
            t_val = t*s_val
            kim = sp(t_val)
            print(t, kim)
            if isinstance(kim, tuple):
                kcoef, koff = kim
            # else:
            #     print("not kim", t)
        print("found match t =", t, kcoef)
        print("c0=", kcoef[0][1], "-", kcoef[-1][1]/2)
        c0 = int(kcoef[-1][1] / 2)
        r = (kcoef[0][1] - c0)
        print("r =", r, t, t_val, type(t_val), type(r))
        s0 = t_val - r
        print("s0 =", s0)
        s2 = (s0*s0)
        its = iter(s2)
        t0 = next(its)
        s = t0[1]
        try:
            _ = next(its)
            done = False
            koff = 0
        except StopIteration:
            done = True
        msr = 1
        st = 2
    while abs(s) >= st*st:
        while s % st*st == 0:
            msr *= st
            s /= st*st
            s = int(s)
        if st == 2:
            st += 1
        else:
            st += 2
    if t < 0:
        r, msr, t = -r, -msr, -t
    return (r, msr, s, t)


class surd(K):
    """Quadratic surd class."""

    def __init__(self, v0, v1=None, v2=None, v3=None):
        """Initialize a qudratic surd."""
        """
          4 variables  x = surd(r,ms,s,t) = (r+ms*sqrt(s))/t
          3 variables  x = surd(a,b,c) = solution to ax^2 + bx + c = 0
          2 variables  x = surd(a,b) = (a+sqrt(b))
          1 variable
           x = surd(K) - attemp to find surd representing a continued fraction
        """
        if v3 is not None:
            self.sqrt = False
            self.r, self.ms, self.s, self.t = v0, v1, v2, v3
            self.fix_rst()
            self.set_abc()
            self.val = (self.r + self.ms * sqrt(self.s))/self.t
        elif v2 is not None:
            self.sqrt = False
            self.a, self.b, self.c = v0, v1, v2
            self.set_rst()
            self.fix_rst()
        elif v1 is not None:
            self.sqrt = False
            self.r, self.ms, self.s, self.t = v0, 1, v1, 1
            self.fix_rst()
            self.set_abc()
        else:
            try_simple = surd_rst(+v0)
            if try_simple is None:
                print("trying square")
                try_square = surd_rst(v0*v0)
                if try_square is None:
                    print("can't find solution")
                else:
                    print("found square solution", try_square)
                    self.r, self.ms, self.s, self.t = try_square
                    self.fix_rst()
                    self.set_abc()
                    self.sqrt = True
                    c = self.r*self.r - self.ms*self.ms*self.s
                    t1 = (self.r + sqrt(c))/2
                    t2 = (self.r - sqrt(c))/2
                    if self.ms < 0:
                        self.aval = sqrt(t1) - sqrt(t2)
                    else:
                        self.aval = sqrt(t1) + sqrt(t2)
                    print(self.aval)
            else:
                self.sqrt = False
                self.r, self.ms, self.s, self.t = try_simple
                self.fix_rst()
                self.set_abc()
        if self.sqrt:
            self.val = sqrt((self.r + self.ms * sqrt(self.s))/self.t)
            self.con = sqrt((self.r - self.ms * sqrt(self.s))/self.t)
        else:
            self.val = (self.r + self.ms * sqrt(self.s))/self.t
            self.con = (self.r - self.ms * sqrt(self.s))/self.t
        self.kid = self.val
        self.kiss = self.kid.kiss
        # self.ival = kp(self.val, 10)

    def set_abc(self):
        """Set values of a, b, and c based on r, ms, s, t values."""
        self.a, self.b, self.c = (
            self.t, -2*self.r, frac(self.r**2-self.s, self.t))

    def set_rst(self):
        """Set values of r, ms, s, and t based on a, b, c values."""
        self.r, self.ms, self.s, self.t = (
            -self.b, 1, (self.b*self.b - 4*self.a*self.c), 2*self.a)

    def fix_rst(self):
        """Remove square terms from s."""
        im = 2
        while abs(self.s) >= im*im:
            while self.s % im*im == 0:
                self.ms *= im
                self.s /= im*im
            if im == 2:
                im += 1
            else:
                im += 2
        g0 = gcd(abs(self.r), abs(self.ms))
        g1 = gcd(g0, self.t)
        if g1 > 1:
            self.r /= g1
            self.ms /= g1
            self.t /= g1
        if self.t < 0:
            self.r *= -1
            self.ms *= -1
            self.t *= -1

    def __repr__(self):
        """Representation of surd."""
        if self.t < 0:
            my_r, my_ms, my_s, my_t = -self.r, -self.ms, self.s, -self.t
        else:
            my_r, my_ms, my_s, my_t = self.r, self.ms, self.s, self.t
        if self.sqrt:
            r_str = "sqrt("
        else:
            r_str = "("
        if my_ms > 0:
            if my_ms != 1:
                if my_r == 0:
                    r_str += "("
                else:
                    r_str += "(" + repr(my_r) + " + "
                r_str += repr(my_ms) + " * sqrt(" + repr(my_s) + "))"
            else:
                if my_r == 0:
                    r_str += "("
                else:
                    r_str += "(" + repr(my_r) + " + "
                r_str += "sqrt(" + repr(my_s) + "))"
        else:
            if my_ms != -1:
                if my_r == 0:
                    r_str += "( - "
                else:
                    r_str += "(" + repr(my_r) + " - "
                r_str += repr(abs(my_ms)) + " * sqrt(" + repr(my_s) + "))"
            else:
                if my_r == 0:
                    r_str += " - sqrt(" + repr(my_s) + ")"
                else:
                    r_str += "(" + repr(my_r) + " - sqrt(" + repr(my_s) + "))"
        if my_t != 1:
            r_str += "/ " + repr(my_t)
        r_str += ")\n"

        r_str += "Solution to :"
        my_b, my_c = frac(self.b, self.a), frac(self.c, self.a)
        if self.sqrt:
            r_str += "x^4"
        else:
            r_str += "x^2"
        if my_b > 0:
            if my_b == 1:
                if self.sqrt:
                    r_str += " + x^2"
                else:
                    r_str += " + x"
            else:
                if self.sqrt:
                    r_str += " +" + repr(my_b) + "x^2"
                else:
                    r_str += " +" + repr(my_b) + "x"
        elif my_b < 0:
            if my_b == -1:
                if self.sqrt:
                    r_str += " - x^2"
                else:
                    r_str += " - x"
            else:
                if self.sqrt:
                    r_str += repr(my_b) + "x^2"
                else:
                    r_str += repr(my_b) + "x"
        if my_c > 0:
            r_str += " +" + repr(my_c)
        elif my_c < 0:
            r_str += " " + repr(my_c)
        r_str += " = 0 \n"
        r_str += "val = " + repr(spat(+self.val)) + '\n'
        r_str += "con = " + repr(spat(+self.con)) + '\n'
        return r_str
