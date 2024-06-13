# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 11:04:23 2021.

@author: Kev
"""
from poly import poly, epoly
from frac import frac
from cx import cxi
from cf import K, kform


def k_p(var, min_rep=5, max_len=50, max_try=0):
    """Return the coeficients and cycle for continued fraction pattern."""
    min_cnt = 10
    if max_try == 0:
        if max_len is None:
            max_try = None
        else:
            max_try = (min_rep+1)*(max_len+1)
    kit = iter(var)
    c_0 = []            # coeficients of var
    g_0 = []            # 'guess vectors' and counts
    while max_try is None or len(c_0) < max_try:
        try:
            c_0.append(next(kit))
        except StopIteration:
            return None
        n_0 = []         # updated guess vector
        for off, cnt in g_0:
            clen = len(off)
            off0 = off[0]
#            print(off,cnt,c_0[-1],c_0[-clen-1])
            if c_0[-1][0] == c_0[-clen-1][0] and \
                    c_0[-1][1] == c_0[-clen-1][1] + off0:
                off = off[1:]
                off.append(off0)
                cnt += 1
            else:
                off = off[1:]
                off.append(c_0[-1][1] - c_0[-clen-1][1])  # corrected off
                cnt = 0
            if cnt >= min_rep*clen and cnt > min_cnt:   # found a good match
                while len(c_0) > (clen + 1) and \
                        c_0[-1][0] == c_0[-1-clen][0] and \
                        c_0[-1][1] == c_0[-1-clen][1]+off[-1]:    # wind back
                    c_0 = c_0[0:-1]
                    off.insert(0, off[-1])
                    off = off[0:-1]
                return (c_0, off)

            n_0.append((off, cnt))
        g_0 = n_0
        if max_len is None or len(c_0) <= max_len:         # add blank guess
            g_0.append(([0]*len(c_0), 0))
    return None


def k_pat(var, min_rep=10, max_len=50, max_try=0):
    """Attempt to identify a pattern in a K variable (simple type)."""
    kim = k_p(var, min_rep, max_len, max_try)
    if isinstance(kim, tuple):
        return K(kim[0], kim[1])
    return var


def s_p(var, min_rep=3, min_rl=10, max_len=50, max_try=150):
    """Attempt to identify a surd pattern (no incriments)."""
    kit = iter(var)
    k_t = []
    pguess = [0, 0]
    n_try = 0
    while max_try == 0 or n_try < max_try:
        k_t.append(next(kit))
        for i in range(2, len(pguess)):
            if k_t[-i] == k_t[-1]:
                pguess[i-1] += 1
                if pguess[i-1] >= min_rep*i and pguess[i-1] >= min_rl:
                    # print("in s_p ", k_t)
                    while len(k_t) > 2 and (k_t[-i] == k_t[-1]):
                        k_t = k_t[:-1]
                    return (k_t, i-1)
            else:
                pguess[i-1] = 0
        if len(pguess) < max_len:
            pguess.append(0)
        n_try += 1
    return None


def s_pf(var, min_rep=3, min_rl=10, max_len=50, max_try=None):
    """Attempt to identify a surd pattern quickly (no incriments)."""
    kit = iter(var)
    k_t = []
    for _ in range(min_rl):
        try:
            k_t.append(next(kit))
        except StopIteration:
            return None
    last_k_t = min_rep*max_len+min_rl
    pguess = [(2, 0)]
    icnt = 3
    while icnt <= last_k_t:
        # print("s_pf icnt =",icnt,pguess)
        try:
            k_t.append(next(kit))
        except StopIteration:
            return None
        new_guess = []
        for goff, gcnt in pguess:
            if k_t[-1] == k_t[-goff]:
                gcnt += 1
                new_guess.append((goff, gcnt))
                if gcnt >= (goff-1)*min_rep and icnt >= min_rl:  # good match
                    while len(k_t) > 2 and (k_t[-1] == k_t[-goff]):
                        k_t = k_t[:-1]
                    return (k_t, goff-1)
        if k_t[-1] == k_t[-icnt]:
            new_guess.append((icnt, 0))
        icnt += 1
        pguess = new_guess
    return None


def s_pat(var, min_rep=3, min_rl=10, max_len=50, max_try=150):
    """Create wrapper for s_p that returns K."""
    my_s_p = s_pf(var, min_rep, min_rl, max_len, max_try)
    if my_s_p is None:
        print("not a surd, try square")
        var2 = var*var
        my_s_p = s_pf(var2, min_rep, min_rl, max_len, max_try)
        if my_s_p is None:
            return None
        print("sqrt( K", my_s_p, ")")
        return sqrt(K(my_s_p[0], my_s_p[1]))
    print("K", my_s_p)
    return K(my_s_p[0], my_s_p[1])
#     if max_try == 0:
#         max_try = min_rep*(max_len+1)
#     kit = iter(var)
#     c_0 = []            # coeficients of var
#     g_0 = []            # 'guess vectors' and counts
#     while len(c_0) < max_try:
#         c_0.append(next(kit))
#         n_0 = []         # updated guess vector
#         for off, cnt in g_0:
#             clen = len(off)
#             off0 = off[0]
# #            print(off,cnt,c_0[-1],c_0[-clen-1])
#             if c_0[-1][0] == c_0[-clen-1][0] and \
#                     c_0[-1][1] == c_0[-clen-1][1] + off0:
#                 off = off[1:]
#                 off.append(off0)
#                 cnt += 1
#             else:
#                 off = off[1:]
#                 off.append(c_0[-1][1] - c_0[-clen-1][1])  # app corrected off
#                 cnt = 0

#             if cnt == min_rep*clen:    # found a good match - wind it back
#                 while len(c_0) > clen and \
#                         c_0[-1][0] == c_0[-1-clen][0] and \
#                         c_0[-1][1] == c_0[-1-clen][1]+off[-1]:
#                     c_0 = c_0[0:-1]
#                     off.insert(0, off[-1])
#                     off = off[0:-1]
#                 return K(c_0, off)

#             n_0.append((off, cnt))
#         g_0 = n0
#         if len(c_0) <= max_len:             # add new blank guess
#             g_0.append(([0]*len(c_0), 0))
#     return var


def pi_gen():
    """Iterate terms of continued fractions for pi."""
    yield ((1, 0))
    yield ((4, 1))
    n_i = 1
    d_i = 3
    while True:
        yield ((n_i, d_i))
        n_i += d_i
        d_i += 2


def pi_gens():
    """Iterate terms of continued fractions for pi."""
    yield ((1, 0))
    yield ((4, 1))
    idx = 1
    while True:
        yield ((idx**2, 2*idx+1))
        idx += 1


pigs = K(pi_gens)
pig = K(pi_gen)
pi = +pig
e = K([2, 1, 2, 1], [0, 2, 0])

exp = k_pat(kform(K(epoly([1])), True))

p_sin = epoly([0, 1, 0, -1])
p_cos = epoly([1, 0, -1, 0])
sin = K(p_sin)
cos = K(p_cos)

tan = k_pat(kform(K(p_sin/p_cos), True))
# cot = k_pat(kform(1/tan, True))

p_sinh = epoly([0, 1])
p_cosh = epoly([1, 0])
sinh = K(p_sinh)
cosh = K(p_cosh)

tanh = k_pat(K(p_sinh/p_cosh))
# coth = k_pat(1/tanh)


def sqrt_it(var):
    """Create iterator for square root of a variable."""
    look_ahead = 20
    # print("in sqrt var=", var.out(10))
    if var < 0:
        imag = True
        var *= -1
    else:
        imag = False
    if isinstance(var, int):
        sn_0 = 0
        while (sn_0+1)**2 <= var:
            sn_0 += 1
        if var == sn_0**2:
            if imag:
                yield (1, cxi(0, sn_0))
                # return K(cxi(0, sn_0))
            else:
                yield (1, sn_0)
                # return K([sn_0])
        else:
            if imag:
                yield (1, cxi(0, sn_0))
                rterm = (sn_0**2-var, cxi(0, 2*sn_0))
                while True:
                    yield rterm
            else:
                yield (1, sn_0)
                rterm = (var - sn_0**2, 2*sn_0)
                while True:
                    yield rterm
    elif isinstance(var, frac):
        t0 = 0
        while var.den*(t0+1)**2 <= var.num:
            t0 += 1
        if var == t0**2:
            if imag:
                yield (1, cxi(0, t0))
                # return K(cxi(0, sn_0))
            else:
                yield (1, t0)
                # return K([sn_0])
        else:
            if imag:
                yield (1, cxi(0, t0))
                rterm = var.num - var.den * t0 * t0
                print("in sqrt", var, t0, rterm)
                while True:
                    yield (cxi(0, rterm), cxi(0, 2*var.den*t0))
                    yield (cxi(0, rterm), cxi(0, 2*t0))
            else:
                yield (1, t0)
                rterm = var.num - var.den*t0
                print("in sqrt", var, t0, rterm)
                while True:
                    yield (rterm, 2*var.den*t0)
                    yield (rterm, 2*t0)
    else:
        sn_0 = 0
        while (sn_0+1)**2 <= var:
            sn_0 += 1
        sn_0_sig = [sn_0]+look_ahead*[1]
        itc = 0
        done = False
        while not done:
            sn_0 = K(sn_0_sig)
            sn1 = (sn_0 + var/sn_0)/2
            sn1_sig = sn1.sig(itc + look_ahead)
            # print("in sqrt itc=", itc, sn_0_sig, sn1_sig)
            while len(sn_0_sig) > itc and len(sn1_sig) > itc and \
                    sn_0_sig[itc] == sn1_sig[itc]:
                if isinstance(sn_0_sig[itc], tuple):
                    yield sn_0_sig[itc]
                else:
                    yield (1, sn_0_sig[itc])
                itc += 1
            sn_0_sig = sn1_sig


class sqrt(K):
    """Square root of continued fraction."""

    def __init__(self, var):
        """Create instance of square root."""
        self.var = var
        self.kid = K(sqrt_it(self.var))
        self.kiss = self.kid.kiss

    def __repr__(self):
        """Representation string for sqrt."""
        rstr = "sqrt("+repr(self.var)+")"
        return rstr
# def sqrt(n):
#     """Return sqare root of n."""
#     return K(sqrt_it(n))
    # if n < 0:
    #     imag = True
    #     n *= -1
    # else:
    #     imag = False
    # if isinstance(n, int):
    #     sn_0 = 0
    #     while (sn_0+1)**2 <= n:
    #         sn_0 += 1
    #     if n == sn_0**2:
    #         if imag:
    #             return K(cxi(0, sn_0))
    #         else:
    #             return K([sn_0])
    #     else:
    #         if imag:
    #             return K([cxi(0, sn_0), (sn_0**2-n, cxi(0, 2*sn_0))], 1)
    #         else:
    #             return K([sn_0, (n-sn_0**2, 2*sn_0)], 1)
    # else:
    #     sn_0 = 0
    #     while (sn_0+1)**2 <= n:
    #         sn_0 += 1
    #     if n == sn_0**2:
    #         if imag:
    #             return K(cxi(0, sn_0))
    #         else:
    #             return K([sn_0])
    #     else:
    #         sn_0 = K([sn_0])
    #         sn1 = (sn_0 + n/sn_0)/2
    #         la = 3
    #         icnt = 0
    #         done = False
    #         while not done:
    #             sig_0 = sn_0.sig(icnt + la)
    #             sig1 = sn1.sig(icnt + la)
    #             if len(sig_0) > icnt and len(sig1) > icnt \
    #                     and sig_0(icnt) == sig1(icnt):
    #                 yield


def atan_gen():
    """Terms of arctan."""
    i = 1
    while True:
        yield 0
        yield frac(1, i)
        i += 2
        yield 0
        yield frac(-1, i)
        i += 2


atan_p = poly(atan_gen)
atan = K(atan_p)


def asin_gen():
    """Terms of arcsin."""
    p_x = 0
    r_x = frac(1, 1)
    while True:
        yield 0
        p_x += 1
        yield r_x/p_x
        r_x *= p_x
        p_x += 1
        r_x /= p_x


asin_p = poly(asin_gen)
asin = K(asin_p)
acos = pi/2 - asin
