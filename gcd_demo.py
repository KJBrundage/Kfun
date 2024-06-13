# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 18:55:26 2023

@author: kjbru_000
"""


def gcd(a, b):
    """Return GCD. Simple GCD."""
    while b != 0:
        a, b = b, a % b
    return a


def show_gcd(a, b=1):
    """Greatest Common Denominator(a, b)."""
    a0, b0 = a, b
    cf_terms = []
    if isinstance(a, int) and isinstance(b, int):
        out_fmt = "{0:10d} \t/ {1:10d} \t={2:5d}\t+ {3:10d}"
        print("      a     \t    b   \t  a//b\t     a%b")
    else:
        out_fmt = "{0:2.10f} \t/ {1:2.10f} \t={2:5d}\t+ {3:0.10f}"
        print("      a     \t\t    b   \t  a//b\t     a%b")
    while b > 1E-8:
        print(out_fmt.format(a, b, int(a // b),  a % b))
        cf_terms.append(int(a//b))
        a, b = b, a % b
    if isinstance(a, int):
        print("GCD =", a)
    if not isinstance(a, int) or a == 1:
        if b0 != 1:
            cf_str = str(a0) + " / " + str(b0) + " = "
        else:
            cf_str = str(a0) + " = "
    else:
        if b0 != 1:
            cf_str = str(a0) + " / " + str(b0) + " = " + str(int(a0/a)) +\
                " / " + str(int(b0/a)) + " = "
        else:
            cf_str = str(a0) + " = " + str(int(a0/a)) +\
                " / " + str(int(b0/a)) + " = "
    for term in cf_terms:
        cf_str += str(term) + " +1/("
    cf_str = cf_str[:-5]
    if isinstance(a, int):
        for _ in cf_terms:
            cf_str += ")"
        cf_str = cf_str[:-1]
    else:
        cf_str += "+..."
    print(cf_str)


def prime(n):
    """Return Generator function for primes < n."""
    yield 2
    next_try = 3
    pp = 1
    while n == 0 or next_try <= n:
        if gcd(pp, next_try) == 1:    # next_try is prime!
            if n == 0 or next_try * next_try < n:
                pp *= next_try
            yield next_try
        next_try += 2
