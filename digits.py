# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 21:25:25 2020.

Digit outputs

@author: Kev
"""
from ops import mxform
from cx import cxi, cxr


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


def digits(kp, base=10):
    """Return generator that iterates through digits of kcon."""
    k_stream = iter(kp.kiss)
    my_rx = mxform()
    r = None
    last = False  # last is returned from the yield to pass flag for last digit
    done = False  # set when last digit is returned
    # integer part
    while r is None:
        try:
            xmb = next(k_stream)
            if xmb == 0:
                my_rx.in_inf()
            else:
                my_rx.in_x(xmb)
        except(StopIteration):
            my_rx.in_inf()
        r = my_rx.drat(last)
    if isinstance(r, int):
        if r < 0:
            pos = False
            my_rx.mul_p(-1)
            r = -r - 1
            # print("negative", r)
        else:
            pos = True
        base_str = i2s(r, base)
        if pos:
            base_str = base_str + '.'
        else:
            base_str = '-' + base_str + '.'
        my_rx.sub_p(r)
        my_rx.mul_p(base)
        yield base_str
        while not done:
            r = None
            while r is None:
                try:
                    xmb = next(k_stream)
                    if xmb is None:
                        raise StopIteration
                    my_rx.in_x(xmb)
                #                    print 'in_x',xmb,my_rx
                except(StopIteration):
                    my_rx.in_inf()
                r = my_rx.drat(last)
#                print 'digits',my_rx.p, my_rx.q, r
            if last:  # if this is the last char (rounded) yield and terminate
                done = True
            else:
                my_rx.sub_p(r)
                my_rx.mul_p(base)
            if r == base:
                last = yield "*"  # round up last digit(s)
            elif r < 10:
                last = yield str(r)
            else:
                last = yield chr(55 + r)
    elif isinstance(r, cxi) or isinstance(r, cxr):
        r0, i0 = r.re, r.im
#        print("complex init",r0,i0)
        if r0 < 0:  # multiply to make real and imaginary parts positive
            if i0 < 0:
                my_rx.mul_p(-1)
                repos = False
                impos = False
                swap = False
            else:
                my_rx.mul_p(cxi(0, -1))
                repos = False
                impos = True
                swap = True
        else:
            if i0 < 0:
                my_rx.mul_p(cxi(0, 1))
                repos = True
                impos = False
                swap = True
            else:
                repos = True
                impos = True
                swap = False
        r = my_rx.drat(last)
        while r is None:
            try:
                xmb = next(k_stream)
                my_rx.in_x(xmb)
                # print( 'in_x',xmb,my_rx)
            except(StopIteration):
                my_rx.in_inf()
            r = my_rx.drat(last)
        re0, im0 = r.re, r.im
#        print("complex in",re0,im0)
        re_base = i2s(re0, base)
        im_base = i2s(im0, base)
        if swap:
            re_base, im_base = im_base, re_base
        if repos:
            re_base = re_base + '.'
        else:
            re_base = '-' + re_base + '.'
        if impos:
            im_base = im_base + '.'
        else:
            im_base = '-' + im_base + '.'
        my_rx.sub_p(cxi(re0, im0))
        my_rx.mul_p(base)
#        print("complex yield",re_base,im_base)
        yield (re_base, im_base)
        while not done:
            #            r = irat(my_rx.p,my_rx.q,last)
            old_r = None
#            new_r = my_rx.drat(last)
            new_r = None
            while (not isinstance(old_r, cxi)) or old_r != new_r:
                try:
                    xmb = next(k_stream)
                    my_rx.in_x(xmb)
#                    print( 'after in_x',xmb,my_rx) #len(my_rx.q[0])
                except(StopIteration):
                    my_rx.in_inf()
                r = my_rx.drat(last)
                old_r, new_r = new_r, r
#                print('new drat', my_rx.p[0]/my_rx.q[0],my_rx.p[1]/my_rx.q[1],
#                     last,old_r,new_r)
#                print my_rx.p[0]*my_rx.q[1]-my_rx.p[1]*my_rx.q[0]
            if last:  # if this is the last char (rounded) yield and terminate
                done = True
            else:
                my_rx.sub_p(r)
                my_rx.mul_p(base)
            re0, im0 = r.re, r.im
            if re0 == base:
                res = '*'
            elif re0 < 10:
                res = str(re0)
            else:
                res = chr(55 + re0)
            if im0 == base:
                ims = '*'
            elif im0 < 10:
                ims = str(im0)
            else:
                ims = chr(55 + im0)
            if swap:
                res, ims = ims, res
#            print("complex yield",res,ims)
            last = yield (res, ims)
