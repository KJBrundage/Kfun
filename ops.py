# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 20:54:42 2020.

@author: Kev
"""
# from numbers import Number
# from poly import poly
from cx import cxi, cint
from poly import poly
from frac import gcd


class mxform():
    """Mobius transform."""

    def __init__(self, xform=[[1, 0], [0, 1]]):
        """Create mobius xform [[cx, c0], [dx, d0]]."""
#        if isinstance(xform, mxform):
        self.mx = xform
        if not isinstance(xform, list):
            self.p = xform.p
            self.q = xform.q
        else:
            self.p = xform[0]
            self.q = xform[1]

    def __iter__(self):
        """Mobius transform iterator."""
        for coef in self.p:
            yield coef
        for coef in self.q:
            yield coef

    def mat(self):
        """Return the list form of this transform."""
        return [self.p, self.q]

    def norm(self):
        """Normalize the mobius transform."""
        g = None
        for coef in self:
            if coef != 0:
                if g is None:
                    g = coef
                else:
                    g = gcd(coef, g)
        if g > 1:
            self.p = [int(self.p[0]/g), int(self.p[1]/g)]
            self.q = [int(self.q[0]/g), int(self.q[1]/g)]

    def __mul__(self, other):
        """Compose twe mobius transforms."""
        if isinstance(other, mxform):
            newx = [[self.p[0]*other.p[0]+self.p[1]*other.q[0],
                     self.p[0]*other.p[1]+self.p[1]*other.q[1]],
                    [self.q[0]*other.p[0]+self.q[1]*other.q[0],
                     self.q[0]*other.p[1]+self.q[1]*other.q[1]]]
            return mxform(newx)
        else:
            return other*self

    def in_x(self, xmb):
        """Absorb term."""
        xm, b = xmb
        # print('xmb',xm,b)
        cx, c0 = self.p
        dx, d0 = self.q
        cx, c0 = cx * b + c0 * xm, cx
        dx, d0 = dx * b + d0 * xm, dx
        # print('in_x [[',cx,',',c0,'],[',dx,',',d0,']]')
        self.p = [cx, c0]
        self.q = [dx, d0]

    def in_inf(self):
        """Terminate input."""
        cx, c0 = self.p
        dx, d0 = self.q
        self.p = [cx, cx]
        self.q = [dx, dx]

    def mul_p(self, xm):
        """Multiply p terms."""
        for i in range(2):
            self.p[i] *= xm

    def mul_q(self, xm):
        """Multiply q terms."""
        for i in range(2):
            self.q[i] *= xm

    def sub_p(self, b):
        """Subtract b*q from p."""
        for i in range(2):
            self.p[i] -= self.q[i] * b

    def sub_q(self, b):
        """Subtract b*p from q."""
        for i in range(2):
            self.q[i] -= self.p[i] * b

    def swap(self):
        """Swap numerator and denominator."""
        self.p, self.q = self.q, self.p

    def rat(self, round=False):
        """Are the two ratios equal."""
        if self.q[0] == 0 or self.q[1] == 0:
            return None
        if isinstance(self.q[0], cxi):
            r0 = cxi(self.p[0]/self.q[0])
            r1 = cxi(self.p[1]/self.q[1])
        elif round:
            r0 = (self.q[0]+2*self.p[0]) // (2*self.q[0])
            r1 = (self.q[1]+2*self.p[1]) // (2*self.q[1])
        else:
            r0 = self.p[0] // self.q[0]
            r1 = self.p[1] // self.q[1]
        # print('rat',type(self.p[0]),r0,r1, self)
        # else :

        if r0 == r1:
            return r0
        else:
            return None

    def drat(self, round=False):
        """Ratio test for digits."""
        if self.q[0] == 0 or self.q[1] == 0:
            return None
        elif isinstance(self.q[0], cxi):
            if round:
                r0 = self.p[0] // self.q[0]
                r1 = self.p[1] // self.q[1]
            else:
                # r0 = int(self.p[0]/self.q[0])
                # r1 = int(self.p[1]/self.q[1])
                r0r = self.p[0]/self.q[0]
                r1r = self.p[1]/self.q[1]
                r0 = cint(r0r)
                r1 = cint(r1r)
        else:
            if round:
                r0 = (self.q[0]+2*self.p[0]) // (2*self.q[0])
                r1 = (self.q[1]+2*self.p[1]) // (2*self.q[1])
            else:
                r0 = self.p[0] // self.q[0]
                r1 = self.p[1] // self.q[1]
#        print('ratio test',round,self.p,self.q,r0,r1)
        if r0 == r1:
            return r0
        else:
            return None

    def __getitem__(self, items):
        """Provide bracket access to mobius transforms."""
        return self.mat().__getitem__(items)

    def __repr__(self):
        """Return representation of mxform."""
        return '[' + repr(self.p) + ', ' + repr(self.q)+']'


class txform():
    """Tensor transform."""  # xform = [[[cxy,cy],[cx,c0]],[[dxy,dy],[dx,d0]]]

    def __init__(self, xform):
        if not isinstance(xform, list):
            self.p = xform.p[:]
            self.q = xform.q[:]
        else:
            self.p = xform[0]
            self.q = xform[1]
        self.norm()
        self.bign = 10000000000  # ???

    def __iter__(self):
        """Iterate over all terms of tensor."""
        for elem in self.p:
            for coef in elem:
                yield coef
        for elem in self.q:
            for coef in elem:
                yield coef

    def __rmul__(self, other):
        """Compose mobius transform with tensor."""
        ax, a0, bx, b0 = other[0][0], other[0][1], other[1][0], other[1][1]
        cxy = ax*self[0][0][0] + a0*self[1][0][0]
        cy = ax*self[0][0][1] + a0*self[1][0][1]
        cx = ax*self[0][1][0] + a0*self[1][1][0]
        c0 = ax*self[0][1][1] + a0*self[1][1][1]
        # print(cxy, cy, cx, c0)
        dxy = bx*self[0][0][0] + b0*self[1][0][0]
        dy = bx*self[0][0][1] + b0*self[1][0][1]
        dx = bx*self[0][1][0] + b0*self[1][1][0]
        d0 = bx*self[0][1][1] + b0*self[1][1][1]
        # print(dxy, dy, dx, d0)
        return txform([[[cxy, cy], [cx, c0]], [[dxy, dy], [dx, d0]]])

    # def __rmul__(self, other):
    #     """Compose mobius transform with tensor. Mobius must be affine."""
    #     if other[1][0] == 0:
    #         ax, a0, b0 = other[0][0], other[0][1], other[1][1]
    #         cxy = ax*self[0][0][0] + a0*self[1][0][0]
    #         cy = ax*self[0][0][1] + a0*self[1][0][1]
    #         cx = ax*self[0][1][0] + a0*self[1][1][0]
    #         c0 = ax*self[0][1][1] + a0*self[1][1][1]
    #         print(cxy, cy, cx, c0)
    #         dxy = b0*self[1][0][0]
    #         dy = b0*self[1][0][1]
    #         dx = b0*self[1][1][0]
    #         d0 = b0*self[1][1][1]
    #         print(dxy, dy, dx, d0)
    #         return txform([[[cxy, cy], [cx, c0]], [[dxy, dy], [dx, d0]]])
    #     elif other[0][0] == 0:
    #         a0, bx, b0 = other[0][1], other[1][0], other[1][1]
    #         cxy = a0*self[1][0][0]
    #         cy = a0*self[1][0][1]
    #         cx = a0*self[1][1][0]
    #         c0 = a0*self[1][1][1]
    #         print(cxy, cy, cx, c0)
    #         dxy = bx*self[0][0][0] + b0*self[1][0][0]
    #         dy = bx*self[0][0][1] + b0*self[1][0][1]
    #         dx = bx*self[0][1][0] + b0*self[1][1][0]
    #         d0 = bx*self[0][1][1] + b0*self[1][1][1]
    #         print(dxy, dy, dx, d0)
    #         return txform([[[cxy, cy], [cx, c0]], [[dxy, dy], [dx, d0]]])

    def comp(self, x_var, y_var):
        """Compose mobius transform with tensor."""
        # print("txform mult", type(self), type(x_var), type(y_var))
        cxy, cy, cx, c0, dxy, dy, dx, d0 = iter(self)
        if x_var != 1:
            r, s, t, u = iter(x_var)
            cxy, cy, cx, c0 = \
                (cxy*r + cy*t), (cxy*s + cy*u), (cx*r + c0*t), (cx*s + c0*u)
            dxy, dy, dx, d0 = \
                (dxy*r + dy*t), (dxy*s + dy*u), (dx*r + d0*t), (dx*s + d0*u)
        if y_var != 1:
            r, s, t, u = iter(y_var)
            cxy, cy, cx, c0 = \
                (cxy*r + cx*t), (cy*r + c0*t), (cxy*s + cx*u), (cy*s + c0*u)
            dxy, dy, dx, d0 = \
                (dxy*r + dx*t), (dy*r + d0*t), (dxy*s + dx*u), (dy*s + d0*u)
        nt = [[[cxy, cy], [cx, c0]], [[dxy, dy], [dx, d0]]]
        return txform(nt)

    def mat(self):
        """Return the list version (2x2x2) of the tensor."""
        return [self.p, self.q]

    def norm(self):
        """Normalize the mobius transform."""
        g = None
        for coef in self:
            if coef != 0:
                if g is None:
                    g = coef
                else:
                    g = gcd(coef, g)
            if g == 1:
                break
        if g > 1:
            self.p = [[int(self.p[0][0]/g), int(self.p[0][1]/g)],
                      [int(self.p[1][0]/g), int(self.p[1][1]/g)]]
            self.q = [[int(self.q[0][0]/g), int(self.q[0][1]/g)],
                      [int(self.q[1][0]/g), int(self.q[1][1]/g)]]

    def in_x(self, mxb):
        """Input x value."""
        # print("in_x with ", mxb)
        mx, b = mxb
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        cx, c0 = c0*mx + cx*b, cx
        dx, d0 = d0*mx + dx*b, dx
        cxy, cy = cy*mx + cxy*b, cxy
        dxy, dy = dy*mx + dxy*b, dxy
        self.p = [[cxy, cy], [cx, c0]]
        self.q = [[dxy, dy], [dx, d0]]

    def in_xinf(self):
        """Terminate x stream."""
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        c0 = cx
        d0 = dx
        cy = cxy
        dy = dxy
        self.p = [[cxy, cy], [cx, c0]]
        self.q = [[dxy, dy], [dx, d0]]

    def in_y(self, mxb):
        """Input x value."""
        mx, b = mxb
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        cy, c0 = c0*mx + cy*b, cy
        dy, d0 = d0*mx + dy*b, dy
        cxy, cx = cx*mx + cxy*b, cxy
        dxy, dx = dx*mx + dxy*b, dxy
        self.p = [[cxy, cy], [cx, c0]]
        self.q = [[dxy, dy], [dx, d0]]

    def in_yinf(self):
        """Terminate y stream."""
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        c0 = cy
        d0 = dy
        cx = cxy
        dx = dxy
        self.p = [[cxy, cy], [cx, c0]]
        self.q = [[dxy, dy], [dx, d0]]

    def out(self, mx, b):
        """Output continued fraction term, b."""
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        c0, d0 = d0,  c0*mx - b*d0
        cx, dx = dx,  cx*mx - b*dx
        cy, dy = dy,  cy*mx - b*dy
        cxy, dxy = dxy, cxy*mx - b*dxy
        self.p = [[cxy, cy], [cx, c0]]
        self.q = [[dxy, dy], [dx, d0]]

    def rat(self, close=False):
        """Are the ratios the same."""
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        err = 0
        cval = 100
        rxy = None
        try:
            if isinstance(dxy, cxi):
                # cvar = True
                rxy = cxi(cxy/dxy)
                rx = cxi(cx/dx)
                ry = cxi(cy/dy)
                r0 = cxi(c0/d0)
#               print('ratfuck',rxy,rx,ry,r0)
            else:
                #  cvar = False
                rxy = cxy // dxy
                rx = cx // dx
                ry = cy // dy
                r0 = c0 // d0
        except ZeroDivisionError:
            if dxy == 0:
                err += 8
            if dx == 0:
                err += 2
            if dy == 0:
                err += 4
            if d0 == 0:
                err += 1
            if err != 0:
                return err
            # print('ratfuck', close, rxy, rx, ry, r0)
        if isinstance(rxy, poly):
            close = False
        else:
            if rxy > cval:
                if rx < -cval or ry < -cval or r0 < -cval:
                    return 0
        if rxy != rx:
            if close:
                if abs(cx - rxy*dx) >= abs(dx):
                    err += 2
            else:
                err += 2    # need in_y
        if rxy != ry:
            if close:
                if abs(cy - rxy*dy) >= abs(dy):
                    err += 4
            else:
                err += 4    # need in_x

        if rxy != r0:
            if close:
                if abs(c0 - rxy*d0) >= abs(d0):
                    err += 1
            else:
                err += 1
        if err == 0:
            if isinstance(rxy, int):
                return (1, rxy)
            else:
                return rxy
        else:
            return err

    def rat2(self, close=False):
        """Are the ratios the same optimized version."""
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        err = 0
        cval = 100
        rxy = None
        if dxy == 0:
            err += 8
        if dx == 0:
            err += 2
        if dy == 0:
            err += 4
        if d0 == 0:
            err += 1
        if err != 0:
            return [err, 0]
        else:
            if isinstance(dxy, cxi):
                # cvar = True
                rxy = cxi(cxy/dxy)
                rx = cxi(cx/dx)
                ry = cxi(cy/dy)
                r0 = cxi(c0/d0)
#               print('ratfuck',rxy,rx,ry,r0)
            else:
                #  cvar = False
                rxy = cxy // dxy
                rx = cx // dx
                ry = cy // dy
                r0 = c0 // d0
            # print('rat2fuck', close, rxy, rx, ry, r0)
            if isinstance(rxy, poly):
                close = False
            else:
                if rxy > cval:
                    if rx < -cval or ry < -cval or r0 < -cval:
                        return [0, 0]

            min, max = rxy, rxy

            if rxy != rx:
                if rx < min:
                    min = rx
                else:
                    max = rx
                if close:
                    if abs(cx - rxy*dx) >= abs(dx):
                        err += 2
                else:
                    err += 2    # need in_y
            if rxy != ry:
                if ry < min:
                    min = ry
                elif ry > max:
                    max = ry
                if close:
                    if abs(cy - rxy*dy) >= abs(dy):
                        err += 4
                else:
                    err += 4    # need in_x

            if rxy != r0:
                if r0 < min:
                    min = r0
                elif r0 > max:
                    max = r0
                if close:
                    if abs(c0 - rxy*d0) >= abs(d0):
                        err += 1
                else:
                    err += 1
        if err == 0:
            if isinstance(rxy, int):
                return (1, rxy)
            else:
                return rxy
        else:
            return [err, max-min]

    def inf(self):
        """Is the tensor diverging."""
        cxy, cy = self.p[0]
        cx, c0 = self.p[1]
        dxy, dy = self.q[0]
        dx, d0 = self.q[1]
        # print("test for inf in ops", cxy, dxy)
        if not isinstance(cxy, poly) and dxy == 0 or abs(cxy/dxy) > self.bign:
            # dxy*d0 < 0 and
            # print("infinity in ops")
            return True
        else:
            return False

    def __getitem__(self, items):
        """Provide bracket access to tensor items."""
        return self.mat().__getitem__(items)

    def __repr__(self):
        """Representation of txform."""
        return "[" + repr(self.p) + "," + repr(self.q) + "]"


class tensor():
    """Two variable tensor continued fraction operator."""

# tx = txform([[cxy,cy],[cx,c0]],[[dxy,dy],[dx,d0]])

    def __init__(self, tx, x_var, y_var):
        self.xv = x_var
        self.yv = y_var
        self.xs = iter(x_var.kiss)
        self.ys = iter(y_var.kiss)
        self.tx = tx
        self.xy = True

    def in_x(self):
        """Input x term."""
        try:
            xmb = next(self.xs)
            self.tx.in_x(xmb)
        except StopIteration:
            self.tx.in_xinf()

    def in_y(self):
        """Input y term."""
        try:
            xmb = next(self.ys)
            self.tx.in_y(xmb)
        except StopIteration:
            self.tx.in_yinf()

    def in_r(self):
        """Input x and y."""
        self.xy = not self.xy
        if self.xy:
            self.in_x()
        else:
            self.in_y()

    def in_2r(self):
        """Input x and y."""
        self.xy = not self.xy
        if self.xy:
            self.in_x()
            self.in_x()
        else:
            self.in_y()
            self.in_y()

    def rat(self, round=False):
        """Ratio test."""
        return self.tx.rat(round)

    def rat2(self, round=False):
        """Ratio test with range."""
        return self.tx.rat2(round)

    def inf(self):
        """Infinity."""
        return self.tx.inf()

    def out(self, mxb):
        """Output transform."""
        mx, b = mxb
        self.tx.out(mx, b)

    def __repr__(self):
        """Return representation for tensor."""
        """representation of tensor"""
        return "tensor(" + repr(self.tx) + ". " + \
            repr(self.xv) + ", " + repr(self.yv) + ")"
