# Factor RSA numbers using Gidney-Ekera algorithm
# https://arxiv.org/abs/1905.09749

from typing import Dict
from qualtran import Bloq, BloqBuilder, Signature, Register, SoquetT, QUInt
from qualtran.bloqs.data_loading.qrom import QROM

class RSA_GE(Bloq):
    """
        e,x -> e,x*g^e mod N
        N: an RSA integer to be factored
        g: random integer between 2 and N-1
        ew: exponent window size
        mw: multiplication window size
    """
    def __init__(self, N: int, g: int, ew: int, mw: int):
        self.N = N
        self.N_bitsize = N.bit_length()
        self.g = g
        self.ew = ew
        self.mw = mw
        
    @property
    def signature(self):
        return Signature([Register('ew', QUInt(self.ew)), Register('x', QUInt(self.N_bitsize))])
    
    def build_composite_bloq(self, bb: 'BloqBuilder', ew: 'Soquet', x: 'Soquet') -> Dict[str, 'SoquetT']:
        print("bitsize of", self.N, "is", self.N_bitsize)
        m = self.N_bitsize//2 + 1
        e1_bitsize = 2*m
        print("first component bitsize = ", e1_bitsize)
        ew, x = self.ekera_hastad_component(bb, self.g, e1_bitsize, ew, x)
        s = 1
        e2_bitsize = m//s
        print("second component bitsize = ", e2_bitsize)
        y = pow(self.g, self.N+1, self.N)
        ew, x = self.ekera_hastad_component(bb, y, e2_bitsize, ew, x)
        return {'ew': ew, 'x': x}

    def ekera_hastad_component(self, bb: 'BloqBuilder', g: int, e_bitsize: int, ew: 'Soquet', x: 'Soquet'):
        e_bits = list(range(e_bitsize))
        for i in range(0, e_bitsize, self.ew):
            gi = pow(g, 2**(e_bitsize-1-i), self.N)
            if len(e_bits[i:i+self.ew]) == self.ew:
                ew, x = bb.add(TimesExpMod(g=gi, N=self.N, e_bitsize=len(e_bits[i:i+self.ew]), x_bitsize=self.N_bitsize), e=ew, x=x)
            ew = bb.add(SemiclassicalQFT(e_bitsize, self.ew, i), m=ew)
        return ew, x

class SemiclassicalQFT(Bloq):
    """
        n: total number of qubits
        m: number of measured qubits in each part
        i: i-th part of the semiclassical QFT
        Note: user must explicitly recycle qubits or call with_recycling()
    """
    def __init__(self, n: int, m: int, i: int):
        self.n = n
        self.m = m
        self.i = i
        
    @property
    def signature(self) -> 'Signature':
        return Signature([Register('m', QUInt(self.m))])

    def build_composite_bloq(self, bb: 'BloqBuilder', m: 'Soquet') -> Dict[str, 'SoquetT']:
        # add Rz(m_previous), Hadamard, measure Z outcome m
        from qualtran.bloqs.basic_gates import Identity
        m = bb.add(Identity(self.m))
        return {'m', m}

    def with_recycling(self):
        return NotImplemented

class TimesExpMod(Bloq):
    """
        e,x -> e,x*g^e mod N for classical g,N
    """

    def __init__(self, g:int, N:int, e_bitsize:int, x_bitsize:int):
        self.g = g
        self.N = N
        self.e_bitsize = e_bitsize
        self.x_bitsize = x_bitsize

    @property
    def signature(self) -> 'Signature':
        return Signature([Register('e', QUInt(self.e_bitsize)), Register('x', QUInt(self.x_bitsize))])

    def build_composite_bloq(self, bb: 'BloqBuilder', e: 'Soquet', x: 'Soquet') -> Dict[str, 'SoquetT']:
        e = bb.split(e)
        g = self.g % self.N
        for j in range(self.e_bitsize - 1, 0 - 1, -1):
            e[j], x = bb.add(CModMulK(QUInt(self.x_bitsize), k=g, mod=self.N), ctrl=e[j], x=x)
            g = (g * g) % self.N

        return {'e': bb.join(e, dtype=QUInt(self.e_bitsize)), 'x': x}

if __name__ == "__main__":
    from qualtran.drawing import show_bloq
    test = RSA_GE(N=89*97, g=7, ew=3, mw=2).decompose_bloq()
    show_bloq(test)