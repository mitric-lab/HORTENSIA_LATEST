#!/usr/bin/env python3

from configparser import ConfigParser

import numpy as np
from scipy import special as ssp

import hortensia_latest.calcD as calcD

config = ConfigParser()
config.read('config.ini')

class ExactInt:
    def __init__(self, atS, grid):

        self.setupGrid(grid)

        with open("basis") as f:
            basis  = f.readlines()

        symTol   = {"S":0, "P":1, "D":2}
        symToml  = {"S":1, "P":3, "D":6}
        lTolx    = {0: [0], 1: [1,0,0], 2: [2,0,0,1,1,0]}
        lToly    = {0: [0], 1: [0,1,0], 2: [0,2,0,1,0,1]}
        lTolz    = {0: [0], 1: [0,0,1], 2: [0,0,2,0,1,1]}

        llist    = []
        atlist   = []
        explist1 = []
        explist2 = []
        tmpexp1  = []
        tmpexp2  = []
        tmpl     = []
        first    = True
        for i,j in enumerate(basis):
            if len(j.split()) == 3:
                tmpl.append(symTol[j.split()[0]])
                tmpexp1.append(float(basis[i+1].split()[0]))

                for k in range(symToml[j.split()[0]]):
                    tmpexp2.append(float(basis[i+1].split()[0]))

            if len(j.split()) == 2 and j.split()[1] == "0":
                atlist.append(basis[i].split()[0])

                if not first:
                    llist.append(tmpl)
                    explist1.append(tmpexp1)
                    explist2.append(tmpexp2)
                    tmpexp1, tmpexp2, tmpl = [], [], []

                first = False

        llist.append(tmpl)
        explist1.append(tmpexp1)
        explist2.append(tmpexp2)

        self.expalt  = []
        self.powers  = [[],[],[]]
        self.xyzlist = []

        nxyz    = 0
        for i in range(len(atS)):
            atom = atlist.index(atS[i])

            for j in explist2[atom]:
                self.expalt.append(j)
                self.xyzlist.append(nxyz)

            for j in llist[atom]:
                for k in lTolx[j]: self.powers[0].append(k)
                for k in lToly[j]: self.powers[1].append(k)
                for k in lTolz[j]: self.powers[2].append(k)

            nxyz += 1

        self.expalt = np.asarray(self.expalt)
        self.powers = np.asarray(self.powers)
        self.nbfs   = len(self.expalt)
        self.maxl   = np.max(self.powers)

        norms       = 1./np.sqrt(self.expalt**sum(self.powers) *
                                 ssp.factorial2(2*self.powers[0]-1) *
                                 ssp.factorial2(2*self.powers[1]-1) *
                                 ssp.factorial2(2*self.powers[2]-1) )

        self.cshift = np.zeros(np.shape(self.powers), dtype=int)
        self.Nshift = np.ones(len(self.powers[0]))

        for i in range(len(self.powers[-1])):
            if self.powers[0,i] == 2:
                self.cshift[0,i] = 1
                self.Nshift[i]   = norms[i]

            elif self.powers[1,i] == 2:
                self.cshift[1,i] = -1
                self.Nshift[i]   = norms[i]

            elif self.powers[2,i] == 2:
                self.cshift[2,i] = -2
                self.Nshift[i]   = norms[i]


        pk       = np.arange(self.kpoints)
        px       = np.arange(3)
        pb       = np.arange(self.nbfs)
        pI       = np.arange(3*self.maxl+1)
        pJ       = np.arange(self.maxl+1)

        self.ccd = np.sqrt(2 * self.expalt[:,None] * self.expalt[None,:] /
                           (self.expalt[:,None] + self.expalt[None,:]))
        bcd      = self.expalt[:,None] / \
                         (self.expalt[:,None] + self.expalt[None,:])
        self.b   = (self.expalt[None,:,None] + self.expalt[None,None,:]) * \
                    self.expalt[:,None,None] / (self.expalt[:,None,None] + \
                    self.expalt[None,:,None] + self.expalt[None,None,:])

        Llen   = 3*self.maxl+1
        jTok   = np.ones((Llen,Llen))
        self.jTok = np.tril(jTok)
        self.jTok[0] = 0
        self.jTok[:,0] = 0


        Ltoj = np.zeros((Llen,Llen), dtype=int)
        for lx in range(Llen):
            for lyz in range(Llen-lx):
                Ltoj[lx,lyz] = 1


        JtoS = np.zeros((3,3,3,7,3,3))
        for l1 in range(3):
            for l3 in range(3):
                for l4 in range(3):
                    for I in range(7):
                        for J in range(l1+1):
                            if (max(0,I-J-l4) < min(I,l3)+1) and (max(0,I-J-l4) < 3):
                                JtoS[l1,l3,l4,I,J,max(0,I-J-l4):min(I,l3)+1] = 1

        a,x,c,d      = np.meshgrid(pb,px,pb,pb)
        JtoS         = JtoS[self.powers[x,a],self.powers[x,c],self.powers[x,d]]


        ci,i,di,I,J,s = np.meshgrid(pb, px, pb, pI, pJ, pJ)
        f            = ssp.binom(self.powers[i,ci],s) * \
                       ssp.binom(self.powers[i,di],I-J-s) * \
                       bcd[ci,di]**(self.powers[i,ci]-s) * \
                       (-bcd[di,ci])**(self.powers[i,di]+s-I+J)
        self.f       = (np.einsum("icdIJs,iacdIJs->iacdIJ", f, JtoS))

        x,k,a,J      = np.meshgrid(px, pk, pb, pJ)
        self.Hxak    = (-1)**self.powers[x,a] * (-1j * self.kbasis[k,x])**J * \
                       ssp.binom(self.powers[x,a],J)


        self.I1bPre  = np.exp(-bcd[pb[:,None],pb[None,:]] * \
                              self.expalt[pb[None,:]])
        c,a,d        = np.meshgrid(pb, pb, pb)
        self.I1pre   = norms[a] * norms[c] * norms[d] * \
                       (2*self.expalt[c]*self.expalt[d] / \
                       (np.pi*self.expalt[a]))**0.75 / \
                       (self.expalt[c]+self.expalt[d])**1.5


        k,L,I        = np.meshgrid(pI,pI,pI)
        B1           = (-1)**k * ssp.factorial(L-I)
        Bf           = ssp.factorial(L-I-k) * ssp.factorial(2*k-L+I)
        Bf[Bf==0]    = 1
        B2           = 1 / (Bf * 2.**(L-I-k))
        B            = B1 * B2

        a,x,c,d      = np.meshgrid(pb, px, pb, pb)
        powtmp       = self.powers[x,a] + self.powers[x,c] + self.powers[x,d]
        self.B       = B[powtmp,:,:]


        c,a,d,j      = np.meshgrid(pb, pb, pb, pI)
        self.ksumPre = (2*self.b[a,c,d])**j / ssp.factorial2(2*j-1)
        self.Gpre1   = ssp.factorial2(2*j-1)
        self.Gpre2   = np.sqrt(np.pi*self.b[a,c,d])


        lB, lA       = np.meshgrid(pI, pI)
        self.lApB    = (lA + lB) * Ltoj


    def setupGrid(self, grid):
        q   = (1+np.sqrt(5))/2
        ilen  = np.sqrt(1+q**2)

        kdict = {2 : np.array([[ 1.0, 0.0, 0.0]]),
                 6 : np.array([[ 0.0, 0.0, 1.0],[ 0.0, 1.0, 0.0],
                               [ 1.0, 0.0, 0.0]]),
                 8 : np.array([[ 1.0, 1.0, 1.0],[ 1.0, 1.0,-1.0],
                               [ 1.0,-1.0, 1.0],[ 1.0,-1.0,-1.0]]) / np.sqrt(3),
                 12: np.array([[ 0.0, 1.0, q],[ 0.0, 1.0,-q],
                               [ 1.0, q, 0.0],[ 1.0,-q, 0.0],
                               [ q, 0.0, 1.0],[ q, 0.0,-1.0]]) / ilen}

        sdict = {12: np.array([[np.arctan(1/q)        , np.pi/2               ],
                               [np.pi - np.arctan(1/q), np.pi/2               ],
                               [np.pi/2               , np.arctan(q)          ],
                               [np.pi/2               , 2*np.pi - np.arctan(q)],
                               [np.arctan(q)          , 0.0                   ],
                               [np.arctan(q)          , np.pi                 ],
                               [np.pi - np.arctan(1/q), 3/2*np.pi             ],
                               [np.arctan(1/q)        , 3/2*np.pi             ],
                               [np.pi/2               , np.pi + np.arctan(q)  ],
                               [np.pi/2               , -np.arctan(q) + np.pi ],
                               [np.pi - np.arctan(q)  , np.pi                 ],
                               [np.pi - np.arctan(q)  , 0.0]])}

        try:
            self.intNk  = int(config['2e-Integrals']['intNk'])
            points = kdict[self.intNk]
            if self.intNk == 12:
                self.sph = sdict[self.intNk]
        except:
            quit("Error: Number of intNk not supported! (possible: 2,6,8,12)")

        self.nEk = int(config['Continuum']['nEk'])
        self.nkfib = int(config['Continuum']['nkfib'])
        self.intkSkip = int(config['2e-Integrals']['intkSkip'])

        if self.nEk < self.intkSkip:
            if grid.kl[0] == grid.kl[-1]:
                kl = [grid.kl[0]]
            else:
                kl = [grid.kl[0], grid.kl[-1]]
        else:
            kl = []
            for i,j in enumerate(grid.kl[::self.nkfib]):
                if i%self.intkSkip == 0:
                    kl.append(j)
            if (self.nEk-1)%self.intkSkip != 0:
                kl.append(grid.kl[-1])

        self.kpoints = len(kl) * self.intNk//2
        self.kbasis  = np.zeros((self.kpoints,3))
        self.kl      = np.zeros(self.kpoints)

        for i in range(len(kl)):
            self.kbasis[i*self.intNk//2:(i+1)*self.intNk//2] = points*kl[i]
            self.kl[i*self.intNk//2:(i+1)*self.intNk//2]     = kl[i]
        self.klsingle = kl

        return self.kbasis


    def twoElInt(self, coord):
        bfscoord = np.zeros((self.nbfs,3))
        for i in range(self.nbfs):
            bfscoord[i] = coord[self.xyzlist[i]]

        ### Calculation of D
        D, C, Cc = calcD.calcD(self.nbfs, bfscoord, self.expalt, self.kbasis,
                           self.powers, self.ccd, self.Hxak, self.f, self.B)
        Dx = D[0]
        Dy = D[1]
        Dz = D[2]

        ### Calculation of pre-sumation factor of integral
        I1a = np.exp(-self.kl[:,None]**2/(4*self.expalt[None,:]) - \
              1j*np.einsum("kx,ax->ka", self.kbasis, bfscoord))

        Rdc = bfscoord[None,:,:] - bfscoord[:,None,:]

        I1b = self.I1bPre**np.einsum("cdx,cdx->cd",Rdc,Rdc)

        self.I1 = self.I1pre[None,:,:,:]*I1a[:,:,None,None]*I1b[None,None,:,:]

        ### Generation of integral and additional terms for Cartesian d-orbitals
        j       = np.arange(3*self.maxl+1)
        Cc2     = Cc[:,:,:,:,None]**(2*j[None,None,None,None,:])
        erfTerm = ssp.erf(np.sqrt(self.b[None,:,:,:]) * Cc[:,:,:,:])
        eTerm   = np.exp(-self.b[None,:,:,:] * Cc[:,:,:,:]**2) / \
                   self.Gpre2[None,:,:,:,0]
        ksum    = self.ksumPre[None,:,:,:,:] * Cc2 / Cc[:,:,:,:,None]
        ksum    = np.einsum("Kacdk,jk->Kacdj", ksum, self.jTok)

        Int = calcD.calcInt(Dx, Dy, Dz, self.I1,
                                 self.powers, self.cshift, self.Nshift,
                                 self.nbfs, self.maxl, self.b, Cc, ksum,
                                 self.Gpre1, erfTerm, eTerm, Cc2)

        I = np.zeros((2*self.kpoints,self.nbfs,self.nbfs,self.nbfs),
                          dtype=complex)

        for ki in range(2*self.kpoints//self.intNk):
            I[ki*self.intNk:ki*self.intNk+self.intNk//2] = \
                Int[ki*self.intNk//2:(ki+1)*self.intNk//2]
            I[ki*self.intNk+self.intNk//2:(ki+1)*self.intNk] = \
                np.conjugate(Int[ki*self.intNk//2:(ki+1)*self.intNk//2])

        return I


    def interpolate(self, grid, Int):
        exlen = len(self.klsingle)
        fiblen = grid.kl[::self.nkfib]
        lfilen = len(fiblen)
        fibint = np.zeros(len(grid.kl), dtype=complex)

        betwI  = []
        betwH  = []
        tmplen = 0.0
        for l in range(exlen):
            theta = self.sph[:,0]
            phi   = self.sph[:,1]

            ttmp  = []
            ptmp  = []
            iRtmp = []
            iItmp = []
            for i in range(len(theta)):
                ii = l*self.intNk + i
                if theta[i] < np.pi/2:
                    if phi[i] < np.pi:
                        ttmp.append(theta[i])
                        ptmp.append(phi[i]+2*np.pi)
                        ttmp.append(-theta[i])
                        ptmp.append(phi[i]+2*np.pi)
                    elif phi[i] > np.pi:
                        ttmp.append(theta[i])
                        ptmp.append(phi[i]-2*np.pi)
                        ttmp.append(-theta[i])
                        ptmp.append(phi[i]-2*np.pi)
                    else:
                        continue
                    ttmp.append(-theta[i])
                    ptmp.append(phi[i])

                    for _ in range(3):
                        iRtmp.append(Int[ii].real)
                        iItmp.append(Int[ii].imag)

                elif theta[i] > np.pi/2:
                    if phi[i] < np.pi:
                        ttmp.append(theta[i])
                        ptmp.append(phi[i]+2*np.pi)
                        ttmp.append(2*np.pi-theta[i])
                        ptmp.append(phi[i]+2*np.pi)
                    elif phi[i] > np.pi:
                        ttmp.append(theta[i])
                        ptmp.append(phi[i]-2*np.pi)
                        ttmp.append(2*np.pi-theta[i])
                        ptmp.append(phi[i]-2*np.pi)
                    else:
                        continue
                    ttmp.append(2*np.pi-theta[i])
                    ptmp.append(phi[i])

                    for _ in range(3):
                        iRtmp.append(Int[ii].real)
                        iItmp.append(Int[ii].imag)

            tnew = np.zeros(len(ttmp)+len(theta))
            pnew = np.zeros(len(ptmp)+len(phi))
            iRnew = np.zeros(len(iRtmp)+len(theta))
            iInew = np.zeros(len(iItmp)+len(theta))
            tnew[:len(ttmp)] = np.asarray(ttmp)
            tnew[len(ttmp):] = theta
            pnew[:len(ttmp)] = np.asarray(ptmp)
            pnew[len(ttmp):] = phi
            iRnew[:len(ttmp)] = np.asarray(iRtmp)
            iRnew[len(ttmp):] = Int[l*self.intNk:(l+1)*self.intNk].real
            iInew[:len(ttmp)] = np.asarray(iItmp)
            iInew[len(ttmp):] = Int[l*self.intNk:(l+1)*self.intNk].imag

            from scipy.interpolate import griddata
            testR = griddata((tnew, pnew), iRnew, (grid.theta, grid.phi),
                             method="linear")
            testI = griddata((tnew, pnew), iInew, (grid.theta, grid.phi),
                             method="linear")
            #print(testR)

            for i in range(lfilen):
                if fiblen[i] == self.klsingle[l]:
                    fibint[i*self.nkfib:(i+1)*self.nkfib] = testR + 1.0j*testI
                    betwH.append(i)
                elif fiblen[i] > tmplen and fiblen[i] < self.klsingle[l]:
                    betwI.append([i,l])
                elif fiblen[i] > self.klsingle[l]:
                    tmplen = self.klsingle[l]
                    break

        for i in betwI:
            for j in range(1,len(betwH)):
                if i[0] > betwH[j-1] and i[0] < betwH[j]:

                    Ilow  = fibint[betwH[j-1] * \
                            self.nkfib:(betwH[j-1]+1)*self.nkfib]
                    Ihigh = fibint[betwH[j] * \
                            self.nkfib:(betwH[j]+1)*self.nkfib]

                    fact  = (fiblen[i[0]]-fiblen[betwH[j-1]]) / \
                            (fiblen[betwH[j]]-fiblen[betwH[j-1]])
                    fibint[i[0]*self.nkfib:(i[0]+1)*self.nkfib] = Ilow + \
                        (Ihigh - Ilow) * fact

        return fibint
