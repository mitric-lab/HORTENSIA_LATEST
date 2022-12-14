cimport cython

import numpy as np
cimport numpy as np
import cmath
import scipy.special as ssp
from scipy.interpolate import griddata

import hortensia_latest.misc as misc

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calcD(int nbfs, np.ndarray[double, ndim=2] coord,
    np.ndarray[double, ndim=1] exp, np.ndarray[double, ndim=2] kbasis,
    np.ndarray[long, ndim=2] powers, np.ndarray[double, ndim=2] ccd,
    np.ndarray[complex, ndim=4] Hxak, np.ndarray[double, ndim=6] f,
    np.ndarray[double, ndim=6] B):

        cdef int maxl = np.max(powers)
        cdef int klen = len(Hxak)
        cdef int Ilen = 3*maxl+1
        cdef int Jlen = maxl+1
        pb   = np.arange(nbfs)
        pI   = np.arange(Ilen)
        pJ   = np.arange(Jlen)

        cdef np.ndarray[complex, ndim=6] D  = np.zeros((3,klen,nbfs,nbfs,nbfs,
                                                        Ilen),dtype=complex)
        cdef np.ndarray[double , ndim=4] A1 = np.zeros((nbfs,nbfs,Ilen,Jlen))
        cdef complex A

        cdef complex tmp
        cdef int tmpI
        cdef int i
        cdef int ki
        cdef int ai
        cdef int ci
        cdef int di
        cdef int li
        cdef int Ii
        cdef int Ji

        cdef double t0
        cdef double t1
        cdef double t2
        cdef double t3
        cdef double t4

        cdef np.ndarray[complex,ndim=5] C  = np.zeros((klen,3,nbfs,nbfs,nbfs),
                                                      dtype=complex)

        d,c,I,J    = np.meshgrid(pb,pb,pI,pJ)
        ImJ        = I-J
        ImJ[ImJ<0] = 0

        for i in range(3):
            A1 = (-1)**(powers[i,d]+J) * ccd[c,d]**(I-J) * \
                 ssp.eval_hermitenorm(ImJ, ccd[c,d] *
                 (coord[d,i] - coord[c,i]))

            for ki in range(klen):
                for ai in range(nbfs):
                    for ci in range(nbfs):
                        for di in range(ci,nbfs):
                            C[ki,i,ai,ci,di] = (exp[ci] * coord[ci,i] + \
                                exp[di] * coord[di,i]) / (exp[ci] + exp[di]) - \
                                coord[ai,i] + 1.j*kbasis[ki,i]/(2*exp[ai])
                            C[ki,i,ai,di,ci] = C[ki,i,ai,ci,di]

                            tmpI = powers[i,ai] + powers[i,ci] + powers[i,di]

                            for li in range(tmpI+1):
                                tmp = 0.0+0.0j
                                for Ii in range(max(0,tmpI-2*li), tmpI-li+1):
                                    A = 0.0 + 0.0j

                                    for Ji in range(max(0,Ii-(powers[i,ci]+powers[i,di])),
                                        min(Ii,powers[i,ai])+1):
                                            A = A + A1[ci,di,Ii,Ji] * \
                                                Hxak[ki,i,ai,Ji] * \
                                                f[i,ai,ci,di,Ii,Ji]

                                    tmp = tmp + A * B[i,ai,ci,di,li,Ii] * \
                                           C[ki,i,ai,ci,di]**(2*li+Ii-tmpI)

                                D[i,ki,ai,ci,di,li] = tmp
                                D[i,ki,ai,di,ci,li] = tmp

        Cc = np.sqrt(np.einsum("kxacd,kxacd->kacd",C,C))

        return D, C, Cc


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[complex,ndim=4] calcInt(np.ndarray[complex, ndim=5] Dx,
    np.ndarray[complex, ndim=5] Dy, np.ndarray[complex, ndim=5] Dz,
    np.ndarray[complex, ndim=4] I1, np.ndarray[long, ndim=2] powers,
    np.ndarray[long, ndim=2] shift, np.ndarray[double, ndim=1] norms,
    int nbfs,
    int maxl, np.ndarray[double, ndim=3] b, np.ndarray[complex, ndim=4] Cc,
    np.ndarray[complex, ndim=5] ksum, np.ndarray[double, ndim=4] Gpre1,
    np.ndarray[complex, ndim=4] erf, np.ndarray[complex, ndim=4] eTerm,
    np.ndarray[complex, ndim=5] Cc2):

        cdef int klen = len(Dx)
        cdef np.ndarray[complex,ndim=4] Int  = np.zeros((klen,nbfs,nbfs,nbfs),
                                                         dtype=complex)

        cdef np.ndarray[long, ndim=1] shiftsum = sum(shift)
        cdef int nd
        cdef complex tmp2 = 0.0+0.0j

        cdef float sqrt3 = np.sqrt(3)

        cdef int ki
        cdef int ai
        cdef int ci
        cdef int di
        cdef int i
        cdef int lx
        cdef int ly
        cdef int lz
        cdef int Lx
        cdef int Ly
        cdef int Lz
        cdef int ai0
        cdef int ai1
        cdef int ai2
        cdef int ci0
        cdef int ci1
        cdef int ci2
        cdef int di0
        cdef int di1
        cdef int di2

        cdef int j        = 0
        cdef int k        = 0
        cdef int maxj     = 0
        cdef complex[7] G = [0.0+0.0j, 0.0+0.0j, 0.0+0.0j, 0.0+0.0j,
                              0.0+0.0j, 0.0+0.0j, 0.0+0.0j]
        cdef np.ndarray[double, ndim=3] bpisqrt = np.sqrt(b/np.pi)

        for ki in range(klen):
            for ai in range(nbfs):
                for ci in range(nbfs):
                    for di in range(ci,nbfs):
                        Lx = powers[0,ai] + powers[0,ci] + powers[0,di]
                        Ly = powers[1,ai] + powers[1,ci] + powers[1,di]
                        Lz = powers[2,ai] + powers[2,ci] + powers[2,di]

                        maxj = Lx+Ly+Lz

                        if maxj == 0:
                            G[0] = erf[ki,ai,ci,di] / Cc[ki,ai,ci,di]
                        else:
                            for j in range(maxj+1):
                                if abs(Cc2[ki,ai,ci,di,j]) < 1e-10:
                                    G[j] = bpisqrt[ai,ci,di] * \
                                        (2*b[ai,ci,di])**j/(j+0.5)
                                else:
                                    G[j] = Gpre1[ai,ci,di,j]
                                    G[j] = G[j] / Cc2[ki,ai,ci,di,j] / Cc[ki,ai,ci,di]
                                    G[j] = G[j] * (erf[ki,ai,ci,di] - \
                                        eTerm[ki,ai,ci,di] * ksum[ki,ai,ci,di,j])

                        tmp2 = 0.0+0.0j
                        for lx in range(Lx+1):
                            for ly in range(Ly+1):
                                for lz in range(Lz+1):
                                    tmp2 += G[lx+ly+lz] * \
                                            Dx[ki,ai,ci,di,lx] * \
                                            Dy[ki,ai,ci,di,ly] * \
                                            Dz[ki,ai,ci,di,lz]

                        Int[ki,ai,ci,di] = tmp2 * I1[ki,ai,ci,di]

                        nd = 0
                        if shiftsum[ai] != 0:
                            nd += 1
                        if shiftsum[ci] != 0:
                            nd += 1
                        if shiftsum[di] != 0:
                            nd += 1

                        if nd == 0:
                            Int[ki,ai,di,ci] = Int[ki,ai,ci,di]
                            continue

                        ai0,ai1,ai2 = ai+shift[0,ai], ai+shift[1,ai], ai+shift[2,ai]
                        ci0,ci1,ci2 = ci+shift[0,ci], ci+shift[1,ci], ci+shift[2,ci]
                        di0,di1,di2 = di+shift[0,di], di+shift[1,di], di+shift[2,di]

                        Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                        Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                        Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                        tmp2 = 0.0+0.0j
                        for lx in range(Lx+1):
                            for ly in range(Ly+1):
                                for lz in range(Lz+1):
                                    tmp2 += G[lx+ly+lz] * \
                                            Dx[ki,ai0,ci0,di0,lx] * \
                                            Dy[ki,ai1,ci1,di1,ly] * \
                                            Dz[ki,ai2,ci2,di2,lz]

                        if nd == 1:
                            Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                2/3.**0.5 * I1[ki,ai,ci,di] * \
                                tmp2 / (norms[ai] * norms[ci] * norms[di])
                            Int[ki,ai,di,ci] = Int[ki,ai,ci,di]
                            continue

                        elif nd == 2:
                            Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                4/3. * I1[ki,ai,ci,di] * \
                                tmp2 / (norms[ai] * norms[ci] * norms[di])

                        elif nd == 3:
                            Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                8/27.**0.5 * I1[ki,ai,ci,di] * \
                                tmp2 / (norms[ai] * norms[ci] * norms[di])

                        if shiftsum[ai] != 0:
                            ai0,ai1,ai2 = ai, ai, ai
                            ci0,ci1,ci2 = ci+shift[0,ci], ci+shift[1,ci], ci+shift[2,ci]
                            di0,di1,di2 = di+shift[0,di], di+shift[1,di], di+shift[2,di]

                            Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                            Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                            Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                            tmp2 = 0.0+0.0j
                            for lx in range(Lx+1):
                                for ly in range(Ly+1):
                                    for lz in range(Lz+1):
                                        tmp2 += G[lx+ly+lz] * \
                                                Dx[ki,ai0,ci0,di0,lx] * \
                                                Dy[ki,ai1,ci1,di1,ly] * \
                                                Dz[ki,ai2,ci2,di2,lz]
                            if nd == 2:
                                Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                    2/3.**0.5 * I1[ki,ai,ci,di] * \
                                    tmp2 / (norms[ci] * norms[di])
                            elif nd == 3:
                                Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                    4/3. * I1[ki,ai,ci,di] * \
                                    tmp2 / (norms[ci] * norms[di])

                        if shiftsum[ci] != 0:
                            ai0,ai1,ai2 = ai+shift[0,ai], ai+shift[1,ai], ai+shift[2,ai]
                            ci0,ci1,ci2 = ci, ci, ci
                            di0,di1,di2 = di+shift[0,di], di+shift[1,di], di+shift[2,di]

                            Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                            Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                            Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                            tmp2 = 0.0+0.0j
                            for lx in range(Lx+1):
                                for ly in range(Ly+1):
                                    for lz in range(Lz+1):
                                        tmp2 += G[lx+ly+lz] * \
                                                Dx[ki,ai0,ci0,di0,lx] * \
                                                Dy[ki,ai1,ci1,di1,ly] * \
                                                Dz[ki,ai2,ci2,di2,lz]
                            if nd == 2:
                                Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                    2/3.**0.5 * I1[ki,ai,ci,di] * \
                                    tmp2 / (norms[ai] * norms[di])
                            elif nd == 3:
                                Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                    4/3. * I1[ki,ai,ci,di] * \
                                    tmp2 / (norms[ai] * norms[di])

                        if shiftsum[di] != 0:
                            ai0,ai1,ai2 = ai+shift[0,ai], ai+shift[1,ai], ai+shift[2,ai]
                            ci0,ci1,ci2 = ci+shift[0,ci], ci+shift[1,ci], ci+shift[2,ci]
                            di0,di1,di2 = di, di, di

                            Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                            Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                            Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                            tmp2 = 0.0+0.0j
                            for lx in range(Lx+1):
                                for ly in range(Ly+1):
                                    for lz in range(Lz+1):
                                        tmp2 += G[lx+ly+lz] * \
                                                Dx[ki,ai0,ci0,di0,lx] * \
                                                Dy[ki,ai1,ci1,di1,ly] * \
                                                Dz[ki,ai2,ci2,di2,lz]
                            if nd == 2:
                                Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                    2/3.**0.5 * I1[ki,ai,ci,di] * \
                                    tmp2 / (norms[ai] * norms[ci])

                            elif nd == 3:
                                Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                    4/3. * I1[ki,ai,ci,di] * \
                                    tmp2 / (norms[ai] * norms[ci])

                        if nd == 2:
                            Int[ki,ai,di,ci] = Int[ki,ai,ci,di]
                            continue

                        # single d->s terms, if all orbs are d-type
                        ai0,ai1,ai2 = ai+shift[0,ai], ai+shift[1,ai], ai+shift[2,ai]
                        ci0,ci1,ci2 = ci, ci, ci
                        di0,di1,di2 = di, di, di

                        Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                        Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                        Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                        tmp2 = 0.0+0.0j
                        for lx in range(Lx+1):
                            for ly in range(Ly+1):
                                for lz in range(Lz+1):
                                    tmp2 += G[lx+ly+lz] * \
                                            Dx[ki,ai0,ci0,di0,lx] * \
                                            Dy[ki,ai1,ci1,di1,ly] * \
                                            Dz[ki,ai2,ci2,di2,lz]

                        Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                2/3.**0.5 * I1[ki,ai,ci,di] * \
                                tmp2 / norms[ai]

                        ai0,ai1,ai2 = ai, ai, ai
                        ci0,ci1,ci2 = ci+shift[0,ci], ci+shift[1,ci], ci+shift[2,ci]
                        di0,di1,di2 = di, di, di

                        Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                        Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                        Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                        tmp2 = 0.0+0.0j
                        for lx in range(Lx+1):
                            for ly in range(Ly+1):
                                for lz in range(Lz+1):
                                    tmp2 += G[lx+ly+lz] * \
                                            Dx[ki,ai0,ci0,di0,lx] * \
                                            Dy[ki,ai1,ci1,di1,ly] * \
                                            Dz[ki,ai2,ci2,di2,lz]

                        Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                2/3.**0.5 * I1[ki,ai,ci,di] * \
                                tmp2 / norms[ci]

                        ai0,ai1,ai2 = ai, ai, ai
                        ci0,ci1,ci2 = ci, ci, ci
                        di0,di1,di2 = di+shift[0,di], di+shift[1,di], di+shift[2,di]

                        Lx = powers[0,ai0] + powers[0,ci0] + powers[0,di0]
                        Ly = powers[1,ai1] + powers[1,ci1] + powers[1,di1]
                        Lz = powers[2,ai2] + powers[2,ci2] + powers[2,di2]
                        tmp2 = 0.0+0.0j
                        for lx in range(Lx+1):
                            for ly in range(Ly+1):
                                for lz in range(Lz+1):
                                    tmp2 += G[lx+ly+lz] * \
                                            Dx[ki,ai0,ci0,di0,lx] * \
                                            Dy[ki,ai1,ci1,di1,ly] * \
                                            Dz[ki,ai2,ci2,di2,lz]

                        Int[ki,ai,ci,di] = Int[ki,ai,ci,di] + \
                                2/3.**0.5 * I1[ki,ai,ci,di] * \
                                tmp2 / norms[di]
                        Int[ki,ai,di,ci] = Int[ki,ai,ci,di]

        return Int
