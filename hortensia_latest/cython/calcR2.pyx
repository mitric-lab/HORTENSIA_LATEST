cimport cython

import numpy as np
cimport numpy as np

cpdef double getR2(int nbas, np.ndarray[double, ndim=1] C, 
    np.ndarray[long, ndim=1] contl, np.ndarray[double, ndim=1] pl, 
    np.ndarray[double, ndim=1] nrm, np.ndarray[complex, ndim=1] bl, 
    np.ndarray[long, ndim=2] Ll, np.ndarray[complex, ndim=1] at, 
    np.ndarray[double, ndim=2] Rl):
        cdef np.ndarray[complex, ndim=1] blc = np.conjugate(bl)

        cdef np.ndarray[complex, ndim=2] ap = np.conjugate(at)[:, None] + \
                                              at[None, :]
        cdef np.ndarray[double, ndim=3] AmB = Rl[:, None, :] - Rl[None, :, :]
        cdef np.ndarray[complex, ndim=2] EIJ = \
            np.exp(-np.conjugate(at)[:, None] * at[None, :] / ap * \
            np.einsum("ijk,ijk->ij", AmB, AmB))

        cdef np.ndarray[complex, ndim=3] P = \
            (np.conjugate(at)[:, None, None] * Rl[:, None, :] + 
             at[None, :, None] * Rl[None, :, :]) / \
            (np.conjugate(at)[:, None, None] + at[None, :, None])

        cdef double r2 = 0.0
        cdef complex IR2 = 0.0 + 0.0j

        cdef int i = 0
        cdef int j = 0
        for i in range(nbas):
            r2 = r2 + abs(C[contl[i]] * pl[i] * nrm[i] * bl[i])**2 * \
                Iabr2(Ll[i], Ll[i], ap[i, i], Rl[i], Rl[i], 
                      EIJ[i, i], P[i, i]).real

            for j in range(i+1, nbas):
                # value for comparison: nrm[i]*nrm[j]*intr2[i,j]
                IR2 = C[contl[i]] * C[contl[j]] * pl[i] * pl[j] * \
                    nrm[i] * nrm[j] * blc[i] * bl[j] * \
                    Iabr2(Ll[i], Ll[j], ap[i, j], Rl[i], Rl[j], 
                          EIJ[i, j], P[i, j])
                r2 = r2 + 2 * IR2.real
        
        return r2

@cython.boundscheck(False)
@cython.wraparound(False)
cdef complex Iabr2(np.ndarray[long, ndim=1] l1, np.ndarray[long, ndim=1] l2, 
    complex ap, np.ndarray[double, ndim=1] A, np.ndarray[double, ndim=1] B, 
    complex EIJ, np.ndarray[complex, ndim=1] P):
        cdef complex retval = 0.0 + 0.0j
        cdef complex val = 0.0 + 0.0j

        cdef double pi = 3.141592653589793

        cdef int x = 0
        cdef int y = 0 
        cdef int z = 0

        for x in range(3):
            y, z = (x+1)%3, (x+2)%3
            val = 0.0 + 0.0j
            # for x2, we only need L=0, M=0 and N = 0, 1, 2, 
            # thus we sum explicitly

            # N = L = M == 0
            val = val + (P[x]**2 + 1/(2*ap)) * \
                  get_d(0, l1[x], l2[x], A[x], B[x], P[x], ap)
            # N = 1, L = M = 0
            if (l1[x] + l2[x]) > 0:
                val = val + 2 * P[x] * \
                    get_d(1, l1[x], l2[x], A[x], B[x], P[x], ap)
                # N = 2, L = M = 0
                if (l1[x] + l2[x]) > 1:
                    val = val + 2 * \
                        get_d(2, l1[x], l2[x], A[x], B[x], P[x], ap)

            val = val * (pi/ap)**(1.5) * \
                    get_d(0, l1[y], l2[y], A[y], B[y], P[y], ap) * \
                    get_d(0, l1[z], l2[z], A[z], B[z], P[z], ap)

            retval = retval + EIJ * val

        return retval # UNNORMALIZED!!

@cython.boundscheck(False)
@cython.wraparound(False)
cdef complex get_d(int N, int n1, int n2, double A, double B, complex P, 
    complex ap):
        # A, B, P should be the relevant components of the respective vectors
        cdef complex PA = P-A
        cdef complex PB = P-B
        if (N < 0) or (N > (n1 + n2)): 
            return 0.0

        elif (N == n1 == n2 == 0): 
            return 1.0

        elif (n1 == 0):
            retval = 1/(2*ap) * \
                get_d(N-1, n1, n2-1, A, B, P, ap) + \
                PB * \
                get_d(N, n1, n2-1, A, B, P, ap) + \
                (N + 1) * \
                get_d(N+1, n1, n2-1, A, B, P, ap)

            return retval

        else:
            retval = 1/(2*ap) * \
                get_d(N-1, n1-1, n2, A, B, P, ap) + \
                PA * \
                get_d(N, n1-1, n2, A, B, P, ap) + \
                (N + 1) * \
                get_d(N+1, n1-1, n2, A, B, P, ap)

            return retval