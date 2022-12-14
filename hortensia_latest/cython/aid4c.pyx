cimport cython

import numpy as np
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[complex,ndim=1] int2eParallel(np.ndarray[double, ndim=3] A,
    np.ndarray[double, ndim=3] Ap, np.ndarray[complex, ndim=2] B,
    np.ndarray[double, ndim=1] nrm, int nbfs, mol, np.ndarray[int, ndim=1] ii,
    int j):
        '''
        function for parallel calculation of four-center integrals
        '''

        cdef int klen = len(B[0])

        cdef np.ndarray[double, ndim=4] int4c

        cdef double tmp3m4
        cdef np.ndarray[complex, ndim=1] term3m4 = np.zeros(klen, dtype=complex)

        cdef int lenl1 = 0
        cdef int l1 = 0
        cdef int l2 = 0
        cdef int l3 = 0
        cdef int l4 = 0
        cdef int k  = 0

        cdef int i = 0

        for i in range(ii[0],ii[1]):
            int4c  = mol.intor('int2e_cart', shls_slice = (
                               i, i+1, 0, mol.nbas, 0, mol.nbas, 0, mol.nbas))

            lenl1  = len(int4c)

            for l1 in range(lenl1):
                for l2 in range(nbfs):
                    for l3 in range(nbfs):
                        for l4 in range(nbfs):
                            int4c[l1,l2,l3,l4] = int4c[l1,l2,l3,l4] * \
                                               nrm[j+l1]*nrm[l2]*nrm[l3]*nrm[l4]

                tmp3m4 = 0.0
                for l2 in range(nbfs):
                    for l3 in range(nbfs):
                        for l4 in range(nbfs):
                            tmp3m4 += int4c[l1,l3,l2,l4] * \
                                      (A[l2,l3,l4] + Ap[l2,l3,l4]) - \
                                      int4c[l1,l4,l2,l3] * \
                                      A[l2,l3,l4]

                for k in range(klen):
                    term3m4[k] = term3m4[k] + B[j+l1,k] * tmp3m4

            j += lenl1

        return term3m4

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[complex, ndim=2] calcSkMOa(np.ndarray[complex, ndim=2] preSkMO,
    np.ndarray[double, ndim=2] orbs, int nk, int nbfs, int nmo):
        cdef np.ndarray[complex, ndim=2] SkMOa = np.zeros((nmo, nk),
                                                          dtype=complex)

        for m in range(nmo):
            for a in range(nbfs):
                for k in range(nk):
                    SkMOa[m, k] = SkMOa[m, k] + preSkMO[k, a] * orbs[a, m]

        return SkMOa
