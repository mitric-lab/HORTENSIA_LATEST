#!/usr/bin/env python3

import numpy as np
import numpy.linalg as la
from scipy.special import factorial,factorial2

def enumerate_microstates(na, nb, norb):
    """
    enumerate the occupations for all singly excited Slater determinants
    (microstates)

    Parameters
    ----------
    na   : number of occupied alpha orbitals
    nb   : number of occupied beta orbitals, total number of electrons is na+nb
    norb : total number of molecular orbitals, occupied + virtual

    Returns
    -------
    iterator to microstates, each microstate is a list of 0's and 1's indicating
    whether the orbitals is occupied or not in the microstate. The first half
    are the occupations for the alpha orbitals, the second half are for the
    beta orbitals. Examples:
         [1, 0, 1, 0]  ground state, HOMO is double occuped
         [1, 0, 0, 1]  singly excited Slater determinant  HOMO(beta)->LUMO(beta)
    The total number of microstates (Slater determinants) is given by
         1 + na * (norb-na) + nb * (norb-nb)
    """

    # ground state
    alpha0 = [1] * na + [0] * (norb - na)
    beta0  = [1] * nb + [0] * (norb - nb)
    micro = alpha0 + beta0
    yield np.array(micro)

    # enumerate all single excitations from the GS which conserve Sz=(na-nb)
    for i in range(0, na):  # enumerate occupied orbitals
        for a in range(na, norb):  # enumerate virtual orbitals
            # promotion of an alpha orbital  i->a
            alpha = alpha0[:]
            alpha[i] = 0
            alpha[a] = 1
            micro = alpha + beta0
            yield np.array(micro)

    for i in range(0, nb):  # enumerate occupied orbitals
        for a in range(nb, norb):  # enumerate virtual orbitals
            # promotion of a beta orbital  i->a
            beta = beta0[:]
            beta[i] = 0
            beta[a] = 1
            micro = alpha0 + beta
            yield np.array(micro)


def AminusB_diagonal(orbe_alpha, orbe_beta, na, nb, norb):
    """
    compute diagonal elements of A-B matrix for a pure DFT functional
    """

    dim = 1 + (norb - na) * na + (norb - nb) * nb
    AmB = np.zeros(dim, dtype=float)
    AmB[0] = 1.0

    # alpha orbitals
    for occa in range(0, na):
        for virta in range(na, norb):
            indexA = 1 + (norb - na) * occa + virta - na
            AmB[indexA] = orbe_alpha[virta] - orbe_alpha[occa]

    # beta orbitals
    for occb in range(0, nb):
        for virtb in range(nb, norb):
            indexB = 1 + (norb - nb) * occb + virtb - nb + (norb - na) * na
            AmB[indexB] = orbe_beta[virtb] - orbe_beta[occb]

    return AmB


def basis_overlap(basis1, basis2):
    """
    compute the overlap matrix between two sets of uncontracted basis functions
    """

    l1, m1, n1 = basis1.powersD[0], basis1.powersD[1], basis1.powersD[2]
    l2, m2, n2 = basis2.powersD[0], basis2.powersD[1], basis2.powersD[2]
    alpha1     = basis1.exponentsD
    alpha2     = basis2.exponentsD
    A          = basis1.centersD
    B          = basis2.centersD
    S          = overlap(alpha1, l1, m1, n1, A, alpha2, l2, m2, n2, B)

    return S

###########################################################################
#                                                                         #
# INTEGRALS                                                               #
#                                                                         #
# routines for computing overlaps and dipole matrix elements between      #
# primitive Gaussian basis functions, stolen from 'pyints.py' in PyQuante #
#                                                                         #
# see THO paper                                                           #
#   H.Taketa, S. Huzinga, K. O-ohata,                                     #
#  'Gaussian-Expansion Methods for Molecular Integrals'                   #
#   J. Phys. Soc. Japan, 21, 2313, 1966.                                  #
#                                                                         #
###########################################################################


def norm(alpha, l, m, n):
    # eqn. (2.2) from THO paper
    if isinstance(l, np.int64):
        nrm = np.sqrt(pow(2, 2 * (l+m+n) + 1.5) *
                      pow(alpha, l+m+n + 1.5) /
                      factorial2(2*l - 1) /
                      factorial2(2*m - 1) /
                      factorial2(2*n - 1) /
                      pow(np.pi, 1.5))
    else:
        f2l = factorial2(2*l - 1)
        f2m = factorial2(2*m - 1)
        f2n = factorial2(2*n - 1)
        f2term = f2l * f2m * f2n
        nrm = (2/np.pi)**0.75 * (4**(l+m+n) * alpha**(l+m+n+1.5) / f2term)**0.5

    return nrm


def overlap(alpha1, l1, m1, n1, A, alpha2, l2, m2, n2, B):
    """
    overlap <g1|g2> between two normalized cartesian Gaussian functions which
    are centered at A and B

    g1(x,y,z) = N(l1,m1,n1) *
                (x-Ax)**l1 * (y-Ay)**m1 * (z-Az)**n1 *
                e**(-alpha1 * (r-A)**2)

    and a similar expression for g2(x,y,z).

    Parameters
    ----------
    alpha1     :  exponent, > 0
    (l1,m1,n1) :  powers, L=l1+m1+n1 is the angular momentum
    A          :  3-dim vector, center of basis function

    and similary for second basis function

    Returns
    -------
    scalar, <g1|g2>
    """

    olap = overlap_unnorm(alpha1, l1, m1, n1, A, alpha2, l2, m2, n2, B)
    norm1 = norm(alpha1, l1, m1, n1)
    norm2 = norm(alpha2, l2, m2, n2)

    return norm1[:, None] * norm2[None, :] * olap


def overlap_unnorm(alpha1, l1, m1, n1, A, alpha2, l2, m2, n2, B):
    """
    overlap without normalization constants
    """

    rab2  = np.einsum("xab->ab", (A[:, :, None] - B[:, None, :])**2)
    gamma = alpha1[:, None] + alpha2[None, :]

    P = (alpha1[None, :, None] * A[:, :, None] +
         alpha2[None, None, :] * B[:, None, :]) / gamma[None, :, :]

    pre = (np.pi/gamma)**1.5 \
        * np.exp(-alpha1[:, None] * alpha2[None, :] * rab2 / gamma)
    wx = overlap_1Darray(l1, l2, P[0], A[0], B[0], gamma)
    wy = overlap_1Darray(m1, m2, P[1], A[1], B[1], gamma)
    wz = overlap_1Darray(n1, n2, P[2], A[2], B[2], gamma)

    return pre*wx*wy*wz


def overlap_1Darray(l1, l2, Px, Ax, Bx, gamma):
    """
    overlap of 1D Gaussians for the whole array of basis functions, definitely
    programmed a lot more complicated and confusing than it needs to be,
    but fast at least
    """

    PAx = Px-Ax[:, None]
    PBx = Px-Bx[None, :]

    l1max, l2max = max(l1), max(l2)
    kmax         = int(0.5 * (l1max + l2max))
    tmax         = 2 * kmax

    lTot = np.zeros((l1max + 1, l2max + 1, kmax + 1, tmax + 1))
    for i in range(l1max + 1):
        for j in range(l2max + 1):
            for k in range(1 + int(0.5 * (i + j))):
                for t in range(1 + 2*k):
                    if 2*k - i <= t <= j:
                        lTot[i, j, k, t] = 1.0
    lTot = lTot[l1[:, None], l2[None, :], :, :]

    k = np.arange(kmax + 1)
    t = np.arange(1 + 2*kmax)

    # The following lines calculate the necessary factorials for all
    # combinations of k and t, regardless of the fact that 2k - t cannot be < 0.
    # Instances where 2*k - t < 0 produce fac(x) = 0 and therefore numerical
    # errors when divided by. Since they are sorted out by lTot anyways these
    # instances are set to 1 arbitrarily
    # Sadly in Python it is way faster to calculate everything and throw away
    # half the results than to calculate only the necessary terms
    f1 = factorial(l1)
    f2 = factorial(2*k[:, None] - t[None, :])
    f2 = f2 + np.asarray(f2==0, dtype=float)
    f3 = factorial(l1[:, None, None] - 2*k[None, :, None] + t[None, None, :])
    f3 = f3 + np.asarray(f3==0, dtype=float)
    f4 = factorial(l2)
    f5 = factorial(t)
    f6 = factorial(l2[:, None] - t[None, :])
    f6 = f6 + np.asarray(f6==0, dtype=float)

    expA = lTot * (l1[:, None, None, None] - 2*k[None, None, :, None] +
                   t[None, None, None, :])
    expB = lTot * (l2[None, :, None, None] - t[None, None, None, :])
    PAxterm = PAx[:, :, None, None]**expA
    PBxterm = PBx[:, :, None, None]**expB

    tmp1 = f1[:, None, None, None] / f2[None, None, :, :] / \
           f3[:, None, :, :] * f4[None, :, None, None] / \
           f5[None, None, None, :] / f6[None, :, None, :] * \
           PAxterm * PBxterm
    osum = np.einsum("ijkt,k,ijk,ijkt->ij", tmp1, factorial2(2*k - 1),
                     1 / (2*gamma[:, :, None])**k[None, None, :], lTot)

    return osum


def dist2(A, B):
    return pow(A[0] - B[0], 2) + pow(A[1] - B[1], 2) + pow(A[2] - B[2], 2)


def gaussian_product_center(alpha1, A, alpha2, B):
    """
    see Szabo & Ostlund
    """

    gamma = alpha1 + alpha2
    return (alpha1 * A[0] + alpha2 * B[0]) / gamma,\
           (alpha1 * A[1] + alpha2 * B[1]) / gamma,\
           (alpha1 * A[2] + alpha2 * B[2]) / gamma

################################################################################
#                                                                              #
#                             Dyson Orbitals                                   #
#                                                                              #
################################################################################


def dyson_orbitals(ionized, neutral, fo, perc_microstates=97.0, max_microstates=10):
    """
    compute Dyson orbitals as the overlap between an electronic wavefunction
    of the ionized molecule and a wavefunction of the neutral molecule

    Parameters
    ----------
    ionized: instance of G09ResultsTDDFT object with electronic structure of
             ionized molecule
    neutral: instance of G09ResultsTDDFT object with electronic structure of
             neutral molecule

    Optional
    --------
    perc_microstates : percentage 0 < perc <= 100 controlling how many
                       microstates are included
    max_microstates  : only include the first most important microstates

    Returns
    -------
    dyson_orbs: 3D numpy array with M.O. coefficients of Dyson orbitals for all
                ionization channels, dyson_orbs[:,i,n] are the coefficients for
                ionization from neutral state `n` to ionized state `i`.
    """

    # overlap between AO bases of neutral and ionized molecules
    Sao = basis_overlap(ionized.basis, neutral.basis)

    # overlap between MOs of neutral and ionized molecules
    Smo_alpha = np.dot(ionized.orbs_alpha.T, np.dot(Sao, neutral.orbs_alpha))
    Smo_beta  = np.dot(ionized.orbs_beta.T , np.dot(Sao, neutral.orbs_beta ))

    # spin overlap
    Smo = np.bmat( [[Smo_alpha, np.zeros((ionized.nmo, neutral.nmo))],
                    [np.zeros((ionized.nmo, neutral.nmo)), Smo_beta]] )
    Smo = np.array(Smo)

    nocca , noccb   = neutral.nelec_alpha, neutral.nelec_beta
    nvirta          = neutral.nmo-neutral.nelec_alpha
    nvirtb          = neutral.nmo-neutral.nelec_beta
    groundstate     = [[]]
    excitations     = groundstate + order_excs((nocca, noccb), (nvirta, nvirtb))
    excitations_ion = groundstate

    orbs     = (neutral.orbs_alpha , neutral.orbs_beta )
    nocc     = (neutral.nelec_alpha, neutral.nelec_beta)
    orbs_ion = (ionized.orbs_alpha , ionized.orbs_beta )
    nocc_ion = (ionized.nelec_alpha, ionized.nelec_beta)

    density_matrices = []
    density_matrices_ion = []
    DysonMat = np.zeros((neutral.nstates, 1)).tolist()
    for st in range(neutral.nstates):
        domex = dominant_excitations(neutral.ucis[:, st], excitations,
                                     perc = perc_microstates,
                                     max_microstates = max_microstates)
        P = CIS_density_matrix(orbs, nocc, domex)
        density_matrices.append((nocc[0] + nocc[1], P))

        print("compute Dyson orbital for transition from N-state"+
               " %s to (N-1)-state %s"%(st,0))
        domex_ion = dominant_excitations(ionized.ucis[:, 0], excitations_ion,
                                         perc = perc_microstates,
                                         max_microstates = max_microstates)
        P_ion = CIS_density_matrix(orbs_ion, nocc_ion, domex_ion)
        density_matrices_ion.append((nocc_ion[0] + nocc_ion[1], P_ion))
        print("=====================================")
        print("Dyson orbital between state %s and %s"%(st,0))
        print("-------------------------------------")
        print("neutral state %s is dominated by the excitations"%st)
        neuprec = 0.0
        for (exc, coef) in domex:
            tmp  = "%9.3f %% of "%(coef**2 * 100)
            neuprec += coef**2 * 100

            if exc == []:
                tmp += "[]"
            else:
                tmp += "%2i%s --> %2i%s"%(exc[0][1] + 1, ["A", "B"][exc[0][0]],
                                          exc[0][2] + 1, ["A", "B"][exc[0][0]])
            print(tmp)
        print("  ------------------------")
        print("%9.3f %% total (desired: %7.3f %%)"%(neuprec, perc_microstates))

        print("\nionized state %s is dominated by the excitations"%0)
        ionprec = 0.0
        for (exc, coef) in domex_ion:
            tmp  = "%9.3f %% of "%(coef**2 * 100)
            ionprec += coef**2 * 100

            if exc == []:
                tmp += "[]"
            else:
                tmp += "%2i%s --> %2i%s"%(exc[0][1] + 1, ["A", "B"][exc[0][0]],
                                          exc[0][2] + 1, ["A", "B"][exc[0][0]])
            print(tmp)
        print("  ------------------------")
        print("%9.3f %% total (desired: %7.3f %%)"%(ionprec, perc_microstates))

        print("")
        # energy difference between initial and final state.
        # (energy difference) = (energy of electron) -
        #                       (energy of absorbed photon)
        DysonMat[st][0] = CIS_Dyson(orbs, domex, nocc,
                                    orbs_ion, domex_ion, nocc_ion,
                                    Sao, neutral.basis.nbfsD)

    DysonMat = np.asarray(DysonMat)
    DysonMat = DysonMat.swapaxes(0, 2)

    return DysonMat


def CIS_density_matrix(orbs, nocc, domex):
    nbf = len(orbs[0][:, 0])
    P = [np.zeros((nbf, nbf)), np.zeros((nbf, nbf))]
    # density matrix for spin-up and -down for a CIS state which is a linear
    # combination of singly excited determinants.
    for excs1, coef1 in domex:
        occs1 = occupied_mos(nocc[0], nocc[1], exs=excs1)
        for excs2, coef2 in domex:
            if excs1 == excs2:
                for s in [0, 1]:
                    tmp = np.einsum("mo,no->mn",
                                    orbs[s][:, occs1[s]].conjugate(),
                                    orbs[s][:, occs1[s]])
                    P[s] += coef1.conjugate()*coef2 * tmp

                continue
            (s1, o1, v1) = excs1[0]
            (s2, o2, v2) = excs2[0]
            if o1 == o2 and v1 != v1 and s1 == s2:
                P[s1] = coef1.conjugate() * coef2 \
                    * orbs[s1][:, None, v1].conjugate() * orbs[s1][None, :, v2]
            if o1 != o2 and v1 == v2 and s1 == s2:
                P[s1] = coef1.conjugate() * coef2 \
                    * orbs[s1][:, None, o1].conjugate() * orbs[s1][None, :, o2]

    return P[0] + P[1]


def occupied_mos(nelec_a, nelec_b, exs=[]):
    occs = [list(range(0, nelec_a)), list(range(0, nelec_b))]
    o_last = [-1, -1]
    v_last = [nelec_a-1, nelec_b - 1]
    for (spin, o, v) in exs:
        assert(o > o_last[spin])
        assert(v > v_last[spin])
        occs[spin][o] = v
        o_last[spin] = o
        v_last[spin] = v

    return (occs)


def CIS_Dyson(orbs, domex, nocc, orbs_ion, domex_ion, nocc_ion, Sat, K):
    """
    construct Dyson orbital between neutral and ionized states.
    """

    DysCis_orbs = np.zeros(K)

    for i in range(0, len(domex)):
        (excs, coef) = domex[i]
        occs = occupied_mos(nocc[0], nocc[1], exs=excs)
        for j in range(0, len(domex_ion)):
            (excs_ion, coef_ion) = domex_ion[j]
            occs_ion = occupied_mos(nocc_ion[0], nocc_ion[1], exs=excs_ion)
            d = dyson_orbital(orbs, occs, orbs_ion, occs_ion, Sat, K)
            dDyson = coef.conjugate() * coef_ion * d
            DysCis_orbs += dDyson
            # Dyson orbital is always spin-down

    return (DysCis_orbs)


def dyson_orbital(orbs, occs, orbs_ion, occs_ion, Sat, K):

    nelec_a = len(occs[0])
    nelec_b = len(occs[1])

    # mo coefficient matrices of occupied orbitals only with dimension nelec x K
    # for spin-up and spin-down, Malpha[:,i] are the coefficients for the i-th
    # occupied spin-up molecular orbital in the N-electron state.
    Malpha, Mbeta = orbs[0][:, occs[0]], orbs[1][:, occs[1]]
    Malpha_ion = orbs_ion[0][:, occs_ion[0]]
    Mbeta_ion  = orbs_ion[1][:, occs_ion[1]]
    SmoA = np.dot(np.dot(Malpha_ion.conjugate().T, Sat), Malpha)
    SmoB = np.dot(np.dot(Mbeta_ion.conjugate().T , Sat), Mbeta)

    Dyson_orbs = np.zeros(K)

    if nelec_a == len(occs_ion[0]) + 1:
        # swap alpha and beta
        SmoA, SmoB = SmoB, SmoA
        nelec_a, nelec_b = nelec_b, nelec_a
        Malpha, Mbeta = Mbeta, Malpha

    alpha_overlap = la.det(SmoA)

    for imo in range(0, nelec_b):
        # Minor with dimension K-1 x K-1 obtained from the K x K-1 matrix Smo
        # by deleting the imo-th row, for spin-down.
        minor_beta = np.delete(SmoB, np.s_[imo:imo+1], axis=1)
        if (minor_beta.size == 0):
            Dyson_orbs += alpha_overlap * Mbeta[:, imo]
        else:
            Dyson_orbs += pow(-1,imo) * la.det(minor_beta) * \
                          alpha_overlap * Mbeta[:, imo]

    return (Dyson_orbs)


def order_excs(nocc, nvirt, maxvirts=10e10):
    excs_a = []
    # excitations with spin-up
    for o in range(nocc[0]):
        for v in range(nocc[0], min(maxvirts, nocc[0] + nvirt[0])):
            excs_a.append([(0, o, v)])
    excs_b = []
    # excitations with spin-down
    for o in range(nocc[1]):
        for v in range(nocc[1], min(maxvirts, nocc[1] + nvirt[1])):
            excs_b.append([(1, o, v)])
    excitations = excs_a + excs_b
    return (excitations)


def dominant_excitations(UCISvec, excitations, **opts):
    perc = opts.get("perc", 97.0)
    max_excitations = opts.get("max_excitations", 10)

    indeces = (-abs(UCISvec)).argsort()
    # sort in descending order
    totcontr = 0.0
    domex = []
    for i in range(0, len(indeces)):
        # find dominant states to include to reach desired percentage
        coef = UCISvec[indeces[i]]
        totcontr += pow(coef, 2)
        domex.append((excitations[indeces[i]], coef))
        if (totcontr*100 >= perc):
            break
        if i % 2 == 1:
            # Check every second i to make sure both spin-up and down are added
            if len(domex) > max_excitations:
                break

    return (domex)


def dyson_energies(ionized, neutral):
    """
    energies of Dyson orbitals (ionization energies for N -> I)
    """

    dyson_orbe = np.zeros((ionized.nstates, neutral.nstates), dtype=float)
    for I in range(ionized.nstates):
        for N in range(neutral.nstates):
            dyson_orbe[I,N] = neutral.energies[N] - ionized.energies[I]
    return dyson_orbe


def dyson_norms(dyson_orbs, ionized, neutral, fo):
    """
    compute norms of Dyson orbitals
    """

    # norms of Dyson orbitals
    norms = np.zeros((ionized.nstates, neutral.nstates))
    # overlap between atomic orbitals
    Sao = basis_overlap(neutral.basis, neutral.basis)
    for I in range(ionized.nstates):
        for N in range(neutral.nstates):
            dorb = dyson_orbs[:, I, N]
            norms[I,N] = np.dot(dorb, np.dot(Sao, dorb))
    fo.write("Dyson norms (summed over cation channels)\n")
    for i in np.sum(norms, axis=0):
        fo.write("%.8f "%abs(i))
    fo.write("\n")

    return norms


def dyson_overlap(dyson_orbs1,dyson_orbs2,ionized1,ionized2, neutral1,neutral2):
    """
    compute overlap of Dyson orbitals
    """

    # norms of Dyson orbitals
    norms1   = np.zeros((ionized1.nstates, neutral1.nstates))
    norms2   = np.zeros((ionized2.nstates, neutral2.nstates))
    overlaps = np.zeros((ionized1.nstates, neutral1.nstates))

    # overlap between atomic orbitals
    Sao1 = basis_overlap(ionized1.basis, neutral2.basis)
    Sao2 = basis_overlap(ionized2.basis, neutral1.basis)
    Sao  = basis_overlap(neutral2.basis, neutral1.basis)

    for I in range(ionized1.nstates):
        for N in range(neutral1.nstates):
            dorb1 = dyson_orbs1[:, I, N]
            dorb2 = dyson_orbs2[:, I, N]
            norms1[I, N] = np.sqrt(np.dot(dorb1, np.dot(Sao1, dorb1)))
            norms2[I, N] = np.sqrt(np.dot(dorb2, np.dot(Sao2, dorb2)))
            overlaps[I, N] = np.dot(dorb1, np.dot(Sao, dorb2)) / \
                             (norms1[I,N] * norms2[I,N])

    return overlaps


def dyson_overlap_K(Step1, Step2, nBound, states):
    """
    Dyson overlap between dorbs of different Step instances for all anion states
    """

    dorb1 = Step1.dorbs[:, 0]
    dorb2 = Step2.dorbs[:, 0]
    nrm1 = np.zeros(nBound)
    nrm2 = np.zeros(nBound)
    Sao1 = basis_overlap(Step1.anion.basis, Step1.anion.basis)
    Sao2 = basis_overlap(Step2.anion.basis, Step2.anion.basis)
    SaoD = basis_overlap(Step1.anion.basis, Step2.anion.basis)
    for i, j in enumerate(states):
        nrm1[i] = np.sqrt(la.multi_dot([dorb1[:, j], Sao1, dorb1[:, j]]))
        nrm2[i] = np.sqrt(la.multi_dot([dorb2[:, j], Sao2, dorb2[:, j]]))

    Sdy = np.zeros((nBound, nBound))
    for i, ii in enumerate(states):
        for j, jj in enumerate(states):
            Sdy[i, j] = la.multi_dot([dorb1[:, ii], SaoD, dorb2[:, jj]]) / \
                        (nrm1[i] * nrm2[j])

    return Sdy
