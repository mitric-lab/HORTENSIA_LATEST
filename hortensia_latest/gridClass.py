#!/usr/bin/env python3

import numpy as np
import numpy.linalg as la

import hortensia_latest.misc as misc

################################################################################
#                                                                              #
#                              Plane Wave Grid                                 #
#                                                                              #
################################################################################


class Grid:
    """
    Class containing the discretized plane waves characteristics
    """

    def __init__(self, kDistr, options):
        """
        initializes k grid and attributes thereof
        """

        if kDistr == 'snub' or kDistr == 'fib':
            if kDistr == 'snub':
                nk         = 24
                points, a1 = self.snubCube()
            elif kDistr == 'fib':
                nk         = options[2]
                # Generates nk k-vectors on a sphere in cartesian coordinates
                # using the Fibonacci sphere algorithm
                points, a1 = self.fibSphere(nk)
                # Evaluates the spherical representation of the cartesian
                # vectors
                self.theta, self.phi = self.spherical(points)

            maxEk, nEk  = options[0], options[1]

            constant    = 2*np.pi/3*(1-np.sqrt(1-a1**2))
            self.points = nk * nEk
            self.kbasis = np.zeros((self.points,3))
            self.kfact  = np.zeros(self.points)
            self.kl     = np.zeros(self.points)
            Ediff       = maxEk / float(nEk)

            for i in range(nEk):
                # Calculation of energetic lower (s) and upper (l) bound of
                # volume elements of energetic sphere
                Eks, Ekl = (2*i+1) * Ediff/2, (2*i+3) * Ediff/2
                # Calculation of lower and upper k length of volume element
                rks      = np.sqrt(2*Eks / misc.hartree_to_eV)
                rkl      = np.sqrt(2*Ekl / misc.hartree_to_eV)
                # length of k vectors
                klen     = np.sqrt(2*(i+1) * Ediff/misc.hartree_to_eV)

                # Generation of k vectors through distribution on sphere and
                # k vector lengths
                self.kbasis[i*nk:(i+1)*nk] = points * klen
                # Square of volume element
                self.kfact[i*nk:(i+1)*nk]  = np.sqrt(constant * (rkl**3-rks**3))
                # length of k vectors
                self.kl[i*nk:(i+1)*nk]     = klen

        else:
            # Cubic distribution of k vectors
            # (Seems conceptually intuitive, but results in energetically
            #  non-uniform vector density, therefore not recommended)
            maxEk, nk   = options[0], options[1]
            self.kx     = np.linspace(-maxEk[0], maxEk[0], nk[0])
            self.ky     = np.linspace(-maxEk[1], maxEk[1], nk[1])
            self.kz     = np.linspace(-maxEk[2], maxEk[2], nk[2])
            self.alpha  = 1e-5 * abs((2*min(maxEk))/max(nk))
            self.points = nk[0] * nk[1] * nk[2]
            self.kfact  = np.sqrt(abs(self.kx[0]-self.kx[1]) *
                abs(self.ky[0]-self.ky[1]) * abs(self.kz[0]-self.kz[1]))

            i = 0
            kbasis  = np.zeros((self.points,3))
            self.kl = np.zeros(self.points)
            for x in self.kx:
                for y in self.ky:
                    for z in self.kz:
                        kbasis[i]  = [x,y,z]
                        self.kl[i] = la.norm(kbasis[i])
                        i += 1
            self.kbasis      = kbasis

        ### electron kinetic energy
        self.kinEn = np.zeros(self.points)
        for i,j in enumerate(self.kbasis):
            self.kinEn[i] = la.norm(j)**2/2.


    def snubCube(self):
        """
        calculates vectors of the 24 vertices of a Snub cube with radius 1
        """

        tk     = 1.839286755214161
        points = [
            (-1.0,1/tk,tk),(1.0,-1/tk,tk),(1.0,1/tk,-tk),(-1.0,-1/tk,-tk),
            (1/tk,1.0,tk),(-1/tk,-1.0,tk),(-1/tk,1.0,-tk),(1/tk,-1.0,-tk),
            (tk,1/tk,1.0),(-tk,-1/tk,1.0),(-tk,1/tk,-1.0),(tk,-1/tk,-1.0),
            (1.0,tk,1/tk),(-1.0,-tk,1/tk),(-1.0,tk,-1/tk),(1.0,-tk,-1/tk),
            (-1/tk,tk,1.0),(1/tk,-tk,1.0),(1/tk,tk,-1.0),(-1/tk,-tk,-1.0),
            (-tk,1.0,1/tk),(tk,-1.0,1/tk),(tk,1.0,-1/tk),(-tk,-1.0,-1/tk)]

        points = np.asarray(points) / la.norm(points[1])
        a1     = 1/12 * (5 * la.norm(points[0] - points[4]) +
                             la.norm(points[0] - points[1]))

        return points, a1


    def fibSphere(self, npoints):
        """
        calculates vectors of points on a Fibonacci sphere with radius 1

        https://people.sc.fsu.edu/~jburkardt/py_src/sphere_fibonacci_grid/sphere_fibonacci_grid.html
        """

        phi    = 0.5 * (1 + np.sqrt(5))
        points = np.zeros((npoints, 3))

        for i in range(npoints):
            i2    = float(2*i - (npoints-1))
            theta = 2*np.pi * i2/phi
            sphi  = i2 / npoints
            cphi  = np.sqrt(npoints**2 - i2**2) / npoints

            points[i] = [cphi * np.sin(theta), cphi * np.cos(theta), sphi]

        diff0 = points[0] - points
        norms = la.norm(diff0, axis=1)
        a1    = sum(np.sort(norms)[1:7])/12.

        return points, a1


    def spherical(self, points):
        xy    = points[:,0]**2 + points[:,1]**2
        theta = np.arctan2(np.sqrt(xy), points[:,2])
        phi   = np.arctan2(points[:,1], points[:,0])
        phi[phi<0] += 2*np.pi

        return theta, phi


class GridT:
    def __init__(self, options):
        """
        initializes k grid and attributes thereof
        """

        nk         = 12
        points, a1 = self.icoSphere()
        self.theta, self.phi = self.spherical(points)

        maxEk, nEk  = options[0], options[1]

        constant    = 2*np.pi/3*(1-np.sqrt(1-a1**2))
        self.points = nk * nEk
        self.kbasis = np.zeros((self.points,3))
        self.kfact  = np.zeros(self.points)
        self.kl     = np.zeros(self.points)
        Ediff       = maxEk / float(nEk)

        for i in range(nEk):
            Eks, Ekl = (2*i+1) * Ediff/2, (2*i+3) * Ediff/2
            rks      = np.sqrt(2*Eks / misc.hartree_to_eV)
            rkl      = np.sqrt(2*Ekl / misc.hartree_to_eV)
            klen     = np.sqrt(2*(i+1) * Ediff/misc.hartree_to_eV)

            self.kbasis[i*nk:(i+1)*nk] = points * klen
            self.kfact[i*nk:(i+1)*nk]  = np.sqrt(constant * (rkl**3-rks**3))
            self.kl[i*nk:(i+1)*nk]     = klen


    def icoSphere(self):
        q       = (1+np.sqrt(5))/2
        ilen    = np.sqrt(1+q**2)

        hpoints = np.array([[ 0.0, 1.0, q],[ 0.0, 1.0,-q],
                            [ 1.0, q, 0.0],[ 1.0,-q, 0.0],
                            [ q, 0.0, 1.0],[ q, 0.0,-1.0]]) / ilen

        points  = np.zeros((12,3))
        points[:6] = hpoints
        points[6:] = -hpoints

        a1 = la.norm(points[0]-points[2])/2

        return points, a1


    def spherical(self, points):
        xy    = points[:,0]**2 + points[:,1]**2
        theta = np.arctan2(np.sqrt(xy), points[:,2])
        phi   = np.arctan2(points[:,1], points[:,0])
        phi[phi<0] += 2*np.pi

        return theta, phi
