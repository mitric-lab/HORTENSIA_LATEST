#!/usr/bin/env python

import numpy as np

atomicmasses = {
    "H" :1837.37449 , "B" :19703.62788, "C" :21892.71094, "N" :25532.68181,
    "O" :29163.55664, "F" :34628.81288, "Si":51196.7756 , "P" :56872.47   ,
    "S" :58447.413  , "Cl":64629.86192, "Br":145656.225 , "x" :1.000000000,
    "Au":359047.964 , "Ag":196631.6863, "Li":12653.03019, "Na":41907.77622,
    "Mg":44305.40221, "Sc":81949.5903 , "Ti":87256.362  , "V" :92860.66728,
    "Nb":169357.9289, "Mo":174896.03  ,
    "au":359047.964 , "c" :21892.71094, "o" :29163.55664, "ag":196631.6863,
    "h" :1837.37449 , "n" :25532.68181, "li":12653.03019, "na":41907.77622,
    "mg":44305.40221, "sc":81949.5903 , "ti":87256.362  , "br":145656.225 ,
    "v" :92860.66728, "nb":169357.9289, "si":51196.7756 , "f" :34628.81288,
    "s" :58447.413  , "cl":64629.86192, "mo":174896.03  , "p" :56872.47   }

atomicnumber = {
    "H":  1, "B" :  5, "C" :  6, "N":  7, "O" :  8, "F" :  9, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Br": 35, "h":  1, "b" :  5, "c" :  6, "n" :  7, "o":  8, 
    "f":  9, "si": 14, "p" : 15, "s": 16, "cl": 17, "br": 35}

valence_electrons={"H":1,"C":4,"O":6,"N":5,"Li":1,
		           "h":1,"c":4,"o":6,"n":5,"li":1}

class atom:
    def __init__(self, symbol, coord):
        """
        This is a class defining atoms
        """
        self.symbol = symbol
        self.coord  = coord
        self.basisfunctions = []

    def getMass(self):
        self.mass = atomicmasses[self.symbol]
        return self.mass

    def getAtomicNumber(self):
        self.atomicnumber = atomicnumber[self.symbol]
        return self.atomicnumber

    def getValEl(self):
        return valence_electrons[self.symbol]


if __name__=='__main__':
	testatom = atom("Au", np.array([0.0, 0.0, 0.0]))
