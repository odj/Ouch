{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Atom.hs-boot - a file to define import recursion
    
    Copyright (c) 2010 Orion D. Jankowski
    
    This file is part of Ouch, a chemical informatics toolkit
    written entirely in Haskell.

    Ouch is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Ouch is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Ouch.  If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------
------------------------------------------------------------------------------}

-- This file is needed to terminate recursive import relationships 

module Ouch.Structure.Atom where

import Ouch.Structure.Bond
import Data.Maybe


import Data.Set as Set
import Data.Map as Map
import Data.List as List
data Atom   = Element {atomicNumber::Integer, neutronNumber::Integer, bondList::[Bond], markerSet::(Set Marker)}
            | LonePair {bondList::[Bond], markerSet::(Set Marker)}
            | Electron {bondList::[Bond], markerSet::(Set Marker)}
            | Unfilled {bondList::[Bond], markerSet::(Set Marker)}
            | Unspecified {bondList::[Bond], markerSet::(Set Marker)}
            | Open {bondList::[Bond], markerSet::(Set Marker)}

            
data Chirality = Levo | Dextro 


data Geometry = Cis {geometetryAtom::Atom} | Trans {geometetryAtom::Atom}


data Marker =  Label {labelNumber::Integer}
              | Chiral {chirality::Chirality}
              | GeoIsomer {geoIsomer::Geometry}
              | Aromatic
              | Traversed
              | Substructure {substructureNumber::Integer}
              | ValenceError {valenceError::String}
              | InRing {ringNumber::Integer}
              | Skip
              | Comment {comment::String}
              | Null  -- This is a dummy value for functions that append marker list for simplicity.


connectAtomsWithBond :: Atom -> Atom -> Bond -> (Atom, Atom)
sigmaBondToAtom :: Atom -> Atom -> (Atom, [Atom])
piBondToAtom :: Atom -> [Atom] -> (Atom, [Atom])
addElectron :: Atom -> [Atom] -> (Atom, [Atom])
addLonePair :: Atom -> [Atom] -> (Atom, [Atom])
addUnfilled :: Atom -> [Atom] -> (Atom, [Atom])
addHydrogen :: Atom -> [Atom] -> (Atom, [Atom]) 
checkValence :: Atom -> Bool
getError :: Atom -> Maybe String
fillValence :: Atom -> [Atom] -> (Atom, [Atom])
atomExactMass :: Atom -> [(Integer, Double)]
atomMW :: Atom -> Double
valence :: Atom -> (Integer, Integer)
isHeavyAtom :: Atom -> Bool
isElement :: Atom -> Bool
numberOfBonds :: Atom -> Integer
numberOfBondsToAtoms :: Atom -> Integer
getBondList :: Atom -> [Bond]
