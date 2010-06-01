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
