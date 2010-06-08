{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Molecule.hs-boot - a file to define import recursion
    
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

module Ouch.Structure.Molecule 
    (
       Molecule(..)
     , PartialMolecule(..)
     , addAtom
     , addMolecule
     , numberOfAtoms
     , numberOfHeavyAtoms
     , fillMoleculeValence
     , molecularWeight
     , exactMass
     , numberOfHydrogenBondDonors
     , numberOfHydrogenBondAcceptors
     , numberOfRings
     , numberOfRotatableBonds
     , molecularFormula
     , connectMoleculesAtIndicesWithBond
     ) where

import {-# SOURCE #-} Ouch.Structure.Atom
import Data.Either

data Molecule 
type PartialMolecule 
addAtom :: Molecule -> Atom -> Molecule

connectMoleculesAtIndicesWithBond :: Molecule -> Integer -> Molecule -> Integer -> Bond -> PartialMolecule

addMolecule :: Molecule -> Molecule -> Molecule

makeMoleculeFromAtom:: Atom -> Molecule

numberOfAtoms :: Molecule -> Maybe Integer

numberOfHeavyAtoms :: Molecule -> Maybe Integer

fillMoleculeValence :: Molecule -> Either Molecule Molecule

molecularWeight :: Molecule -> Maybe Double

exactMass :: Molecule -> Maybe [(Integer, Double)]

numberOfHydrogenBondDonors :: Molecule -> Maybe Integer


numberOfHydrogenBondAcceptors :: Molecule -> Maybe Integer


numberOfRings :: Molecule -> Maybe Integer

numberOfRotatableBonds :: Molecule -> Maybe Integer


molecularFormula :: Molecule -> Maybe String


