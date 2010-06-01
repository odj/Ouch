-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-- Project !Ouch
-- No license selected yet-- project still under development

-- Orion D. Jankowski
-- 2010-May-24
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

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


