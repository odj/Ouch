{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Molecule.hs-boot  - a module to manage molecule data types

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
       Molecule(..),
       removeAtoms,
       getName,
       atomMap,
       getBondMap,
       getPropertyForKey,
       connectMoleculesAtIndicesWithBond,
       bondBetweenIndices,
       freeValenceAtIndex,
       fillMoleculeValence,
       makeMoleculeFromAtom,
       getAtomAtIndex


     ) where

import {-# SOURCE #-} Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import {-# SOURCE #-} Ouch.Property.Builder
import Data.Map as Map


data Molecule

instance Show Molecule
instance Eq Molecule

atomMap :: Molecule -> (Map Int Atom)
removeAtoms :: Molecule -> (Atom -> Bool) -> Molecule
getName :: Molecule -> Maybe String
getBondMap :: Molecule -> Map (Int, Int) NewBond
getPropertyForKey :: Molecule -> String -> Maybe Property
getAtomAtIndex :: Molecule -> Int -> Maybe Atom
connectMoleculesAtIndicesWithBond::Molecule -> Int -> Molecule -> Int -> NewBond -> Molecule
freeValenceAtIndex :: Molecule -> Int -> Integer
fillMoleculeValence :: Molecule -> Molecule
makeMoleculeFromAtom:: Atom -> Molecule
bondBetweenIndices :: Molecule -> Int -> Int -> NewBond




