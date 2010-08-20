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


module Ouch.Structure.Molecule(
    Molecule(..)
    ) where

import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Data.Atom
import Ouch.Property.Property

import Data.Either
import Data.Maybe
import Data.Map as Map
import Data.List as List
import Data.Set as Set
import Data.Maybe as Maybe
import Control.Applicative




data Molecule =   Small    {atomMap::(Map Int Atom)
                          , molMarkerSet::(Set MoleculeMarker)
                          , molPropertyMap::(Map Property)}
                | Markush  {molMarkerSet::(Set MoleculeMarker)}
                | Polymer  {molMarkerSet::(Set MoleculeMarker)}
                | Biologic {molMarkerSet::(Set MoleculeMarker)}


instance Show Molecule 
instance Eq Molecule 
