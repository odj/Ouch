{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Atom.hs-boot - a file to define import recursion for the Atom module

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

module Ouch.Structure.Atom (
      Atom(..)
    , atomicSymbolForAtom
    ) where

import Data.Maybe
import {-# SOURCE #-} Ouch.Structure.Bond
import {-# SOURCE #-} Ouch.Structure.Marker
import Data.Set as Set
import Data.Map as Map
import Data.List as List

data Atom   = Element {atomicNumber::Integer
            , neutronNumber::Integer
            , atomBondSet::(Set Bond)
            , atomMarkerSet::(Set AtomMarker)}
            | LonePair {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Electron {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Unfilled {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Unspecified {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Open {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}


atomicSymbolForAtom :: Atom -> String


instance Eq Atom
instance Show Atom
instance Ord Atom
