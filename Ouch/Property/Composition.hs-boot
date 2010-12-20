
{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Property.Composition
--  Maintainer  :  Orion Jankowski
--  Stability   :  Unstable
--  Portability :


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

--------------------------------------------------------------------------------
-------------------------------------------------------------------------------}

module Ouch.Property.Composition
  (
   molecularWeight
 , molecularFormula
 , exactMass
 , atomCount
 , heavyAtomCount
 , hBondAcceptorCount
 , hBondDonorCount
 , netCharge
  ) where

import {-# SOURCE #-} Ouch.Structure.Molecule
import {-# SOURCE #-} Ouch.Property.Builder


exactMass :: Property
heavyAtomCount :: Property
atomCount :: Property
hBondAcceptorCount :: Molecule -> Maybe Property
hBondDonorCount :: Molecule -> Maybe Property
netCharge :: Molecule -> Maybe Property
molecularWeight :: Property
molecularFormula :: Property

