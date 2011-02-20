{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch
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

{-# LANGUAGE RecordWildCards, CPP #-}
#define BUILD_INFO ("Ouch Version 0.1 * Compiled on: " ++ __DATE__ ++ " | " ++ __TIME__)

module Ouch (
    module Ouch.Structure.Molecule
  , module Ouch.Structure.Bond
  , module Ouch.Structure.Atom
  , module Ouch.Structure.Marker
  , module Ouch.Input.Smiles
  , module Ouch.Output.Smiles
  , module Ouch.Enumerate.Method
  , version
  ) where

import Ouch.Structure.Molecule
import Ouch.Structure.Bond
import Ouch.Structure.Atom
import Ouch.Structure.Marker
import Ouch.Input.Smiles
import Ouch.Output.Smiles hiding (position)
import Ouch.Enumerate.Method


version :: String
version = BUILD_INFO




