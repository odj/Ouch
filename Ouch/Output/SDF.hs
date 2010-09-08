{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Output.SDF
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

module Ouch.Output.SDF (
  sdf
) where

import Ouch.Output.Mol
import Ouch.Structure.Molecule


sdf :: [Molecule] -> String
sdf [] = ""
sdf (m:ms) = case (molfile m) >+> (properties m) >+> (Just _SEP) of
  Nothing -> sdf ms
  Just s  -> s ++ sdf ms


(>+>) :: Maybe String -> Maybe String -> Maybe String
(>+>) s1 s2 = case s1 of
  Nothing -> Nothing
  Just s  -> case s2 of
    Nothing -> Nothing
    Just s' -> Just (s ++ s')


properties :: Molecule -> (Maybe String)
properties m = Just ""


_CR = "\n"
_SEP = "$$$$\n"
