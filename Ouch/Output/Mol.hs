{------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Mol - a module to define export into MOL file format

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
-------------------------------------------------------------------------------}

module Ouch.Output.Mol (
     molfile
     ---

   , countsLine
   ) where


import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Text.String
import Data.Maybe
import Data.Char
import Data.Set as Set
import Data.List as List
import Data.Map as Map
import Control.Applicative



{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}


{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

molfile :: Molecule -> Maybe String
molfile m = foldr (\acc s -> (++) <$> acc <*> s) (Just "") lineList
  where lineList = headerBlock m
                ++ countsLine m
                ++ atomBlock m
                ++ propertiesBlock m


countsLine :: Molecule -> [Maybe String]
countsLine m = let
  aaa = padCountsElem $ show $ Map.size $ atomMap m
  bbb = padCountsElem $ show $ Map.size $ getBondMap m
  lll = padCountsElem "0"
  fff = padCountsElem "0"
  ccc | isChiral = padCountsElem "1" | otherwise = padCountsElem "0"
  sss = padCountsElem "0"
  xxx = padCountsElem "0"
  mmm = padCountsElem "999"
  isChiral = Map.fold (\a acc -> (&&) acc $ hasMarker a $ Chiral Dextro) False (atomMap m)
  in [Just $ aaa
          ++ bbb
          ++ lll
          ++ fff
          ++ ccc
          ++ sss
          ++ xxx ++ xxx ++ xxx ++ xxx
          ++ mmm
          ++ _VERSION]

atomBlock :: Molecule -> [Maybe String]
atomBlock m = undefined

propertiesBlock :: Molecule -> [Maybe String]
propertiesBlock m = undefined

headerBlock :: Molecule -> [Maybe String]
headerBlock m = undefined


{------------------------------------------------------------------------------}
{-------------------------------Constants--------------------------------------}
{------------------------------------------------------------------------------}

_CR = "\n"
_VERSION = " V2000"
padCountsElem = padString RightJustify 3 ' '



