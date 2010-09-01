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
   , atomBlock
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
molfile m = foldr (\s acc -> (++) <$> s <*> acc) (Just "") lineList
  where lineList = List.map (>>= \s -> Just $ s ++ _CR)
                 $ headerBlock m
                ++ countsLine m
                ++ atomBlock m
                ++ propertiesBlock m


countsLine :: Molecule -> [Maybe String]
countsLine m = let
  aaa = padCountsElem $ show  $ Map.size $ Map.filter isElement $ atomMap m
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
atomBlock m = let
  atomLine a = Just $ posX ++ posY ++ posZ ++ " "
                    ++ aaa ++ dd ++ ccc ++ sss ++ hhh
                    ++ bbb ++ vvv ++ hHH ++ rrr ++ iii
                    ++ mmm ++ eee
    where hasPos = hasMarker a $ Position (0, 0, 0)
          hasChg = hasMarker a $ Charge 0
          (x, y, z) | hasPos    = position $ fromJust $ getMarker a $ Position (0, 0, 0)
                    | otherwise = (0 ,0, 0)
          ccc       | hasChg    = padAtomLineElem $ showCharge $ charge $ fromJust $ getMarker a $ Charge 0
                    | otherwise = padAtomLineElem "0"
          posX = padPosElem x
          posY = padPosElem y
          posZ = padPosElem z
          aaa  = padAtomLineElem $ show $ atomicNumber a
          dd   = padString RightJustify 2 ' ' $ show 0
          sss  = padAtomLineElem $ show $ 0
          hhh  = padAtomLineElem $ show $ 0
          bbb  = padAtomLineElem $ show $ 0
          vvv  = padAtomLineElem $ show $ 0
          hHH  = padAtomLineElem $ show $ 0
          rrr  = padAtomLineElem ""
          iii  = padAtomLineElem ""
          mmm  = padAtomLineElem ""
          eee  = padAtomLineElem ""
          showCharge i | i ==  3 = "1"
                       | i ==  2 = "2"
                       | i ==  1 = "3"
                       | i == -1 = "5"
                       | i == -2 = "6"
                       | i == -3 = "7"
                       | otherwise = "0"
  in Map.fold (\a acc -> [atomLine a] ++ acc ) [] $ Map.filter isElement $ atomMap m

propertiesBlock :: Molecule -> [Maybe String]
propertiesBlock m = [Just ""]

headerBlock :: Molecule -> [Maybe String]
headerBlock m = [Just ""]


{------------------------------------------------------------------------------}
{-------------------------------Constants--------------------------------------}
{------------------------------------------------------------------------------}

_CR = "\n"
_VERSION = " V2000"
padCountsElem = padString RightJustify 3 ' '
padAtomLineElem = padString RightJustify 3 ' '
padPosElem = padString RightJustify 10 ' ' . formatNumber 4



