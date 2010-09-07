{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Output.Mol
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

module Ouch.Output.Mol (
     molfile
   ) where


import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Text.String
import Ouch.Property.Property
import Data.Maybe
import Data.Char
import Data.Set as Set
import Data.List as List
import Data.Map as Map
import Control.Applicative


{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

molfile :: Molecule -> Maybe String
molfile m = foldr (\s acc -> (++) <$> s <*> acc) (Just "") lineList
  where m' = removeAtoms m isLonePair
        lineList = List.map (>>= \s -> Just $ s ++ _CR)
                 $ headerBlock m'
                ++ countsLine m'
                ++ atomBlock m'
                ++ bondBlock m'
                ++ propertiesBlock m'
                ++ [Just _END]

headerBlock :: Molecule -> [Maybe String]
headerBlock m = let
  line1 = case getName m of Nothing -> [Just ""]; s -> [s]
  line2 = [Just $ initials ++ _PROGRAM ++ date ++ dim]
  initials = "  "
  date    = "          "
  dim     = "2d"
  scaling = "  "
  energy  = "            "
  registry = "      "
  line3 = case getPropertyForKey m "COMMENT" of
    Nothing -> [Just ""]
    Just s -> [Just $ show $ value s]
  in line1 ++ line2 ++ line3



countsLine :: Molecule -> [Maybe String]
countsLine m = let
  aaa = padCountsElem $ show  $ Map.size $ atomMap m
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
          aaa  = padAtomSymbolElem $ atomicSymbolForAtom a
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

bondBlock :: Molecule -> [Maybe String]
bondBlock m = let
  bondMap = Map.toList $ getBondMap m
  bondLine ((a1, a2), nb) = Just $ atom1 ++ atom2
                                ++ ttt ++ sss ++ xxx
                                ++ rrr ++ ccc
    where atom1 = padBondLineElem $ show (a1 + 1)
          atom2 = padBondLineElem $ show (a2 + 1)
          ttt   = padBondLineElem $ showBondType nb
          sss   = padBondLineElem "0"
          xxx   = padBondLineElem ""
          rrr   = padBondLineElem ""
          ccc   = padBondLineElem "0"
          showBondType b' | b' == Single           = "1"
                          | b' == Double           = "2"
                          | b' == Triple           = "3"
                          | b' == AromaticOnly     = "4"
                          | b' == SingleOrDouble   = "5"
                          | b' == SingleOrAromatic = "6"
                          | b' == DoubleOrAromatic = "7"
                          | b' == AnyBond          = "8"
                          | otherwise              = "8"
  in List.foldr (\b acc -> [bondLine b] ++ acc ) [] bondMap


{------------------------------------------------------------------------------}
{-------------------------------Constants--------------------------------------}
{------------------------------------------------------------------------------}

_CR = "\n"
_VERSION = " V2000"
_PROGRAM = "    OUCH"
_END     = "M  END"
padCountsElem     = padString RightJustify 3 ' '
padAtomLineElem   = padString RightJustify 3 ' '
padAtomSymbolElem = padString LeftJustify 3 ' '
padBondLineElem   = padString RightJustify 3 ' '
padPosElem        = padString RightJustify 10 ' ' . formatNumber 4



