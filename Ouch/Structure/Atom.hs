{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Structure.Atom
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

{-# LANGUAGE CPP #-}
{-# LANGUAGE BangPatterns #-}

module Ouch.Structure.Atom (
      Atom(..)
    , getError
    , getIndexForAtom
    , atomExactMass
    , atomMW
    , valence
    , occupiedValence
    , isHeavyAtom
    , isElement
    , isHydrogen
    , isElectron
    , isLonePair
    , isOpen
    , hasMarker
    , incrementAtom
    , markAtom
    , atomicSymbolForAtom
    , getMatchingClosureNumber
    , removeClosureAtomMarker
    , getMatchingClosureBondType
    , getMarker
    , isElementType

    ) where

import Data.Maybe as Maybe
import Ouch.Data.Atom
import Ouch.Structure.Bond
import {-# SOURCE #-} Ouch.Structure.Molecule
import Ouch.Structure.Marker
import Data.Set as Set
import Data.Map as Map
import Data.List as List

--First Int is atomic number, second Int is number of neutrons (if n=0, then
--assume natural abundance ratio) Bond list is all bond FROM atom.  AtomMarker
--list is used to aid in graph traversing and other functions.
{------------------------------------------------------------------------------}
data Atom   = Element {atomicNumber::Integer
            , neutronNumber::Integer
            , atomBondSet::(Set Bond)
            , atomMarkerSet::(Set AtomMarker)}
            | LonePair {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Electron {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Unfilled {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Unspecified {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}
            | Open {atomBondSet::(Set Bond), atomMarkerSet::(Set AtomMarker)}




getMarkerSet :: Atom -> (Set AtomMarker)
getMarkerSet a = atomMarkerSet a



-- getError
-- Returns a description of the valence problem, if one exists
{------------------------------------------------------------------------------}
getError :: Atom -> Maybe String
getError a = Nothing



-- getExactMass
-- Returns zero for things like electrons and lone-pairs
{------------------------------------------------------------------------------}
atomExactMass :: Atom -> [(Integer, Double)]
atomExactMass a = undefined



-- getMW
-- Returns the isotope averaged molecular weight of the atom
{------------------------------------------------------------------------------}
atomMW :: Atom -> Double
atomMW a = case a of
   Element  {} -> mw $ lookupValue (atomicNumber a)
   LonePair {} -> 0
   Electron {} -> 0
   Unfilled {} -> 0
   Open {}     -> 0
   where lookupValue n = (Map.lookup n atomicWeights)
         mw d = case d of
                Just d -> d
                Nothing -> 0

{------------------------------------------------------------------------------}
atomicSymbolForAtom::Atom -> String
atomicSymbolForAtom a = case a of
    Element {} -> symbol $ lookupValue (atomicNumber a)
    LonePair {} -> ""
    Electron {} -> ""
    Unfilled {} -> ""
    Unspecified {} -> ""
    Open {} -> ""
    where lookupValue n = (Map.lookup n atomicSymbols)
          symbol d = case d of
                 Just d -> d
                 Nothing -> ""

occupiedValence :: Atom -> Integer
occupiedValence a = output
  where output = Set.fold (\b acc -> acc + bondOrder b) 0 (atomBondSet a)
        bondOrder b = case b of
          Sigma    {}   -> 1
          Pi       {}   -> 2
          PiPi     {}   -> 3
          Aromatic {} -> 2
          Delta    {} -> 4
          Hbond    {} -> 0
          Ionic    {} -> 0
          Antibond {} -> 1
          Any      {} -> 0


-- valence
-- Returns the normal valence for the element in question as a tuple.
-- First value is the number of normal bonds.  Second value is the number
-- of lone-pairs.  A negative value for the second number indicates unfilled
-- sp3 (i.e. for Boron).
{------------------------------------------------------------------------------}
valence :: Atom -> (Integer, Integer)
valence !a = let   grp1 = [1,2,11,19,37,55]
                   grp2 = [4,12,20,38,56,88]
                   grp13 = [5,13,31,49,81]
                   grp14 = [6,14,32,50,82]
                   grp15 = [7,15,33,51,83]
                   grp16 = [8,16,34,52,84]
                   grp17 = [9,17,35,53,85]
                   grp18 = [2,10,18,36,54,86]
                   per1 = [1,2]
                   per2 = [3..10]
                   per3 = [11..18]
                   per4 = [19,20,31,32,33,34,35,36]
                   per5 = [37,38,49,50,51,52,53,54]
                   per6 = [55,56,81,82,83,84,85,96]
                   bonds i | elem i grp1 = 1
                           | elem i grp2 = 2
                           | elem i grp13 = 3
                           | elem i grp14 = 4
                           | elem i grp15 = 3
                           | elem i grp16 = 2
                           | elem i grp17 = 1
                           | elem i grp18 = 0
                           | otherwise = 0
                   elecs i | elem i grp1 = 1
                           | elem i grp2 = 2
                           | elem i grp13 = 3
                           | elem i grp14 = 4
                           | elem i grp15 = 5
                           | elem i grp16 = 6
                           | elem i grp17 = 7
                           | elem i grp18 = 8
                           | otherwise = 0
                   sorb i | elem i grp1 = 1
                          | elem i grp2 = 2
                          | elem i grp13 = 4
                          | elem i grp14 = 4
                          | elem i grp15 = 4
                          | elem i grp16 = 4
                          | elem i grp17 = 4
                          | elem i grp18 = 4
                          | otherwise = 0
                   -- This is not right but works for now => CORRECTION NEEDED
                   dorb i | elem i per1 = 0
                          | elem i per2 = 0
                          | elem i per3 = 2
                          | elem i per4 = 2
                          | elem i per5 = 2
                          | elem i per6 = 2
                          | otherwise = 0
                   lp i  = (elecs i) - (sorb i)
            in case a of
               Element  {} -> (bonds (atomicNumber a), lp (atomicNumber a))
               LonePair {} -> (1, 0)
               Electron {} -> (1, 0)
               Unfilled {} -> (1, 0)
               Open     {} -> (1, 0)


{------------------------------------------------------------------------------}
incrementAtom :: Atom -> Int -> Atom
incrementAtom !a !i = a {atomBondSet=newBondSet, atomMarkerSet=newMarkerSet}
  where newBondSet = Set.mapMonotonic incrementBond $ atomBondSet a
        newMarkerSet = Set.insert (Label $ (+i) $ fromJust $ getIndexForAtom a ) $
                                  atomMarkerSet a
        incrementBond b = b {bondsTo=(i + bondsTo b)}

-- markAtom
{------------------------------------------------------------------------------}
getIndexForAtom :: Atom -> Maybe Int
getIndexForAtom a = if n == 0 then Nothing else Just $ fromIntegral (labelNumber $ lb!!0)
  where lb = Set.toList $ Set.filter isLabel (atomMarkerSet a)
        n  = length lb
        isLabel l = case l of Label _ -> True; _ -> False





-- markAtom
{------------------------------------------------------------------------------}
markAtom :: Atom -> AtomMarker -> Atom
markAtom a am = newAtom
   where newSet = Set.insert am (atomMarkerSet a)
         newAtom = a {atomMarkerSet=newSet}



-- isHeavyAtom
-- Returns True for anything that has mass and is not a Hydrogen Atom
{------------------------------------------------------------------------------}
isHeavyAtom :: Atom -> Bool
isHeavyAtom a = case a of
   Element i _ _ _ -> i > 1
   _               -> False


{------------------------------------------------------------------------------}
isElementType :: String -> Atom -> Bool
isElementType s a = case a of
  Element {atomicNumber=an, neutronNumber=_, atomBondSet=_, atomMarkerSet=_} ->
    case num of
    Nothing -> False
    Just n  -> an == n
  _           -> False
  where num = Map.lookup s atomicNumberFromSymbol



-- isElement
-- Returns TRUE if atom is a "real" element
{------------------------------------------------------------------------------}
isElement :: Atom -> Bool
isElement a = case a of
  Element _ _ _ _ -> True
  _               -> False



{------------------------------------------------------------------------------}
isElectron :: Atom -> Bool
isElectron a = case a of
    Electron {} -> True
    _           -> False




{------------------------------------------------------------------------------}
isHydrogen :: Atom -> Bool
isHydrogen a = case a of
    Element 1 _ _ _ -> True
    _               -> False



{------------------------------------------------------------------------------}
isLonePair :: Atom -> Bool
isLonePair a = case a of
    LonePair {} -> True
    _           -> False

{------------------------------------------------------------------------------}
isOpen :: Atom -> Bool
isOpen a = case a of
    Open {} -> True
    _       -> False



hasMarker :: Atom -> AtomMarker -> Bool
hasMarker !a !mk = Set.member mk $ atomMarkerSet a

getMarker :: Atom -> AtomMarker -> Maybe AtomMarker
getMarker !a !mk | hasMarker a mk = Just mk' | otherwise = Nothing
  where mk' = Set.findMax $ Set.filter (\mk' -> EQ == compare mk mk') $ atomMarkerSet a


{------------------------------------------------------------------------------}
removeClosureAtomMarker :: Atom -> Int -> Atom
removeClosureAtomMarker !a !i = a{atomMarkerSet = (Set.delete deleteAtomMarker $ atomMarkerSet a)}
    where deleteAtomMarker = Closure i Single



{------------------------------------------------------------------------------}
getMatchingClosureNumber :: Atom -> Atom -> Maybe Int
getMatchingClosureNumber !a1 !a2 = firstCommonAtomMarker
    where markers1 = fst $ Set.partition isClosure (atomMarkerSet a1)
          markers2 = fst $ Set.partition isClosure (atomMarkerSet a2)
          intersectionSet = Set.intersection markers1 markers2
          isClosure mk  = case mk of Closure {} -> True ; _ -> False
          firstCommonAtomMarker = if (Set.null intersectionSet)
                              then Nothing
                              else Just $ labelNumber (Set.findMin intersectionSet)

{------------------------------------------------------------------------------}
getMatchingClosureBondType :: Atom -> Atom -> NewBond
getMatchingClosureBondType !a1 !a2 = newClosureBond
  where markers1 = fst $ Set.partition isClosure (atomMarkerSet a1)
        markers2 = fst $ Set.partition isClosure (atomMarkerSet a2)
        intersectionSet = Set.intersection markers1 markers2
        isClosure mk  = case mk of Closure {} -> True ; _ -> False
        newClosureBond    = if (Set.null intersectionSet)
                            then NoBond
                            else bondType (Set.findMax intersectionSet)


{------------------------------------------------------------------------------}
{-------------------------------Typeclass Intances-----------------------------}
{------------------------------------------------------------------------------}




{------------------------------------------------------------------------------}
instance Show Atom where
    show a = case a of
          Element i _ b m ->  name i ++ "\t" ++ show b ++ "\t" ++ show m ++ "\n"
          LonePair b m -> "LP\t" ++ show b ++ " " ++ show m ++ "\n"
          Electron b m -> "E\t" ++ show b ++ " " ++ show m ++ "\n"
          Unfilled  {} -> "\t"
          Open     b m -> "OPEN\t" ++ show b ++ " " ++ show m ++ "\n"
          where name b = fromJust $ Map.lookup b atomicSymbols

-- Instance Eq and instance Ord are going to be where all the action is.
-- Everything broken until then!
{------------------------------------------------------------------------------}
instance Eq Atom where
    a == b = (getIndexForAtom a) == (getIndexForAtom b)


{------------------------------------------------------------------------------}
instance Ord Atom where
    compare a b | (getIndexForAtom a) >  (getIndexForAtom b) = GT
                | (getIndexForAtom a) <  (getIndexForAtom b) = LT
                | (getIndexForAtom a) == (getIndexForAtom b) = EQ

