{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Atoms - a module to manage atom data
    
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

{-# LANGUAGE CPP #-}
module Ouch.Structure.Atom (
      Atom(..)
    , addHydrogen
    , addLonePair
    , addUnfilled
    , bondsToHeavyAtomsAtIndices
    , checkValence
    , getError
    , fillValence
    , getIndexForAtom
    , atomExactMass
    , atomMW
    , valence
    , isHeavyAtom
    , isElement
    , isSigmaBondToHeavyAtom
    , markAtom
    , numberOfBonds
    , numberOfBondsToAtoms
    , connectAtomsWithBond
    , atomicSymbolForAtom
    , getMatchingClosureNumber
    , removeClosureAtomMarker
    , getMatchingClosureBondType
    , getMarkerSet -- This should be depreciated
    ) where

import Data.Maybe
import Ouch.Data.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Data.Set as Set
import qualified Data.Map as Map
import qualified Data.List as List

-- First Int is atomic number, second Int is number of neutrons (if n=0, then assume natural abundance ratio)
-- Bond list is all bond FROM atom.  
-- AtomMarker list is used to aid in graph traversing and other functions.
{------------------------------------------------------------------------------}
data Atom   = Element {atomicNumber::Integer, neutronNumber::Integer, bondList::[Bond], markerSet::(Set AtomMarker)}
            | LonePair {bondList::[Bond], markerSet::(Set AtomMarker)}
            | Electron {bondList::[Bond], markerSet::(Set AtomMarker)}
            | Unfilled {bondList::[Bond], markerSet::(Set AtomMarker)}
            | Unspecified {bondList::[Bond], markerSet::(Set AtomMarker)}   --Wildcard atom for smiles symbol *
            | Open {bondList::[Bond], markerSet::(Set AtomMarker)}


-- connectAtomsWithBond
{------------------------------------------------------------------------------}
connectAtomsWithBond :: Atom -> Atom -> NewBond -> (Atom, Atom)
connectAtomsWithBond a1 a2 b = (aa1, aa2) where
      aa1 = case a1 of
          Element  {} ->  Element     { atomicNumber=(atomicNumber a1)
                                      , neutronNumber=(neutronNumber a1)
                                      , bondList=((bondList a1) ++ newBond)
                                      , markerSet=(markerSet a1)}
          LonePair {} -> LonePair     { bondList=((bondList a1) ++ newBond)
                                      , markerSet=(markerSet a1)}
          Electron {} -> Electron     { bondList=((bondList a1) ++ newBond)
                                      , markerSet=(markerSet a1)}
          Unfilled {} -> Unfilled     { bondList=((bondList a1) ++ newBond)
                                      , markerSet=(markerSet a1)}
          where newBond = case b of
                              Single -> [Sigma aa2]
                              Double -> [Sigma aa2] ++ [Pi aa2]
                              Triple -> [Sigma aa2] ++ [Pi aa2] ++ [Pi aa2]
                              NoBond -> []
      aa2 = case a2 of
          Element  {} -> Element      { atomicNumber=(atomicNumber a2)
                                      , neutronNumber=(neutronNumber a2)
                                      , bondList=((bondList a2) ++ newBond)
                                      , markerSet=(markerSet a2)}
          LonePair {} -> LonePair     { bondList=((bondList a2) ++ newBond)
                                      , markerSet=(markerSet a2)}
          Electron {} -> Electron     { bondList=((bondList a2) ++ newBond)
                                      , markerSet=(markerSet a2)}
          Unfilled {} -> Unfilled     { bondList=((bondList a2) ++ newBond)
                                      , markerSet=(markerSet a2)}
          where newBond = case b of
                            Single -> [Sigma aa1]
                            Double -> [Sigma aa1] ++ [Pi aa1]
                            Triple -> [Sigma aa1] ++ [Pi aa1] ++ [Pi aa1]
                            NoBond -> []


getMarkerSet :: Atom -> (Set AtomMarker)
getMarkerSet a = markerSet a

-- addLonePair
-- Create a new lone-pair centered on the atom.
-- Return atom and list containing new lone-pair.
{------------------------------------------------------------------------------}
addLonePair :: Atom -> [Atom] -> (Atom, [Atom])
addLonePair a as  = (a', ([as'] ++ as))
   where (a', as') = connectAtomsWithBond a (LonePair [] Set.empty) Single
         val  = (fst $ valence a) + (abs(snd $ valence a))
         nb   = numberOfBonds a
         




{------------------------------------------------------------------------------}
addHydrogen :: Atom -> [Atom] -> (Atom, [Atom]) 
addHydrogen a as = (a', ([as'] ++ as))
    where (a', as') = connectAtomsWithBond a (Element 1 1 [] Set.empty) Single
          val  = fst $ valence a
          nb   = numberOfBondsToAtoms a
    

addElectron :: Atom -> [Atom] -> (Atom, [Atom]) 
addElectron a as  = (a', ([as'] ++ as))
  where (a', as') = connectAtomsWithBond a (Electron [] Set.empty) Single
        val  = fst $ valence a
        nb   = numberOfBondsToAtoms a



-- addUnfilled
-- Create a new unfilled orbital centered on the atom.
-- Return atom and list containing new unfilled orbital.
{------------------------------------------------------------------------------}
addUnfilled :: Atom -> [Atom] -> (Atom, [Atom])
addUnfilled a = undefined
  
  
   
-- checkValence
-- Verify valence rules are met.  True is what you want.
{------------------------------------------------------------------------------}
checkValence :: Atom -> Bool
checkValence a = True



-- getError
-- Returns a description of the valence problem, if one exists
{------------------------------------------------------------------------------}
getError :: Atom -> Maybe String
getError a = Nothing



-- fillValence
-- Populate free valences with hydrogens/lone-pairs.  Return new atom plus an
-- atom list containing all the hydrogens added (adding to second arg).  
-- Lone pairs (i.e. for Nitrgen atoms) and empty orbitals (i.e. on Boron)
-- are also added and incuded in the list.
{------------------------------------------------------------------------------}
fillValence :: Atom -> [Atom] -> (Atom, [Atom])
fillValence a as =  
    let val          = valence a
        hBool        = Set.member (ExplicitHydrogen 0) (markerSet a)
        h            | hBool = numberH $ Set.findMax $ Set.filter (== (ExplicitHydrogen 0)) (markerSet a)
                     | otherwise = 0
        nba          = (numberOfBondsToAtoms a) + (numberOfBondsToRadicals a)
        nb           = numberOfBonds a
        nbh          = numberOfBondsToHydrogens a
        (aH, asH)    = addHydrogen a as
        (aLP, asLP)  = addLonePair a as
        (aEL, asEL)  = addElectron a as
        nbrB         = (numberOfBondsToRadicals a) == 0
        outputXH   | nbh < h                                       = fillValence aH asH
                   | nb >= ((fst val) + (abs(snd val)))            = (a, as) 
                     -- Fill marked aromatics with a radical
                   | nbrB && Set.member AromaticAtom (markerSet a) = fillValence aEL asEL

                     -- Fill empty valences with radical
                   | nba < fst val                                 = fillValence aEL asEL

                     -- Then, fill lone-pairs if needed
                   | nb < ((fst val) + (abs(snd val)))             = fillValence aLP asLP

                     -- Pattern completion
                   | otherwise = (a, as)
        output     | nb >= ((fst val) + (abs(snd val)))            = (a, as) 
                     -- Fill marked aromatics with a radical
                   | nbrB && Set.member AromaticAtom (markerSet a) = fillValence aEL asEL

                     -- Fill empty valences with hydrogen
                   | nba < fst val                                 = fillValence aH asH

                     -- Then, fill lone-pairs if needed
                   | nb < ((fst val) + (abs(snd val)))             = fillValence aLP asLP

                     -- Pattern completion
                   | otherwise = (a, as)
    in if hBool then outputXH else output
  
   
       
       



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
   
   
-- valence
-- Returns the normal valence for the element in question as a tuple.
-- First value is the number of normal bonds.  Second value is the number 
-- of lone-pairs.  A negative value for the second number indicates unfilled
-- sp3 (i.e. for Boron).
{------------------------------------------------------------------------------}
valence :: Atom -> (Integer, Integer)
valence a = let    grp1 = [1,2,11,19,37,55]
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
                    
                       

-- markAtom
{------------------------------------------------------------------------------}
getIndexForAtom :: Atom -> Maybe Int
getIndexForAtom a = if n == 0 then Nothing else Just $ fromIntegral (labelNumber $ lb!!0)
  where lb = Set.toList $ Set.filter (==(Label 0)) (markerSet a)
        n  = length lb





-- markAtom
{------------------------------------------------------------------------------}
markAtom :: Atom -> AtomMarker -> Atom
markAtom a am = newAtom
   where newSet = Set.insert am (markerSet a)
         newAtom = a {markerSet=newSet}


-- isHeavyAtom
-- Returns True for anything that has mass and is not a Hydrogen Atom
{------------------------------------------------------------------------------}
isHeavyAtom :: Atom -> Bool
isHeavyAtom a = case a of
   Element i _ _ _ -> i > 1
   LonePair {} -> False
   Electron {} -> False
   Unfilled {} -> False
  
  
   
-- isElement
-- Returns TRUE if atom is a "real" element
{------------------------------------------------------------------------------}
isElement :: Atom -> Bool
isElement a = case a of
  Element _ _ _ _ -> True
  LonePair {} -> False
  Electron {} -> False
  Unfilled {} -> False

-- isElementOrRadical
-- Returns TRUE if atom is a "real" element
{------------------------------------------------------------------------------}
isElectron :: Atom -> Bool
isElectron a = case a of
    Element _ _ _ _ -> False
    LonePair {} -> False
    Electron {} -> True
    Unfilled {} -> False


-- numberOfBonds
-- Returns number of covalent connections to other atoms in the molecule 
-- graph (i.e. one sigma and two pi bonds count as a 'one' bond)
{------------------------------------------------------------------------------}
numberOfBonds :: Atom -> Integer
numberOfBonds a = fromIntegral $ length $ bondList a




{------------------------------------------------------------------------------}
numberOfBondsToAtoms :: Atom -> Integer
numberOfBondsToAtoms a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = fromIntegral $ length $ List.filter (==True) $ List.map isAnyBondToAtom b 

{------------------------------------------------------------------------------}
numberOfBondsToRadicals :: Atom -> Integer
numberOfBondsToRadicals a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = fromIntegral $ length $ List.filter (==True) $ List.map isAnyBondToRadical b


{------------------------------------------------------------------------------}
numberOfBondsToHydrogens :: Atom -> Integer
numberOfBondsToHydrogens a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = fromIntegral $ length $ List.filter (==True) $ List.map (\a -> isAnyBondToElement a 1) b 



{------------------------------------------------------------------------------}
numberOfAromaticBondsToAtoms :: Atom -> Integer
numberOfAromaticBondsToAtoms a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = fromIntegral $ length $ List.filter (==True) $ List.map isAnyBondToRadical b    

{------------------------------------------------------------------------------}
bondsToHeavyAtomsAtIndices :: Atom -> [Int]
bondsToHeavyAtomsAtIndices atom = mapMaybe (getIndexForAtom . bondsTo) $ List.filter isSigmaBondToHeavyAtom $ bondList atom

-- currentValence
-- Aromatic bonds count as ONE bond, no matter how many there are.  Aromatic atoms that are
-- incorrectly indicated but not part of a ring system are effectively "radical" system, but
-- will not be explicitely depicted as such in the data structure.  Aromatic atoms that are
-- not directly connected to any other aromatic atom will be considered ?? what exactly???
{------------------------------------------------------------------------------}
currentValence :: Atom ->  Integer
currentValence a = undefined



{------------------------------------------------------------------------------}
isAnyBondToAtom :: Bond -> Bool
isAnyBondToAtom b = case b of
    Sigma atom -> isElement atom
    Pi atom -> isElement atom
    Aromatic atom -> isElement atom
    Delta atom -> isElement atom
    Hbond atom -> False
    Ionic atom -> False
    Antibond atom -> False
    Any atom -> True
    
{------------------------------------------------------------------------------}
isAnyBondToRadical :: Bond -> Bool
isAnyBondToRadical b = case b of
    Sigma atom -> isElectron atom
    Pi _ -> False
    Aromatic _ -> False
    Delta _ -> False
    Hbond _ -> False
    Ionic _ -> False
    Antibond _ -> False
    Any atom -> isElectron atom


{------------------------------------------------------------------------------}
isSigmaBondToHeavyAtom :: Bond -> Bool
isSigmaBondToHeavyAtom b =  
    case b of   Sigma atom -> (isHeavyAtom $ bondsTo b)
                _          -> False    



{------------------------------------------------------------------------------}
isAnyBondToElement :: Bond -> Integer -> Bool
isAnyBondToElement b i = let an | isElement (bondsTo b) = atomicNumber (bondsTo b) | otherwise = 0 in 
    case b of   Sigma atom -> an == i
                Pi atom -> an == i
                Aromatic atom -> an == i
                Delta atom -> an == i
                Hbond atom -> an == i
                Ionic atom -> an == i
                Antibond atom -> an == i
                Any atom -> an == i


{------------------------------------------------------------------------------}
isDeltaBondToAtom :: Bond -> Bool
isDeltaBondToAtom b = case b of
    Sigma atom -> False
    Pi atom -> False
    Aromatic atom -> False
    Delta atom -> isElement atom
    Hbond atom -> False
    Ionic atom -> False
    Antibond atom -> False
    Any atom -> True



{------------------------------------------------------------------------------}
isAromaticBondToAtom :: Bond -> Bool
isAromaticBondToAtom b = case b of
    Sigma atom -> False
    Pi atom -> False
    Aromatic atom -> isElement atom
    Delta atom -> False
    Hbond atom -> False
    Ionic atom -> False
    Antibond atom -> False
    Any atom -> True



{------------------------------------------------------------------------------}
isPiBondToAtom :: Bond -> Bool
isPiBondToAtom b = case b of
    Sigma atom -> False
    Pi atom -> isElement atom
    Aromatic atom -> False
    Delta atom -> False
    Hbond atom -> False
    Ionic atom -> False
    Antibond atom -> False
    Any atom -> True


{------------------------------------------------------------------------------}
isSigmaBondToAtom :: Bond -> Bool
isSigmaBondToAtom b = case b of
    Sigma atom -> isElement atom
    Pi atom -> False
    Aromatic atom -> False
    Delta atom -> False
    Hbond atom -> False
    Ionic atom -> False
    Antibond atom -> False
    Any atom -> True


{------------------------------------------------------------------------------}
removeClosureAtomMarker :: Atom -> Integer -> Atom
removeClosureAtomMarker a i = a{markerSet = (Set.delete deleteAtomMarker $ markerSet a)} 
    where deleteAtomMarker = Closure i Single



{------------------------------------------------------------------------------}
getMatchingClosureNumber :: Atom -> Atom -> Maybe Integer
getMatchingClosureNumber a1 a2 = firstCommonAtomMarker
    where markers1 = fst $ Set.partition isClosure (markerSet a1)
          markers2 = fst $ Set.partition isClosure (markerSet a2)
          intersectionSet = Set.intersection markers1 markers2
          isClosure mk  = case mk of Closure {} -> True ; _ -> False
          firstCommonAtomMarker = if (Set.null intersectionSet) 
                              then Nothing 
                              else Just $ labelNumber (Set.findMin intersectionSet)

{------------------------------------------------------------------------------}
getMatchingClosureBondType :: Atom -> Atom -> NewBond
getMatchingClosureBondType a1 a2 = newClosureBond
  where markers1 = fst $ Set.partition isClosure (markerSet a1)
        markers2 = fst $ Set.partition isClosure (markerSet a2)
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
          LonePair b m -> "LP" ++ show b ++ " " ++ show m ++ "\n"
          Electron b m -> "E" ++ show b ++ " " ++ show m ++ "\n"
          Unfilled {} -> ""
          where name b = fromJust $ Map.lookup b atomicSymbols

-- Instance Eq and instance Ord are going to be where all the action is.
-- Everything broken until then!
{------------------------------------------------------------------------------}
instance Eq Atom where
    a == b = True

    
{------------------------------------------------------------------------------}
instance Ord Atom where
    a > b = True
    a < b = False
