-- Project !Ouch
-- No license selected yet-- project still under development

-- Orion D. Jankowski
-- 2010-May-24

{-# LANGUAGE RecordWildCards, CPP #-}
module Ouch.Structure.Atom (
      Atom(..)
    , Chirality(..)
    , Marker(..)
    , Bond(..)
    , NewBond(..)
    , sigmaBondToAtom
    , addHydrogen
    , piBondToAtom
    , addElectron
    , addLonePair
    , addUnfilled
    , checkValence
    , getError
    , fillValence
    , atomExactMass
    , atomMW
    , valence
    , isHeavyAtom
    , isElement
    , numberOfBonds
    , numberOfBondsToAtoms
    , connectAtomsWithBond
    , getBondList
    , atomicSymbolForAtom
    , getMatchingClosureNumber
    , removeClosureMarker
    ) where

import Data.Maybe
import Ouch.Data.Atom


import Data.Set as Set
import Data.Map as Map
import Data.List as List

-- First Int is atomic number, second Int is number of neutrons (if n=0, then assume natural abundance ratio)
-- Bond list is all bond FROM atom.  
-- Marker list is used to aid in graph traversing and other functions.
data Atom   = Element {atomicNumber::Integer, neutronNumber::Integer, bondList::[Bond], markerSet::(Set Marker)}
            | LonePair {bondList::[Bond], markerSet::(Set Marker)}
            | Electron {bondList::[Bond], markerSet::(Set Marker)}
            | Unfilled {bondList::[Bond], markerSet::(Set Marker)}
            | Unspecified {bondList::[Bond], markerSet::(Set Marker)}   --Wildcard atom for smiles symbol *
            | Open {bondList::[Bond], markerSet::(Set Marker)}
            deriving (Eq)
         
data Chirality = Levo | Dextro 
     deriving (Show, Eq, Ord)

data Geometry = Cis {geometetryAtom::Atom} | Trans {geometetryAtom::Atom}
     deriving (Show, Eq, Ord)

data Marker =  Label {labelNumber::Integer}   -- OUCH specific label
              | Closure {labelNumber::Integer, bondType::NewBond}
              | Class {classNumber::Integer}
              | Chiral {chirality::Chirality}
              | GeoIsomer {geoIsomer::Geometry}
              | AromaticAtom
              | Traversed
              | Substructure {substructureNumber::Integer}
              | ValenceError {valenceError::String}
              | InRing {ringNumber::Integer}
              | Skip
              | Comment {comment::String}
              | Null  -- This is a dummy value for functions that append marker list for simplicity.
              deriving (Show, Ord)

data Bond = Sigma {bondsTo::Atom}
          | Pi {bondsTo::Atom} 
          | Aromatic {bondsTo::Atom}
          | Delta {bondsTo::Atom}
          | Hbond {bondsTo::Atom}
          | Ionic {bondsTo::Atom}
          | Antibond {bondsTo::Atom}
          | Any {bondsTo::Atom}
          deriving (Eq)

data NewBond = Single | Double | Triple | NoBond deriving (Show, Eq, Ord)


connectAtomsWithBond :: Atom -> Atom -> NewBond -> (Atom, Atom)
connectAtomsWithBond a1 a2 b = (aa1, aa2)
  where
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



-- sigmaBondToAtom
-- Create a new sigma bond between two atoms.
-- Return atom and list containing second atom.  List is empty if
-- function cannot add bond.   
-- What chemcial error checking should be done here??
sigmaBondToAtom :: Atom -> Atom -> (Atom, [Atom])
sigmaBondToAtom a1 a2 = (aa1, [aa2])
    where
        aa1 = case a1 of
            Element  {} ->  Element     { atomicNumber=(atomicNumber a1)
                                        , neutronNumber=(neutronNumber a1)
                                        , bondList=((bondList a1) ++ [Sigma aa2])
                                        , markerSet=(markerSet a1)}
            LonePair {} -> LonePair     { bondList=((bondList a1) ++ [Sigma aa2])
                                        , markerSet=(markerSet a1)}
            Electron {} -> Electron     { bondList=((bondList a1) ++ [Sigma aa2])
                                        , markerSet=(markerSet a1)}
            Unfilled {} -> Unfilled     { bondList=((bondList a1) ++ [Sigma aa2])
                                        , markerSet=(markerSet a1)}
        aa2 = case a2 of
            Element  {} -> Element      { atomicNumber=(atomicNumber a2)
                                        , neutronNumber=(neutronNumber a2)
                                        , bondList=((bondList a2) ++ [Sigma a1])
                                        , markerSet=(markerSet a2)}
            LonePair {} -> LonePair     { bondList=((bondList a2) ++ [Sigma a1])
                                        , markerSet=(markerSet a2)}
            Electron {} -> Electron     { bondList=((bondList a2) ++ [Sigma a1])
                                        , markerSet=(markerSet a2)}
            Unfilled {} -> Unfilled     { bondList=((bondList a2) ++ [Sigma a1])
                                        , markerSet=(markerSet a2)}

-- piBondToAtom
-- Create a new pi bond between two atoms.
-- Return atom and list containing second atom.  List is empty if
-- function cannot add bond.
piBondToAtom :: Atom -> [Atom] -> (Atom, [Atom])
piBondToAtom a1 a2 = undefined

-- addElectron
-- Create a new radical centered on the atom.
-- Return atom and list containing new radical.
addElectron :: Atom -> [Atom] -> (Atom, [Atom])
addElectron a = undefined

-- addLonePair
-- Create a new lone-pair centered on the atom.
-- Return atom and list containing new lone-pair.
addLonePair :: Atom -> [Atom] -> (Atom, [Atom])
addLonePair a as 
   | (nb < val)    = (a', (as' ++ as))
   | otherwise     = (a, as)
   where (a', as') = sigmaBondToAtom a (LonePair [] Set.empty)
         val  = (fst $ valence a) + (abs(snd $ valence a))
         nb   = numberOfBonds a
         




addHydrogen :: Atom -> [Atom] -> (Atom, [Atom]) 
addHydrogen a as 
    | (nb < val)    = (a', (as' ++ as))
    | otherwise     = (a, as)
    where (a', as') = sigmaBondToAtom a (Element 1 1 [] Set.empty)
          val  = fst $ valence a
          nb   = numberOfBondsToAtoms a
    


-- addUnfilled
-- Create a new unfilled orbital centered on the atom.
-- Return atom and list containing new unfilled orbital.
addUnfilled :: Atom -> [Atom] -> (Atom, [Atom])
addUnfilled a = undefined
   
-- checkValence
-- Verify valence rules are met.  True is what you want.
checkValence :: Atom -> Bool
checkValence a = True

-- getError
-- Returns a description of the valence problem, if one exists
getError :: Atom -> Maybe String
getError a = Nothing

-- fillValence
-- Populate free valences with hydrogens/lone-pairs.  Return new atom plus an
-- atom list containing all the hydrogens added (adding to second arg).  
-- Lone pairs (i.e. for Nitrgen atoms) and empty orbitals (i.e. on Boron)
-- are also added and incuded in the list.
fillValence :: Atom -> [Atom] -> (Atom, [Atom])
fillValence a as  
       -- Everything filled, nothing to do.  Return original
       | nb >= ((fst val) + (abs(snd val)))            = (a, as) 
       
       -- Fill empty valences with hydrogen
       | nba < fst val                                 = fillValence aH asH
       
       -- Then, fill lone-pairs if needed
       | nb < ((fst val) + (abs(snd val)))             = fillValence aLP asLP
       
       -- Pattern completion
       | otherwise = (a, as)
       where
           val          = valence a
           nba          = numberOfBondsToAtoms a
           nb           = numberOfBonds a
           (aH, asH)    = addHydrogen a as
           (aLP, asLP)  = addLonePair a as
    

-- getExactMass
-- Returns zero for things like electrons and lone-pairs
atomExactMass :: Atom -> [(Integer, Double)]
atomExactMass a = undefined

-- getMW
-- Returns the isotope averaged molecular weight of the atom
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
valence :: Atom -> (Integer, Integer)
valence a = case a of
   Element  {} -> (bonds (atomicNumber a), lp (atomicNumber a))
   LonePair {} -> (1, 0)
   Electron {} -> (1, 0)
   Unfilled {} -> (1, 0)
   where 
       bonds i | elem i grp1 = 1
               | elem i grp2 = 2
               | elem i grp13 = 3
               | elem i grp14 = 4
               | elem i grp15 = 3
               | elem i grp16 = 2
               | elem i grp17 = 1
               | elem i grp18 = 0
               | otherwise = 0 
               where
                   grp1 = [1,2,11,19,37,55]
                   grp2 = [4,12,20,38,56,88]
                   grp13 = [5,13,31,49,81]
                   grp14 = [6,14,32,50,82]
                   grp15 = [7,15,33,51,83]
                   grp16 = [8,16,34,52,84]
                   grp17 = [9,17,35,53,85]
                   grp18 = [2,10,18,36,54,86]
       
       elecs i | elem i grp1 = 1
               | elem i grp2 = 2
               | elem i grp13 = 3
               | elem i grp14 = 4
               | elem i grp15 = 5
               | elem i grp16 = 6
               | elem i grp17 = 7
               | elem i grp18 = 8
               | otherwise = 0 
               where
                 grp1 = [1,2,11,19,37,55]
                 grp2 = [4,12,20,38,56,88]
                 grp13 = [5,13,31,49,81]
                 grp14 = [6,14,32,50,82]
                 grp15 = [7,15,33,51,83]
                 grp16 = [8,16,34,52,84]
                 grp17 = [9,17,35,53,85]
                 grp18 = [2,10,18,36,54,86]

       sorb i | elem i grp1 = 1
              | elem i grp2 = 2
              | elem i grp13 = 4
              | elem i grp14 = 4
              | elem i grp15 = 4
              | elem i grp16 = 4
              | elem i grp17 = 4
              | elem i grp18 = 4
              | otherwise = 0 
              where
                  grp1 = [1,2,11,19,37,55]
                  grp2 = [4,12,20,38,56,88]
                  grp13 = [5,13,31,49,81]
                  grp14 = [6,14,32,50,82]
                  grp15 = [7,15,33,51,83]
                  grp16 = [8,16,34,52,84]
                  grp17 = [9,17,35,53,85]
                  grp18 = [2,10,18,36,54,86]
       -- This is not right but works for now => CORRECTION NEEDED
       dorb i | elem i per1 = 0
              | elem i per2 = 0
              | elem i per3 = 2
              | elem i per4 = 2
              | elem i per5 = 2
              | elem i per6 = 2
              | otherwise = 0
              where
                  per1 = [1,2]
                  per2 = [3..10]
                  per3 = [11..18]
                  per4 = [19,20,31,32,33,34,35,36]
                  per5 = [37,38,49,50,51,52,53,54]
                  per6 = [55,56,81,82,83,84,85,96]
    
       lp i  = (elecs i) - (sorb i)


    
-- isHeavyAtom
-- Returns True for anything that has mass and is not a Hydrogen Atom
isHeavyAtom :: Atom -> Bool
isHeavyAtom a = case a of
   Element i _ _ _ -> i > 1
   LonePair {} -> False
   Electron {} -> False
   Unfilled {} -> False
   
-- isElement
-- Returns TRUE if atom is a "real" element
isElement :: Atom -> Bool
isElement a = case a of
  Element i _ _ _ -> True
  LonePair {} -> False
  Electron {} -> False
  Unfilled {} -> False

-- numberOfBonds
-- Returns number of covalent connections to other atoms in the molecule 
-- graph (i.e. one sigma and two pi bonds count as a 'one' bond)
numberOfBonds :: Atom -> Integer
numberOfBonds a = fromIntegral $ length $ getBondList a


numberOfBondsToAtoms :: Atom -> Integer
numberOfBondsToAtoms a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = fromIntegral $ length $ List.filter (==True) $ List.map isAnyBondToAtom b 

numberOfAromaticBondsToAtoms :: Atom -> Integer
numberOfAromaticBondsToAtoms a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = fromIntegral $ length $ List.filter (==True) $ List.map isAromaticBondToAtom b    

-- currentValence
-- Aromatic bonds count as ONE bond, no matter how many there are.  Aromatic atoms that are
-- incorrectly indicated but not part of a ring system are effectively "radical" system, but
-- will not be explicitely depicted as such in the data structure.  Aromatic atoms that are
-- not directly connected to any other aromatic atom will be considered ?? what exactly???
currentValence :: Atom ->  Integer
currentValence a = undefined

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

removeClosureMarker :: Atom -> Integer -> Atom
removeClosureMarker a i = a{markerSet = (Set.delete deleteMarker $ markerSet a)} 
    where deleteMarker = Closure i Single

getMatchingClosureNumber :: Atom -> Atom -> Maybe Integer
getMatchingClosureNumber a1 a2 = firstCommonMarker
    where markers1 = fst $ Set.partition isClosure (markerSet a1)
          markers2 = fst $ Set.partition isClosure (markerSet a2)
          intersectionSet = Set.intersection markers1 markers2
          isClosure mk  = case mk of Closure {} -> True ; _ -> False
          firstCommonMarker = if (Set.null intersectionSet) 
                              then Nothing 
                              else Just $ labelNumber (Set.findMin intersectionSet)

getBondList :: Atom -> [Bond]
getBondList a = case a of
    Element _ _ b _ -> b
    LonePair b _ -> b
    Electron b _ -> b
    Unfilled b _ -> b

instance Show Bond where
    show b = case b of
        Sigma atom -> "Sigma"
        Pi atom -> "Pi"
        Aromatic atom -> "Aromatic"
        Delta atom -> "Delta"
        Hbond atom -> "Hydrogen"
        Ionic atom -> "Ionic"
        Antibond atom -> "AntiBond"

-- This is really ugly, but need to equate closure markers easily, disregarding bond info.
-- This is because closure bond type only needs to be defined on one end of the molecule,
-- and therefore might not match the other closure atom in a valid smile.
instance Eq Marker where
    a == b = case a of 
        Closure {labelNumber=l1, bondType=b1} -> case b of 
            Closure {labelNumber=l2, bondType=b2} -> if (l1 == l2) then True else False
            _ -> False
        Class {classNumber=l1} -> case b of 
            Class {classNumber=l2} -> if (l1 == l2) then True else False    
            _ -> False 
        Chiral {chirality=l1} -> case b of 
            Chiral {chirality=l2} -> if (l1 == l2) then True else False
            _ -> False
        GeoIsomer {geoIsomer=l1} -> case b of 
            GeoIsomer {geoIsomer=l2} -> if (l1 == l2) then True else False
            _ -> False
        AromaticAtom -> case b of 
            AromaticAtom -> True
            _ -> False
        Traversed -> case b of 
            Traversed -> True   
            _ -> False 
        Substructure {substructureNumber=l1} -> case b of 
            Substructure {substructureNumber=l2} -> if (l1 == l2) then True else False
            _ -> False
        ValenceError {valenceError=l1} -> case b of 
            ValenceError {valenceError=l2} -> if (l1 == l2) then True else False
            _ -> False
        InRing {ringNumber=l1} -> case b of 
            InRing {ringNumber=l2} -> if (l1 == l2) then True else False  
            _ -> False   
        Skip -> case b of 
            Skip -> True
            _ -> False
        Comment {comment=l1} -> case b of 
            Comment {comment=l2} -> if (l1 == l2) then True else False    
            _ -> False
        Null -> case b of 
            Null -> True 
            _ -> False
        _ -> case b of
            _  -> False      
        

instance Show Atom where
    show a = case a of
          Element i _ b m ->  name i ++ show b ++ " " ++ show m
          LonePair {} -> "LP"
          Electron {} -> "•"
          Unfilled {} -> ""
          where name b = fromJust $ Map.lookup b atomicSymbols
{-
-- Instance Eq and instance Ord are going to be where all the action is.
-- Everything broken until then!
instance Eq Atom where
    a == b = True

instance Eq Bond where
    a == b = True

-}
    
instance Ord Atom where
    a > b = True
    a < b = False
