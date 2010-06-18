{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    AtomMarker - a module to manage atom markers
    
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

module Ouch.Structure.Marker (
      Chirality(..)
    , AtomMarker(..)
    , MoleculeMarker(..)
    , NewBond(..)
    , Geometry(..)
    ) where 

import {-# SOURCE #-} Ouch.Structure.Atom
        
{------------------------------------------------------------------------------}
data AtomMarker =   Label {labelNumber::Integer}   -- OUCH specific label maintained to match atom map held by Molecule
                  | Charge {charge::Integer}
                  | Position {position::(Double, Double, Double)}  -- x, y, z vector
                  | Closure {labelNumber::Integer, bondType::NewBond}
                  | Class {classNumber::Integer}
                  | Chiral {chirality::Chirality}
                  | GeoIsomer {geoIsomer::Geometry}
                  | AromaticAtom
                  | Traversed {order::Integer}
                  | ExplicitHydrogen {numberH::Integer}
                  | Substructure {substructureNumber::Integer}
                  | ValenceError {valenceError::String}
                  | InRing {ringNumber::Integer}
                  | Skip
                  | PGraph {reaches::Integer, pathList::[Integer]}
                  | Comment {comment::String}
                  | Null  -- This is a dummy value for functions that append marker list for simplicity.
                  deriving (Show)

data MoleculeMarker =   Info     {molMarker::String}
                      | Name     {molMarker::String}
                      | Warning  {molMarker::String} 
                      | MError {molMarker::String}
                      deriving (Eq)
                      

{------------------------------------------------------------------------------}
data Chirality = Levo | Dextro | UnknownChirality
   deriving (Show, Eq, Ord)


{------------------------------------------------------------------------------}
data Geometry = Cis {geometetryAtom::Atom} | Trans {geometetryAtom::Atom} | ProCis | ProTrans
   deriving (Show, Eq, Ord)

-- NewBond
-- This type communicates what 'should be' rather than what 'is'.  It is a marker used in
-- smiles parsing to help convey information before a bond is actually made.
data NewBond = Single | Double | Triple | NoBond deriving (Show, Eq, Ord)


{------------------------------------------------------------------------------}
{-------------------------------Typeclass Intances-----------------------------}
{------------------------------------------------------------------------------}


instance Show MoleculeMarker where
  show m = case m of
      Info s    -> "Info: " ++ s ++ "\n"
      Name s    -> "Name: " ++ s ++ "\n"
      Warning s -> "WARNING: " ++ s ++ "\n"
      MError s  -> "ERROR: " ++ s ++ "\n"



-- This is REALLY REALLY ugly, but need to equate closure markers easily in sets, disregarding bond info.
-- This is because closure bond type only needs to be defined on one end of the molecule,
-- and therefore might not match the other closure atom in a valid smile.
{------------------------------------------------------------------------------}


instance Eq AtomMarker where
   a == b = case a of 
      Position {position=l1} -> case b of 
          Position {position=l2} -> if (l1 == l2) then True else False
          _ -> False
      Charge {charge=l1} -> case b of 
          Charge {charge=l2} -> if (l1 == l2) then True else False
          _ -> False
      Closure {labelNumber=l1, bondType=b1} -> case b of 
          Closure {labelNumber=l2, bondType=b2} -> if (l1 == l2) then True else False
          _ -> False
      Class {classNumber=l1} -> case b of 
          Class {classNumber=l2} -> if (l1 == l2) then True else False    
          _ -> False 
      Chiral {} -> case b of 
          Chiral {} -> True
          _ -> False
      GeoIsomer {} -> case b of 
          GeoIsomer {} -> True
          _ -> False
      AromaticAtom -> case b of 
          AromaticAtom -> True
          _ -> False
      Traversed {} -> case b of 
          Traversed { }-> True   
          _ -> False 
      ExplicitHydrogen {} -> case b of 
          ExplicitHydrogen {}-> True   
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
      Label {labelNumber=l1} -> case b of
            Label {labelNumber=l2} -> if (l1 == l2) then True else False
            _ -> False



instance Ord AtomMarker where
   compare a b =  case a of 
          Label {}   -> case b of
              Label {}              -> EQ
              _                     -> GT
   
          Charge {} -> case b of 
             Label {}              -> LT
             Charge {}             -> EQ
             _                     -> GT

          Position {} -> case b of 
             Label {}              -> LT
             Charge {}             -> LT
             Position {}           -> EQ
             _                     -> GT
 
             
          Closure {labelNumber=l1} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {labelNumber=l2} -> a
                     where a | (l1 == l2) = EQ
                             | (l1 > l2)  = GT
                             | (l1 < l2)  = LT
                             | otherwise  = GT
              _                     -> GT

          Class {} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> EQ
              _                     -> GT

          Chiral {} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> EQ
              _                     -> GT

          GeoIsomer {} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> EQ
              _                     -> GT

          AromaticAtom -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> EQ
              _                     -> GT

          Traversed {order=l1} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {order=l2}  -> a
                where a | l1 > l2  = GT
                        | l1 < l2  = LT
                        | l1 == l2 = EQ
              _                     -> GT
 
          ExplicitHydrogen {} -> case b of 
                Label {}              -> LT
                Charge {}             -> LT
                Position {}           -> LT
                Closure {}            -> LT
                Class {}              -> LT
                Chiral {}             -> LT
                GeoIsomer {}          -> LT
                AromaticAtom          -> LT
                Traversed {}          -> LT
                ExplicitHydrogen {}   -> EQ
                _                     -> GT
 
          Substructure {substructureNumber=l1} -> case b of 
                  Label {}              -> LT
                  Charge {}             -> LT
                  Position {}           -> LT
                  Closure {}            -> LT
                  Class {}              -> LT
                  Chiral {}             -> LT
                  GeoIsomer {}          -> LT
                  AromaticAtom          -> LT
                  Traversed {}          -> LT
                  ExplicitHydrogen {}   -> LT
                  Substructure {substructureNumber=l2} -> a
                          where a | (l1 == l2) = EQ
                                  | (l1 > l2)  = GT
                                  | (l1 < l2)  = LT
                                  | otherwise  = EQ
                  _                     -> GT


          ValenceError {valenceError=l1} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {}          -> LT
              ExplicitHydrogen {}   -> LT
              Substructure {}       -> LT
              ValenceError {valenceError=l2} -> a
                      where a | (l1 == l2) = EQ
                              | (l1 > l2)  = GT
                              | (l1 < l2)  = LT
                              | otherwise  = EQ
              _                     -> GT

              
          InRing {ringNumber=l1} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {}          -> LT
              ExplicitHydrogen {}   -> LT
              Substructure {}       -> LT
              ValenceError {}       -> LT
              InRing {ringNumber=l2} -> a
                      where a | (l1 == l2) = EQ
                              | (l1 > l2)  = GT
                              | (l1 < l2)  = LT
                              | otherwise  = EQ
              _                     -> GT

 
          Skip -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {}          -> LT
              ExplicitHydrogen {}   -> LT
              Substructure {}       -> LT
              ValenceError {}       -> LT
              InRing {}             -> LT
              Skip   {}             -> EQ
              _                     -> GT

          PGraph {reaches=l1} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {}          -> LT
              ExplicitHydrogen {}   -> LT
              Substructure {}       -> LT
              ValenceError {}       -> LT
              InRing {}             -> LT
              Skip   {}             -> LT
              PGraph {reaches=l2}   -> a
                  where a   | (l1 == l2) = EQ
                            | (l1 > l2)  = GT
                            | (l1 < l2)  = LT
                            | otherwise  = EQ
              _                     -> GT


          Comment {comment=l1} -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {}          -> LT
              ExplicitHydrogen {}   -> LT
              Substructure {}       -> LT
              ValenceError {}       -> LT
              InRing {}             -> LT
              Skip   {}             -> LT
              PGraph {}             -> LT
              Comment {comment=l2} -> a
                      where a | (l1 == l2) = EQ
                              | (l1 > l2)  = GT
                              | (l1 < l2)  = LT
                              | otherwise  = EQ
              Null                  -> LT          

          Null -> case b of 
              Label {}              -> LT
              Charge {}             -> LT
              Position {}           -> LT
              Closure {}            -> LT
              Class {}              -> LT
              Chiral {}             -> LT
              GeoIsomer {}          -> LT
              AromaticAtom          -> LT
              Traversed {}          -> LT
              ExplicitHydrogen {}   -> LT
              Substructure {}       -> LT
              ValenceError {}       -> LT
              InRing {}             -> LT
              Skip   {}             -> LT
              PGraph {}             -> LT
              Comment {}            -> LT
              Null                  -> EQ

instance Ord MoleculeMarker where
  compare a b =  case a of 
      MError {}  -> case b of
          MError {}            -> EQ
          Name {}              -> GT
          Warning {}           -> GT
          Info {}              -> GT
      Name {}  -> case b of
          MError {}          -> LT
          Name {}              -> EQ
          Warning {}           -> GT
          Info {}              -> GT
      Warning {molMarker=l1}  -> case b of
          MError {}          -> LT
          Name {}              -> LT
          Warning {molMarker=l2} ->  a
              where a | (l1 == l2) = EQ
                        | (l1 > l2)  = GT
                        | (l1 < l2)  = LT
                        | otherwise  = EQ
          Info {}              -> GT

      Info {molMarker=l1}  -> case b of
          MError {}          -> LT
          Name {}              -> LT
          Warning {}           -> LT
          Info {molMarker=l2}   ->  a
              where a | (l1 == l2) = EQ
                      | (l1 > l2)  = GT
                      | (l1 < l2)  = LT
                      | otherwise  = EQ
          
          
          
          
          