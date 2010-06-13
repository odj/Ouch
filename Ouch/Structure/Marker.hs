{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Marker - a module to manage atom markers
    
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
    , Marker(..)
    , NewBond(..)
    , Geometry(..)
    ) where 

import {-# SOURCE #-} Ouch.Structure.Atom
        
{------------------------------------------------------------------------------}
data Marker =  Label {labelNumber::Integer}   -- OUCH specific label
              | Charge {charge::Integer}
              | Position {position::(Double, Double, Double)}  -- x, y, z vector
              | Closure {labelNumber::Integer, bondType::NewBond}
              | Class {classNumber::Integer}
              | Chiral {chirality::Chirality}
              | GeoIsomer {geoIsomer::Geometry}
              | AromaticAtom
              | Traversed {order::Integer}
              | Substructure {substructureNumber::Integer}
              | ValenceError {valenceError::String}
              | InRing {ringNumber::Integer}
              | Skip
              | Comment {comment::String}
              | Null  -- This is a dummy value for functions that append marker list for simplicity.
              deriving (Show)

{------------------------------------------------------------------------------}
data Chirality = Levo | Dextro | UnknownChirality
   deriving (Show, Eq, Ord)


{------------------------------------------------------------------------------}
data Geometry = Cis {geometetryAtom::Atom} | Trans {geometetryAtom::Atom}
   deriving (Show, Eq, Ord)

-- NewBond
-- This type communicates what 'should be' rather than what 'is'.  It is a marker used in
-- smiles parsing to help convey information before a bond is actually made.
data NewBond = Single | Double | Triple | NoBond deriving (Show, Eq, Ord)



-- This is really ugly, but need to equate closure markers easily, disregarding bond info.
-- This is because closure bond type only needs to be defined on one end of the molecule,
-- and therefore might not match the other closure atom in a valid smile.
{------------------------------------------------------------------------------}
instance Eq Marker where
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
      Chiral {chirality=l1} -> case b of 
          Chiral {chirality=l2} -> if (l1 == l2) then True else False
          _ -> False
      GeoIsomer {geoIsomer=l1} -> case b of 
          GeoIsomer {geoIsomer=l2} -> if (l1 == l2) then True else False
          _ -> False
      AromaticAtom -> case b of 
          AromaticAtom -> True
          _ -> False
      Traversed {} -> case b of 
          Traversed { }-> True   
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

instance Ord Marker where
   compare a b =  case a of 
          Charge {} -> case b of 
             Charge {} -> EQ
             _ -> LT
          Position {} -> case b of 
             Position {} ->  EQ
             _ -> LT 
          Closure {labelNumber=l1, bondType=b1} -> case b of 
              Closure {labelNumber=l2, bondType=b2} -> a
                  where a | (l1 == l2) = EQ
                          | (l1 > l2)  = GT
                          | (l1 < l2)  = LT
                          | otherwise  = EQ
              _ -> LT
          Class {} -> case b of 
              Class {} -> EQ
              _ -> LT
          Chiral {} -> case b of 
              Chiral {} ->  EQ
              _ -> LT
          GeoIsomer {} -> case b of 
              GeoIsomer {} -> EQ
              _ -> LT
          AromaticAtom -> case b of 
              AromaticAtom -> EQ
              _ -> LT
          Traversed {} -> case b of 
              Traversed {} -> EQ   
              _ -> LT 
          Substructure {substructureNumber=l1} -> case b of 
              Substructure {substructureNumber=l2} -> a
                    where a | (l1 == l2) = EQ
                            | (l1 > l2)  = GT
                            | (l1 < l2)  = LT
                            | otherwise  = EQ
              _ -> LT
          ValenceError {valenceError=l1} -> case b of 
              ValenceError {valenceError=l2} -> a
                    where a | (l1 == l2) = EQ
                            | (l1 > l2)  = GT
                            | (l1 < l2)  = LT
                            | otherwise  = EQ
              _ -> LT
          InRing {ringNumber=l1} -> case b of 
              InRing {ringNumber=l2} -> a
                    where a | (l1 == l2) = EQ
                            | (l1 > l2)  = GT
                            | (l1 < l2)  = LT
                            | otherwise  = EQ
              _ -> LT   
          Skip -> case b of 
              Skip -> EQ
              _ -> LT
          Comment {comment=l1} -> case b of 
              Comment {comment=l2} -> a
                    where a | (l1 == l2) = EQ
                            | (l1 > l2)  = GT
                            | (l1 < l2)  = LT
                            | otherwise  = EQ  
              _ -> LT
          Null -> case b of 
              Null -> EQ 
              _ -> LT
          _ -> case b of
              _  -> LT
