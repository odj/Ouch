{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Property.Ring
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

module Ouch.Property.Ring
  (
    PGraph(..)
  , pathLength
  , inPath
  , hasOverlap
  ) where



import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import {-# SOURCE #-} Ouch.Structure.Molecule
import Data.Either
import Data.Maybe
import Data.Map as Map
import Data.List as List
import Data.Set as Set


data PGraph = PGraph {molecule::Molecule, vertexList::[Int]}

hasOverlap :: PGraph -> PGraph -> Bool
hasOverlap p1 p2 = let
  l1 = vertexList p1
  l2 = vertexList p2
  intersectList = List.intersect l1 l2
  in List.length intersectList > 0

inPath :: PGraph -> Int -> Bool
inPath pg i = List.elem i $ vertexList pg

pathLength :: PGraph -> Integer
pathLength p = fromIntegral $ List.length $ vertexList p

instance Show PGraph where
  show p = (show $ vertexList p) ++ "\n"

{--
createPGraph :: Molecule -> PGraph
createPGraph m = PGraph m' v
    where atoms = Map.filter isHeavyAtom $ atomMap m
          v = Map.keys atoms
          newAtomMap = Map.adjust
          addPEdge a = (\a )(PEdge i [])


                  | PEdge {reaches::Integer, pathList::[Integer]}

markAtom :: Atom -> AtomMarker -> Atom
--}










