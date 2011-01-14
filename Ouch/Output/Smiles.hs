{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Output.Smiles
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

module Ouch.Output.Smiles (
     SmiWriterState (..)
   , writeSmiles
   ) where

import {-# SOURCE #-} Ouch.Property.Extrinsic.Fingerprint
import {-# SOURCE #-} Ouch.Structure.Molecule
import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Enumerate.Method
import Ouch.Property.Ring
import Ouch.Data.Atom
import Ouch.Data.Bond
import Ouch.Structure.Marker
import Data.Maybe
import Data.Set as Set
import Data.Map as Map
import Data.List as List


-- | Stores state information used while writing SMILES
data SmiWriterState = SmiWriterState
  { closureMap :: Map Int Int             -- ^ Information on how to match closures
                                          --   Key is closure label
                                          --   Value is atom number that was assigned the closure
  , smilogger  :: [String]                -- ^ Free-form logging string
  } deriving (Show)

-- | Look at the SMILES writer state and generate the required closure label
getClosureLabel :: SmiWriterState         -- ^ The state
                -> Int                    -- ^ The atom number being written
                -> Int                    -- ^ The atom number to connect to
                -> (Int, SmiWriterState)  -- ^ The closure label to use and new state
getClosureLabel state@SmiWriterState {closureMap=m} atomNum connectNum = let
  connectionMap = Map.filter (==connectNum) m
  isPair = 1 == Map.size connectionMap

  -- If there is already a pair, find it and remove its closure
  withClosure = fst $ Map.findMax connectionMap
  newMap = Map.delete withClosure m

  -- If not, find a new closure
  withClosure' = nextNum m
  newMap' = Map.insert withClosure' atomNum m
  output | isPair = (withClosure, state {closureMap=newMap})
         | otherwise = (withClosure', state {closureMap=newMap'})
  in output

-- | Generate the next closure number
nextNum :: (Map Int Int) -> Int
nextNum m | 0 == Map.size m = 1
          | otherwise = Set.findMin $ Set.difference ks ks'
  where ks = Set.fromList $ Map.keys m
        ks' = Set.fromList [1..maxKey]
        maxKey = (+1) $ Set.findMax ks


{------------------------------------------------------------------------------}
-- | Writes atom and bond information from position in a path
writeAtom :: PGraph -> Int -> String
writeAtom g i = let
  mol = molecule g
  atom = fromJust $ getAtomAtIndex mol ((vertexList g)!!i)
  hasExplicitH = Set.member (ExplicitHydrogen 0) (atomMarkerSet atom)
  hasCharge = Set.member (Charge 0) (atomMarkerSet atom)
  isAnyCharge charge = case charge of Charge {} -> True; _ -> False
  h | hasExplicitH = numberH $ Set.findMax $ Set.filter (== (ExplicitHydrogen 0)) (atomMarkerSet atom)
    | otherwise = 0
  c | hasCharge = charge $ Set.findMax $ Set.filter isAnyCharge (atomMarkerSet atom)
    | otherwise = 0
  symbol = atomicSymbolForAtom atom
  hString | not hasExplicitH = ""
          | h == 1 = "H"
          | h == 0 = ""
          | otherwise = "H" ++ (show h)
  cString | not hasCharge = ""
          | c == 1 = "+"
          | c == (-1) = "-"
          | c > 0 = "+" ++ (show c)
          | otherwise = (show c) -- Will show minus sign by default
  output | hasExplicitH || hasCharge = "[" ++ symbol  ++ hString ++ cString ++ "]"
         | otherwise = symbol
  in output


{------------------------------------------------------------------------------}
-- | Writes the SMILES string for a given path with a provided atom rendering function

{------------------------------------------------------------------------------}
-- | Writes the SMILES string for a given molecule
writeSmiles :: Molecule -> String
writeSmiles m = writeCanonicalPathWithStyle writeAtom m'
  where m':_ = [m] >#> stripMol


