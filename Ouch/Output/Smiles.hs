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
import Ouch.Structure.Marker hiding (position)
import Data.Maybe
import Data.Set as Set
import Data.Map as Map
import Data.List as List


newtype Logger = Logger {logger :: [String]} deriving (Show)
newtype Pair = Pair {pair :: (Int, Int)} deriving (Show, Ord, Eq)

-- | Stores state information used while writing SMILES
data SmiWriterState = SmiWriterState
  { smiString  :: String
  , style      :: SmiStyle
  , closureMap :: Map Int Pair
  , position   :: Int
  , traversing :: PGraph
  , tranversed :: [PGraph]
  , smiLogger  :: Logger
  }


data SmiStyle = SmiStyle
  { atomStyle :: Atom -> String
  , bondStyle :: Bond -> String
  }

instance Show SmiWriterState where
  show state = smiString state


startAtom :: Pair -> Int
startAtom p = fst $ pair p


endAtom :: Pair -> Int
endAtom p = snd $ pair p


logString :: Logger -> String -> Logger
logString l s =  Logger $ s:(logger l)


smiStart :: Molecule -> SmiWriterState
smiStart m = SmiWriterState
  { smiString  = ""
  , style      = writeAtom
  , closureMap = Map.empty
  , position   = 0
  , traversing = longestLeastPath m
  , tranversed = []
  , smiLogger  = Logger []
  }



-- | Look at the SMILES writer state and generate the required closure label
getClosureLabel :: SmiWriterState         -- ^ The state
                -> Int                    -- ^ The atom number being written
                -> Int                    -- ^ The atom number to connect to
                -> (Int, SmiWriterState)  -- ^ The closure label to use and new state
getClosureLabel state@SmiWriterState {closureMap=m, smiLogger=l} fromAtom toAtom = let
  newPair = Pair (fromAtom, toAtom)
  pairComplement = Pair (toAtom, fromAtom)
  isPair = 1 == (Map.size $ Map.filter (==pairComplement) m)
  isRedundant = 1 == (Map.size $ Map.filter (==newPair) m)

  -- If there is already a pair, find it and remove its closure
  withMatch = fst $ Map.findMax $ Map.filter (==pairComplement) m
  newMapMatch = Map.delete withMatch m

  withRedundant = fst $ Map.findMax $ Map.filter (==newPair) m
  newMapRedundant = m
  logError = logString l $ "Closure already requested for pair: " ++ (show newPair)

  -- If not, find a new closure
  withNew = nextNum m
  newMapNew = Map.insert withNew newPair m
  output | isPair      = (withMatch,     state {closureMap=newMapMatch})
         | isRedundant = (withRedundant, state {closureMap=newMapRedundant, smiLogger=logError})
         | otherwise   = (withNew,       state {closureMap=newMapNew})
  in output



-- | Generate the next closure number
nextNum :: Map Int Pair -> Int
nextNum m | 0 == Map.size m = 1
          | otherwise = Set.findMin $ Set.difference ks' ks
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
  h | hasExplicitH = numberH $ Set.findMax $ Set.filter (==(ExplicitHydrogen 0)) (atomMarkerSet atom)
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
-- | Writes the SMILES string for a given molecule
writeSmiles :: Molecule -> String
writeSmiles m = writeCanonicalPathWithStyle writeAtom m'
  where m':_ = [m] >#> stripMol


