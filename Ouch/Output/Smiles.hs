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
  , traversed  :: [PGraph]
  , smiLogger  :: Logger
  }

(<+>) :: SmiWriterState -> SmiWriterState -> SmiWriterState
(<+>) s1 s2 = s1 { smiString=(smiString s1) ++ (smiString s2)
                 , smiLogger=Logger $ (logger $ smiLogger s1) ++ (logger $ smiLogger s2)
                 , traversed=(traversed s1 ++ traversed s2)
                 , closureMap=(Map.union (closureMap s1) (closureMap s2))
                 }

data SmiStyle = SmiStyle
  { atomStyle :: Atom -> String
  , bondStyle :: NewBond -> String
  }

instance Show SmiWriterState where
  show state = "SMILES: " ++ smiString state ++ "\n"
            ++ "POSITION: " ++ (show $ position state) ++ "\n"
            ++ "TRAVERSING: " ++ (show $ traversing state)
            ++ "LOGS: " ++ (show $ smiLogger state) ++ "\n"


startAtom :: Pair -> Int
startAtom p = fst $ pair p


endAtom :: Pair -> Int
endAtom p = snd $ pair p


logString :: Logger -> String -> Logger
logString l s =  Logger $ s:(logger l)

smiNewSub :: SmiWriterState -> PGraph -> SmiWriterState
smiNewSub state@SmiWriterState { style=st
                               , traversed=tv
                               , closureMap=cmap} path =
  SmiWriterState { smiString = "("
                 , style = st
                 , closureMap = cmap
                 , position = 0
                 , traversing = path
                 , traversed = tv
                 , smiLogger = Logger []
                 }


smiStart :: Molecule -> SmiWriterState
smiStart m = SmiWriterState
  { smiString  = ""
  , style      = SmiStyle { atomStyle=writeAtom
                          , bondStyle=writeBond
                          }
  , closureMap = Map.empty
  , position   = 0
  , traversing = longestLeastPath $ head ([m] >#> stripMol)
  , traversed  = []
  , smiLogger  = Logger []
  }


advanceSS :: SmiWriterState -> SmiWriterState
advanceSS state = let
  render = (renderAtomSS state) ++ (renderBondSS state)
  renderedSubpaths = List.map runSS $ findSubpathsSS state
  advancedState = forceRenderSS (advancePosSS $ advanceClosuresSS $ state) render
  in List.foldl (<+>) advancedState renderedSubpaths


advancePosSS :: SmiWriterState -> SmiWriterState
advancePosSS s = s { position = (position s) + 1 }


logSS :: SmiWriterState -> String -> SmiWriterState
logSS s str = s { smiLogger = logString (smiLogger s) str }


renderAtomSS :: SmiWriterState -> String
renderAtomSS s@SmiWriterState {traversing=path, style=st, position=p_i} = let
  atom = fromJust $ getAtomAtIndex (molecule path) (pathIndex path p_i)
  in atomStyle st $ atom


renderBondSS :: SmiWriterState -> String
renderBondSS s@SmiWriterState {traversing=path, style=st, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i
  index' = pathIndex path (p_i+1)
  atom = fromJust $ getAtomAtIndex mol index
  hasNext = p_i + 1 < (fromInteger $ pathLength path)
  bond | hasNext = bondStyle st $ (bondBetweenIndices mol index index' )
       | otherwise = ""
  in bond

advanceClosuresSS :: SmiWriterState -> SmiWriterState
advanceClosuresSS state = state

findClosuresSS :: SmiWriterState -> [Pair]
findClosuresSS state = undefined



atEndSS :: SmiWriterState -> Bool
atEndSS state = (position state) == (fromIntegral $ pathLength $ traversing state)


runSS :: SmiWriterState -> SmiWriterState
runSS state | atEndSS state = if (head $ smiString state) == '('
                              then forceRenderSS state ")"
                              else state
            | otherwise = runSS $ advanceSS state


forceRenderSS :: SmiWriterState -> String -> SmiWriterState
forceRenderSS state str = state {smiString = (smiString state) ++ str}



findSubpathsSS :: SmiWriterState -> [SmiWriterState]
findSubpathsSS state@SmiWriterState {traversing=path, traversed=paths, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex mol index
  pathIndexSet = List.foldr (\a acc -> Set.union acc $ Set.fromList $ vertexList a)
                            Set.empty (path:paths)
  validIndexList = Set.toList $ Set.difference bondIndexSet pathIndexSet
  pathIndexList = Set.toList pathIndexSet
  branchPaths = List.sort $ List.map (\a -> longestLeastAnchoredPath path {vertexList=pathIndexList, root=(Just index)} a)
                         validIndexList
  in List.map (\p -> state {smiString = "("
                          , traversing=p
                          , traversed=(path:paths)
                          , position=0})  branchPaths


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
writeAtom :: Atom -> String
writeAtom atom = let
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
writeSmiles m = undefined



-- | Basic rendering of bond information for SMILES.  Similar to Show.
writeBond :: NewBond  -- ^ The bond to render
          -> String   -- ^ The rendered string
writeBond nb = case nb of
  Single -> ""
  Double -> "="
  Triple -> "#"
  NoBond -> "."


db m = putStrLn $ debugShow $ head ([m] >#> stripMol)
 {--

{------------------------------------------------------------------------------}
-- | Writes the SMILES string for a given path with a provided atom rendering function
writePath :: (PGraph -> Int -> String)    -- ^ The method used to render atoms to text
          -> [PGraph]                     -- ^ The list of subgraphs we've already traversed
          -> PGraph                       -- ^ The subgraph we are currently traversing
          -> Int                          -- ^ The position in our current subgraph
          -> Bool                         -- ^ Are we part of a SMILES substructure
          -> String                       -- ^ The string being rendered
writePath style gx g@PGraph {molecule=m, vertexList=l} i subStructure = let
  s = style g i
  endOfPath = i >= (fromInteger $ pathLength g)
  hasNext = i + 1 < (fromInteger $ pathLength g)
  nb | hasNext = writeBond (bondBetweenIndices m  (l!!i) $ l!!(i+1) )
     | otherwise = ""
  output | endOfPath && subStructure = ")"
         | endOfPath = ""
         | otherwise =  s ++ writeSubpath style gx g i
                          ++ nb ++ writePath style gx g (i+1) subStructure
  in output


{------------------------------------------------------------------------------}
-- | Writes the SMILES strings for all subpaths (if any exist) at position i in a
-- given path g, excluding travesal through any atoms in the paths gx
writeSubpath :: (PGraph -> Int -> String) -- ^ The method used to render atoms to text
             -> [PGraph]                  -- ^ The list of subgraphs we've already traversed
             -> PGraph                    -- ^ The subgraph we are currently traversing
             -> Int                       -- ^ The position in our current subgraph
             -> String                    -- ^ The string being rendered
writeSubpath outputStyle gx g@PGraph {molecule=m, vertexList=l} i = let
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex m (l!!i)
  pathIndexSet = List.foldr (\a acc -> Set.union acc $ Set.fromList $ vertexList a)
                            Set.empty (g:gx)
  validIndexList = Set.toList $ Set.difference bondIndexSet pathIndexSet
  pathIndexList = Set.toList pathIndexSet
  branchPaths = List.map (\a -> longestLeastAnchoredPath g {vertexList=pathIndexList, root=(Just $ pathIndex g i)} a)
                         validIndexList
  nextBranch =  findLongestLeastPath branchPaths 0
  s = outputStyle g i
  output | (pathLength nextBranch) > 0 =  "(" ++ writePath outputStyle (g:gx) nextBranch 0 True
                                              ++ writeSubpath outputStyle (nextBranch:gx) g i
         | otherwise =  ""
  in output



{------------------------------------------------------------------------------}
-- | Writes the SMILES string for a given molecule
writeCanonicalPath :: Molecule -> String
writeCanonicalPath m = writeCanonicalPathWithStyle writeAtomOnly m


{------------------------------------------------------------------------------}
-- | Writes the SMILES string for a given molecule with specified atom rendering function
writeCanonicalPathWithStyle :: (PGraph -> Int -> String)  -- ^ The method used to render atoms to text
                            -> Molecule                   -- ^ The molecule to process
                            -> String                     -- ^ The output string
writeCanonicalPathWithStyle style m = let
  backbone = longestLeastPath m
  in writePath style [] backbone 0 False





{------------------------------------------------------------------------------}
-- | Writes atom information from position i in a path.   A debugging render
-- function.
writeAtomOnly :: PGraph   -- ^ The path
              -> Int      -- ^ The path position to render
              -> String   -- ^ The output String
writeAtomOnly g i = let
  mol = molecule g
  atom = fromJust $ getAtomAtIndex mol ((vertexList g)!!i)
  in atomicSymbolForAtom atom

bip :: PGraph
         -> Int
         -> String
bip  _ _ = "*"

writeBip :: Molecule -> String
writeBip m = writeCanonicalPathWithStyle bip m'
  where m':_ = [m] >#> stripMol






  --}
