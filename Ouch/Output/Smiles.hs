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

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}

module Ouch.Output.Smiles (
     SmiWriterState (..)
   , writeSmiles
   ) where

import {-# SOURCE #-} Ouch.Property.Extrinsic.Fingerprint
import {-# SOURCE #-} Ouch.Structure.Molecule
import Ouch.Property.Graph
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
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V


newtype Logger = Logger {logger :: [String]} deriving (Show)
newtype Pair = Pair {pair :: (Int, Int)} deriving (Show, Ord, Eq)



{------------------------------------------------------------------------------}
-- | Writes the SMILES string for a given molecule
writeSmiles :: Molecule -> String
writeSmiles = smiData . runSS . smiStart



-- | Stores state information used while writing SMILES
data SmiWriterState a = SmiWriterState
  { smiData  :: !a                -- ^ The SMILES String
  , preprocessMethod :: Maybe Method
  , selectStrategy :: [(Molecule -> Int -> Int)]
  , searchStrategy :: [(Molecule -> Int -> Int)]
  , ordStrategy :: [(PGraph -> PGraph -> Int -> Ordering)]
  , style      :: SmiStyle a              -- ^ The Style to render with
  , closureMap :: Map Int Pair          -- ^ Closure map for matching rings
  , position   :: !Int                   -- ^ The position on the current path
  , traversing :: !PGraph                -- ^ The path currently being rendered
  , traversed  :: ![PGraph]              -- ^ Paths that have already been rendered
  , smiLogger  :: Logger                -- ^ Logger for errors and warnings
  }




-- | Concatonates two SMILES Writer states where the second state is
-- typically a rendered substructure of the first state.
(<+>) :: SmiWriterState String -> SmiWriterState String -> SmiWriterState String
(<+>) !s1 !s2 = s1 { smiData=(smiData s1) ++ (smiData s2)  
                 , smiLogger=Logger $ (logger $ smiLogger s1) ++ (logger $ smiLogger s2)
                 , traversed=List.filter (/= traversing s1) $ (traversing s2):((traversed s1) `List.union` (traversed s2))
                 , closureMap=closureMap s2
                 }

(<@>) :: Map Int Pair -> Map Int Pair -> Map Int Pair
(<@>) cMap sub_cMap = let
  removed = Map.foldWithKey (\k a acc -> case Map.lookup k sub_cMap of 
                              Nothing -> Map.delete k acc
                              Just a  -> acc) cMap cMap
  {-in Map.union cMap sub_cMap-}
  in Map.union removed sub_cMap

(>>>) :: SmiWriterState a -> (SmiWriterState a -> SmiWriterState a) -> SmiWriterState a
(>>>) !state stateFun = stateFun state

-- | Data type to pass higher-order functions to writer
data SmiStyle a = SmiStyle 
  { atomStyle :: Atom -> a
  , bondStyle :: NewBond -> a
  }


-- | Show instance for debugging writer state
instance Show (SmiWriterState String) where
  show state = "SMILES: "     ++ smiData state           ++ "\n"
            ++ "POSITION: "   ++ (show $ position state)   ++ "\n"
            ++ "TRAVERSING: " ++ (show $ traversing state) ++ "\n"
            ++ "TRAVERSED:\n" ++ (show $ traversed state)  ++ "\n"
            ++ "CLOSURES: "   ++ (show $ closureMap state) ++ "\n"
            ++ "LOGS: "       ++ (show $ smiLogger state)  ++ "\n"


-- | Convenience method for newtype Pair
startAtom :: Pair -> Int
startAtom !p = fst $ pair p


-- | Convenience method for newtype Pair
endAtom :: Pair -> Int
endAtom !p = snd $ pair p


-- | Adds a new string to the newtype Logger
logString :: Logger -> String -> Logger
logString l s =  Logger $ s:(logger l)


-- | Creates a new state from an existing state
-- to be treated as a substructure.
smiNewSub :: SmiWriterState String -> PGraph -> SmiWriterState String
smiNewSub !state !path = 
  state  { smiData = "("
         , position = 0
         , traversing = path
         , smiLogger = Logger []
         }


-- | Creates an initial state from a Molecule with settings
-- that should give typically good results
smiStart :: Molecule -> SmiWriterState String
smiStart !m = loadWriter fastCanonicalWriter m

loadWriter :: SmiWriterState a -> Molecule -> SmiWriterState a
loadWriter !state !m = let
  processedMolecule = head $ [m] >#> (preprocessMethod state)
  initialPath = initialPathForStrategy processedMolecule 
                                       (selectStrategy state) 
                                       (searchStrategy state) 
                                       (ordStrategy state) 
  in state { position = 0
           , traversing = initialPath
           , traversed = [] 
           , smiLogger = Logger []
           }


emptySmiState = SmiWriterState 
  { smiData = undefined
  , preprocessMethod = Nothing
  , selectStrategy = []
  , searchStrategy = []
  , ordStrategy = []
  , closureMap = Map.empty
  , style = undefined
  , position   = 0
  , traversing = emptyPath
  , traversed  = []
  , smiLogger  = Logger []
  }

fastCanonicalWriter = SmiWriterState 
  { smiData  = ""
  , preprocessMethod = stripMol
  , selectStrategy = [ bottomVertextStrategy
                     , topBondStrategy
                     , atomTypeStrategy
                     ]
  , searchStrategy = [ topVertexStrategy
                     , atomTypeStrategy 
                     , topBondStrategy
                     ]
  , ordStrategy = []
  , style      = SmiStyle { atomStyle=writeAtom
                          , bondStyle=writeBond
                          }
  , closureMap = Map.empty
  , position   = 0
  , traversing = emptyPath
  , traversed  = []
  , smiLogger  = Logger []
  }

-- | Advances the rendering of the state by one atom along the path
advanceSS :: SmiWriterState String -> SmiWriterState String
advanceSS !state= state >>> renderRootBondSS
                        >>> renderAtomSS
                        >>> renderClosuresSS
                        >>> advanceClosuresSS
                        >>> renderSubpathsSS
                        >>> renderBondSS
                        >>> advancePosSS

-- | Advances the state position by one atom along the path without rendering
advancePosSS :: SmiWriterState a -> SmiWriterState a
advancePosSS !s = s { position = (position s) + 1 }


-- | Logs a string to the state logger
logSS :: SmiWriterState a -> String -> SmiWriterState a
logSS s str = s { smiLogger = logString (smiLogger s) str }


-- | Renders the atom at the current state position according to the state's style
renderAtomSS :: SmiWriterState String -> SmiWriterState String
renderAtomSS !s@SmiWriterState {traversing=path, style=st, position=p_i} = let
  atom = fromJust $ getAtomAtIndex (molecule path) (pathIndex path p_i)
  in forceRender s (atomStyle st $ atom)


-- | Renders the bond at the current state position according to the state's style
renderBondSS :: SmiWriterState String -> SmiWriterState String
renderBondSS !s@SmiWriterState {traversing=path, style=st, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i
  index' = pathIndex path (p_i+1)
  atom = fromJust $ getAtomAtIndex mol index
  hasNext = p_i + 1 < (fromInteger $ pathLength path)
  bond | hasNext = bondStyle st $ (bondBetweenIndices mol index index' )
       | otherwise = ""
  in forceRender s bond

renderRootBondSS :: SmiWriterState String -> SmiWriterState String
renderRootBondSS !state@SmiWriterState {traversing=path, style=st, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i
  r = root path
  output | U.length r == 0 = ""
         | otherwise = renderWithRoot (U.head r)
  renderWithRoot m_i | p_i /= 0 = ""
                     | List.head (smiData state) /= '(' = ""
                     | otherwise = bondStyle st $ (bondBetweenIndices mol index m_i)
  in forceRender state output


pairBondType :: Molecule
         -> Pair
         -> NewBond
pairBondType !mol !pair = bondBetweenIndices mol (startAtom pair) (endAtom pair)

-- | Renders the closures at the current state position according to the state's style
renderClosuresSS :: SmiWriterState String -> SmiWriterState String
renderClosuresSS !state@SmiWriterState {traversing=path, closureMap=cMap, style=st} = let
  mol = molecule path
  renderClosure (n, pair) | Map.member n cMap = renderLabel n
                          | otherwise = (bondStyle st $ pairBondType mol pair)
                                     ++ (renderLabel n)
  renderLabel n | (length $ show n) > 1 = '%':(show n)
                | otherwise = show n
  pairs = findClosuresSS $ (state >>> renderSubpathsSS) {closureMap = cMap}
  closures = fst $ getClosureLabels cMap pairs
  zipped = List.sortBy (\a b -> compare (fst b) (fst a)) $ List.zip closures pairs
  in forceRender state (List.foldr (\z acc -> acc ++ (renderClosure z)) "" zipped)


-- | Advances the closure state
advanceClosuresSS :: SmiWriterState String -> SmiWriterState String
advanceClosuresSS !state@SmiWriterState {closureMap=cMap} = state {closureMap=newMap}
  where newMap = snd $ getClosureLabels cMap (findClosuresSS $ (state >>> renderSubpathsSS) {closureMap = cMap})


-- | Generates a list of closure Pairs required at this path position, but does
-- not advance the state
findClosuresSS :: SmiWriterState a -> [Pair]
findClosuresSS state@SmiWriterState {traversing=path, traversed=paths, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i

  -- The span set makes sure that we do not add a closure label to the atom we
  -- were just on, or the atom we are about to go to.
  spanSet | atStartSS state && atLastSS state = Set.empty
          | atStartSS state = Set.fromList [pathIndex path (p_i + 1)]
          | atLastSS state = Set.fromList [pathIndex path (p_i - 1)]
          | otherwise = Set.fromList [pathIndex path (p_i + 1), pathIndex path (p_i - 1)]

  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex mol index

  -- All indices we have seen or are going to see on this path
  pathIndexSet = List.foldr (\a acc -> Set.union acc $ vToSet $ vertexList a)
                            Set.empty $ (path:paths) 

  -- If we have a root path, remove the first index from consideration.  It is 
  -- where we started from.  But ONLY if we are at the beginning of our path. 
  -- Otherwise, we mess up spiro systems.
  pathIndexSet_noRoot | U.length (root path) == 0 = pathIndexSet
                      | p_i == 0 = Set.delete (U.head $ root path) pathIndexSet
                      | otherwise = pathIndexSet

  -- Remove the heads from any visited paths IF their root bond is the current position
  pathIndexSet_noHeads = List.foldr (\a acc -> Set.delete (U.head $ vertexList a) acc) pathIndexSet_noRoot rootPaths
    where rootPaths = List.filter (\p -> equalsRoot p) paths
          equalsRoot p | U.length (root p) == 0 = False
                       | otherwise = (==) index $ U.head $ root p
  validIndexSet = Set.difference pathIndexSet_noHeads spanSet
  closureIndexList = Set.toList $ Set.intersection bondIndexSet validIndexSet
  in List.map (\p -> Pair (index, p)) closureIndexList


-- | Tests to see if the state is at the end of its current path
atEndSS :: SmiWriterState a -> Bool
atEndSS state = (position state) == (fromIntegral $ pathLength $ traversing state)

-- | Tests to see if the state is at the end of its current path
atLastSS :: SmiWriterState a -> Bool
atLastSS state = (position state) + 1 == (fromIntegral $ pathLength $ traversing state)


-- | Tests to see if the state is at the beginning of its current path
atStartSS :: SmiWriterState a -> Bool
atStartSS state = (position state) == 0


runSS :: SmiWriterState String -> SmiWriterState String
runSS state | atEndSS state = if (head $ smiData state) == '('
                              then forceRender state ")"
                              else state
            | otherwise = state >>> advanceSS >>> runSS


forceRender :: SmiWriterState String -> String -> SmiWriterState String
forceRender state str = state `seq` state {smiData = (smiData state) ++ str}

renderSubpathsSS :: SmiWriterState String -> SmiWriterState String
renderSubpathsSS state@SmiWriterState {traversing=path, traversed=paths, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex mol index
  pathIndexSet = List.foldr (\a acc -> Set.union acc $ vToSet $ vertexList a)
                            Set.empty (path:paths)
  validIndexList = Set.toList $ Set.difference bondIndexSet pathIndexSet
  nextIndex = head $ List.sort $ validIndexList
  pathIndexList = Set.toList pathIndexSet
  branchPaths = List.sort $ List.map (\a -> longestLeastAnchoredPath path {vertexList=(U.fromList pathIndexList)
                                                                         , root=(U.fromList (index:pathIndexList))} a) validIndexList
  nextPath = head branchPaths
  branchState = runSS $ state { smiData = "("
                              , traversing=nextPath
                              , traversed=(path:paths)
                              , position=0
                              } 
  newState = state <+> branchState 
  output | List.length branchPaths > 0 = newState <+> renderSubpathsSS newState {smiData = "" }
         | otherwise = state
  in output


findSubpathsSS :: SmiWriterState String -> [SmiWriterState String]
findSubpathsSS state@SmiWriterState {traversing=path, traversed=paths, position=p_i} = let
  mol = molecule path
  index = pathIndex path p_i
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex mol index
  pathIndexSet = List.foldr (\a acc -> Set.union acc $ vToSet $ vertexList a)
                            Set.empty (path:paths)
  validIndexList = Set.toList $ Set.difference bondIndexSet pathIndexSet
  pathIndexList = Set.toList pathIndexSet
  branchPaths = List.sort $ List.map (\a -> longestLeastAnchoredPath path {vertexList=(U.fromList pathIndexList)
                                                                         , root=(U.fromList (index:pathIndexList))} a)
                         validIndexList
  in List.map (\p -> state {smiData = "("
                          , traversing=p
                          , traversed=(path:paths)
                          , position=0})  branchPaths


getClosureLabels :: Map Int Pair
                 -> [Pair]
                 -> ([Int], Map Int Pair)
getClosureLabels cMap pairs = List.foldl'
  (\(closures, cMap') pair -> (closures ++ [fst $ getClosureLabel cMap' closures pair]
                             , snd $ getClosureLabel cMap' closures pair) ) ([], cMap) pairs

-- | Look at the SMILES writer state and generate the required closure label
getClosureLabel :: Map Int Pair
                -> [Int]
                -> Pair
                -> (Int, Map Int Pair)
getClosureLabel cMap exclude pair = let
  pairComplement = Pair (endAtom pair, startAtom pair)
  isPair = 1 == (Map.size $ Map.filter (==pairComplement) cMap)
  isRedundant = 1 == (Map.size $ Map.filter (==pair) cMap)

  -- If there is already a pair, find it and remove its closure
  withMatch = fst $ Map.findMax $ Map.filter (==pairComplement) cMap
  newMapMatch = Map.delete withMatch cMap

  -- The original pair was already there.  This portends bad things will happen
  -- but for now do nothing.
  withRedundant = fst $ Map.findMax $ Map.filter (==pair) cMap
  newMapRedundant = cMap

  -- If not, find a new closure
  withNew = nextNum cMap exclude
  newMapNew = Map.insert withNew pair cMap
  output | isPair      = (withMatch, newMapMatch)
         | isRedundant = (withRedundant, newMapRedundant)
         | otherwise   = (withNew, newMapNew)
  in output



-- | Generate the next closure number
nextNum :: Map Int Pair   -- ^ The map to check against
        -> [Int]          -- ^ Keys to exclude from consideration
        -> Int            -- ^ The next key
nextNum m exclude | 0 == Set.size ks = 1
                  | otherwise = Set.findMin $ Set.difference ks' ks
  where ks = Set.fromList $ (exclude ++ Map.keys m)
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


-- | Basic rendering of bond information for SMILES.  Similar to Show.
writeBond :: NewBond  -- ^ The bond to render
          -> String   -- ^ The rendered string
writeBond nb = case nb of
  Single -> ""
  Double -> "="
  Triple -> "#"
  NoBond -> "."

db m = putStrLn $ debugShow $ head ([m] >#> stripMol)

