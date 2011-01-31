{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Property.Extrinsic.FingerPrint
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

module Ouch.Property.Extrinsic.Fingerprint (
    PGraph(..)
  , pathLength
  , inPath
  , hasOverlap
  , atomBits_OUCH
  , bondBits_OUCH
  , molBits_OUCH
  , molBits_N
  , molBits_ID
  , pathBits
  , findPaths
  , allPaths
  , pathIndex
  , findLongestLeastPath
  , longestPaths
  , longestLeastPath
  , longestLeastAnchoredPath
  , comparePaths
  , vToSet
  , (.||.)
  , (.|||.)


) where

import Data.Binary.Get as G
import Data.ByteString.Lazy as L
import Data.Binary.Builder as B
import Data.Bits
import Data.Word
import Ouch.Structure.Atom
import Ouch.Structure.Bond
import {-# SOURCE #-} Ouch.Structure.Molecule
import Ouch.Enumerate.Method
import Ouch.Property.Ring
import Ouch.Data.Atom
import Ouch.Data.Bond
import Ouch.Structure.Marker
import Ouch.Input.Smiles
import Data.Maybe
import Data.Set as Set
import Data.List as List hiding (intersect)
import Data.Map as Map
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Debug.Trace (trace)


data PGraph = PGraph { molecule   :: Molecule    -- ^ The molecule to apply the path
                     , vertexList :: U.Vector Int       -- ^ The path
                     , root       :: U.Vector Int       -- ^ The root mol paths to break recurson
                     }


intersect :: (U.Unbox a, Eq a) => U.Vector a -> U.Vector a -> U.Vector a
intersect v1 v2 = U.filter (\e -> U.elem e v2) v1

vToSet :: (U.Unbox a, Ord a) => U.Vector a -> Set a
vToSet v = Set.fromList $ U.toList v


hasOverlap :: PGraph -> PGraph -> Bool
hasOverlap p1 p2 = let
  l1 = vertexList p1
  l2 = vertexList p2
  intersectList = intersect l1 l2
  in U.length intersectList > 0

inPath :: PGraph -> Int -> Bool
inPath pg i = U.elem i $ vertexList pg

pathLength :: PGraph -> Integer
pathLength p = fromIntegral $ U.length $ vertexList p

instance Eq PGraph where
  (==) a b = (==) (vertexList a) (vertexList b)

instance Show PGraph where
  show p = (show $ vertexList p) ++ "\n"
        ++ "With root :" ++ (show $ root p)

instance Ord PGraph where
  compare a b = comparePaths a b


data Fingerprint = SimpleFingerPrint

{------------------------------------------------------------------------------}
atomBits_OUCH :: Molecule
              -> Atom
              -> Builder
atomBits_OUCH m a = let
  n = fromIntegral $ numberBondsToHeavyAtomsAtIndex m $ fromJust $ getIndexForAtom a
  i = fromInteger $ atomicNumber a
  in case a of
    Element {} ->  B.putWord64le ((bit (mod i 64)) :: Word64)
         `B.append` B.singleton ((bit n) :: Word8)
    Open {}    ->  B.putWord64le ((bit 63) :: Word64)
         `B.append` B.singleton ((bit n) :: Word8)
    _          -> B.putWord64le (0 :: Word64)
         `B.append` B.singleton (0 :: Word8)
{------------------------------------------------------------------------------}
atomBits_RECURSIVE :: Int
                   -> Molecule
                   -> Atom
                   -> Builder
atomBits_RECURSIVE depth  m a =  let
  i = fromJust $ getIndexForAtom a
  f_pb = pathBits atomBits_OUCH bondBits_OUCH
  paths = findPaths depth (PGraph m U.empty U.empty) i
  in V.foldr (\p b -> f_pb p .||. b) B.empty paths



{------------------------------------------------------------------------------}
bondBits_OUCH :: Molecule
              -> Bond
              -> Builder
bondBits_OUCH m b = B.singleton (bit $ bondKey b::Word8)


{------------------------------------------------------------------------------}
molBits_N :: Int
          -> Molecule
          -> Builder
molBits_N depth m = B.append (sizeBits_OUCH m)
                             (moleculeBits atomBits_OUCH bondBits_OUCH depth m)


{------------------------------------------------------------------------------}
molBits_ID :: Int
           -> Molecule
           -> Builder
molBits_ID depth m = B.append (sizeBits_OUCH m)
                              (moleculeBits atomBits_R bondBits_OUCH depth m)
  where atomBits_R = atomBits_RECURSIVE depth



{------------------------------------------------------------------------------}
molBits_OUCH :: Molecule
             -> Builder
molBits_OUCH m = B.append (sizeBits_OUCH m)
                          (moleculeBits atomBits_OUCH bondBits_OUCH 7 m)


{------------------------------------------------------------------------------}
sizeBits_OUCH :: Molecule
              -> Builder
sizeBits_OUCH m = B.singleton (n::Word8)
  where n = fromIntegral $ Map.size $ atomMap m


{------------------------------------------------------------------------------}
moleculeBits :: (Molecule -> Atom -> Builder)
             -> (Molecule -> Bond -> Builder)
             -> Int
             -> Molecule
             -> Builder
moleculeBits atomB bondB depth m = let
  pathBits_OUCH = pathBits atomB bondB
  paths = allPaths depth m
  in V.foldr (\p b -> b .||. pathBits_OUCH p) B.empty paths


{------------------------------------------------------------------------------}
pathBits :: (Molecule -> Atom -> Builder)
         -> (Molecule -> Bond -> Builder)
         -> PGraph
         -> Builder
pathBits atomB bondB p@PGraph {molecule=m, vertexList=v} = let
  x = U.head v
  xs = U.tail v
  atomBits = atomB m atom
  atom = fromJust $ getAtomAtIndex m x
  bond = Set.findMax $ Set.filter (\b -> (bondsTo b) == U.head xs)
                                  $ atomBondSet atom
  bondBits = bondB m bond
  bits = B.append atomBits bondBits
  output | U.length v == 0 = B.empty
         | U.length v == 1 = atomBits
         | otherwise = B.append bits $ pathBits atomB bondB p {vertexList=xs}
  in output

-- Logical OR where bytes expand to the length of the longest pair
{------------------------------------------------------------------------------}
(.||.) :: Builder
       -> Builder
       -> Builder
(.||.) b1 b2 = let
  bytes1 = L.unpack $ B.toLazyByteString b1
  bytes2 = L.unpack $ B.toLazyByteString b2
  l1 = List.length bytes1
  l2 = List.length bytes2
  bytes1' | l1 > l2 = bytes1
          | otherwise = bytes1 ++ (List.replicate (l2 - l1) (0::Word8))
  bytes2' | l2 > l1 = bytes2
          | otherwise = bytes2 ++ (List.replicate (l1 - l2) (0::Word8))
  zipped = List.zip bytes1' bytes2'
  logicalOrList = List.map (\(a1, a2) -> a1 .|. a2) zipped
  in  B.fromLazyByteString $ L.pack logicalOrList

-- Logical OR where list contracts to the length of the shortest pair
{------------------------------------------------------------------------------}
(.|||.) :: Builder
        -> Builder
        -> Builder
(.|||.) b1 b2 = let
  bytes1 = L.unpack $ B.toLazyByteString b1
  bytes2 = L.unpack $ B.toLazyByteString b2
  l1 = List.length bytes1
  l2 = List.length bytes2
  zipped = List.zip bytes1 bytes2
  logicalOrList = List.map (\(a1, a2) -> a1 .|. a2) zipped
  in  B.fromLazyByteString $ L.pack logicalOrList

-- Logical OR where list expands to the length of the RIGHT argument
(.||>.) :: Builder
        -> Builder
        -> Builder
(.||>.) b1 b2 = let
  bytes1 = L.unpack $ B.toLazyByteString b1
  bytes2 = L.unpack $ B.toLazyByteString b2
  l1 = List.length bytes1
  l2 = List.length bytes2
  bytes1' | l1 > l2 = bytes1
          | otherwise = bytes1 ++ (List.replicate (l2 - l1) (0::Word8))
  zipped = List.zip bytes1' bytes2
  logicalOrList = List.map (\(a1, a2) -> a1 .|. a2) zipped
  in  B.fromLazyByteString $ L.pack logicalOrList

{------------------------------------------------------------------------------}
{-- allPaths --}
-- Returns all paths up to a given depth.  Always an even numbered list.
allPaths :: Int
         -> Molecule
         -> (V.Vector PGraph)
allPaths depth m = V.foldr (\i p -> p `seq` p V.++ findPaths depth (PGraph m U.empty U.empty) i)
                              V.empty indexList
  where indexList = V.fromList $ Map.keys $ atomMap m

-- | Returns all paths up to a depth that start at an external position of the
-- molecule.
allTerminalPaths :: Int
                 -> Molecule
                 -> V.Vector PGraph
allTerminalPaths depth m = V.foldr (\i p -> p `seq`
                                       p V.++ findPaths depth
                                       (PGraph m U.empty U.empty) i)
                                       V.empty (maxAtomicNumber m $ allTerminalVertices m)


-- | Returns a list of all terminal atoms in a molecule (those that have the
-- LEAST number of vertices
allTerminalVertices :: Molecule
                    -> V.Vector Int
allTerminalVertices m = V.fromList $ List.foldr (\a acc -> if ((snd a) == minValence)
                                              then (fst a):acc
                                              else acc) [] zippedIndex
  where indexList = Map.keys $ atomMap m
        valenceList = List.map (\m_i -> Set.size $ atomBondSet
                                                 $ fromJust
                                                 $ getAtomAtIndex m m_i)
                                                 indexList
        minValence = List.minimum $ valenceList
        zippedIndex = List.zip indexList valenceList

maxAtomicNumber :: Molecule
                -> V.Vector Int
                -> V.Vector Int
maxAtomicNumber m m_is = let
  atomicNums = V.map (\i -> atomicNumber $ fromJust $ getAtomAtIndex m i) m_is
  maxAN = V.maximum atomicNums
  zippedIndex = V.zip m_is atomicNums
  in V.fromList $ V.foldr (\a acc -> if ((snd a) == maxAN)
                                     then (fst a):acc
                                     else acc) [] zippedIndex

{------------------------------------------------------------------------------}
{-- longestPaths --}
-- Returns the longest paths found in a molecule.  Because paths can go
-- in either direction, this will always be an even numbered list.
longestPaths :: Molecule
             -> V.Vector PGraph
longestPaths m = let
  depth = Map.size (atomMap m)
  paths = allPaths depth m
  -- paths = allTerminalPaths depth m
  maxLength = V.maximum $ V.map pathLength paths
  longest = V.filter ((==maxLength) . pathLength) paths
  in longest

{------------------------------------------------------------------------------}
-- | Takes a molecule and a starting position that is connected to the
-- first PGraph and finds the longest chain of connections in the molecule
-- choosing only atoms that are NOT in the first PGraph.  If more than
-- one exists, returns the LONGEST LEAST of these.
longestLeastAnchoredPath :: PGraph
                         -> Int
                         -> PGraph
longestLeastAnchoredPath exclude@PGraph{molecule=m, vertexList=l, root=r} anchor = let
  depth = fromIntegral $ pathLength exclude
  paths = findPathsExcluding (vToSet $ r U.++ l) depth (exclude {vertexList=U.empty}) anchor
  nonExcludedPaths = V.filter (\a -> False == hasOverlap exclude a) paths
  rootedPaths = V.map (\p -> p {root=(U.singleton anchor)}) nonExcludedPaths
  output | V.length nonExcludedPaths == 0 = exclude {vertexList=U.empty}
         | otherwise = findLongestLeastPath nonExcludedPaths 0
  in output

{------------------------------------------------------------------------------}
-- | Takes a list of paths of the same length and from the same molecule and
-- returns the "least" path according to atom ordering rules.  Used in selecting a
-- path for canonicalization.
-- Not used directly.
findLongestLeastPath :: V.Vector PGraph   -- ^ The paths to select from
                     -> Int        -- ^ The current index being compared
                     -> PGraph     -- ^ The longest least path from starting list
--findLongestLeastPath gs i | (trace $ show gs) False = undefined
--findLongestLeastPath gs i | (trace $ "#Paths: " ++ (show $ List.length gs)) False = undefined
--findLongestLeastPath V.empty i = PGraph emptyMolecule U.empty U.empty
findLongestLeastPath gs i = let
  gsL = pLongest gs
  mol = molecule (V.head gs)
  ranks r acc | acc==LT || r  ==LT = LT
              | acc==EQ && r  ==GT = GT
              | acc==EQ && r  ==EQ = EQ
              | acc==GT            = GT
  foldRanks g = V.foldl (\acc a -> ranks (ordAtom g a i) acc ) EQ gsL
  mapRanks = V.map (\a -> foldRanks a) gsL
  topRank = V.maximum mapRanks
  gs' = V.filter ((==topRank) . foldRanks) gsL
  output | V.length gs == 0 = PGraph emptyMolecule U.empty U.empty
         | V.length gsL == 0 =  PGraph emptyMolecule U.empty U.empty
         | V.length gsL == 1 = V.head gsL
         | pathLength (V.head gsL) == (fromIntegral i) = V.head gsL
         | otherwise = gs' `seq` findLongestLeastPath gs' (i+1)
  in  output



comparePaths :: PGraph
             -> PGraph
             -> Ordering
comparePaths p1 p2 = let
  mol = molecule p1
  ranks r acc | acc==LT            = LT
              | acc==EQ && (r==LT) = LT
              | acc==EQ && (r==GT) = GT
              | acc==EQ && (r==EQ) = EQ
              | acc==GT            = GT
  rankMap = List.map (\i -> ordAtom p1 p2 i) [0..fromInteger ((pathLength p1) -1 )]
  output | (pathLength p1) > (pathLength p2)  = GT
         | (pathLength p1) < (pathLength p2)  = LT
         | (pathLength p1) == (pathLength p2) = List.foldl (\acc a -> ranks a acc) EQ rankMap
  in output


{------------------------------------------------------------------------------}
-- | Finds the longest least path in a molecule.  Used for canonicalization.
longestLeastPath :: Molecule  -- ^ The Molecule
                 -> PGraph    -- ^ The longest least path
longestLeastPath m = let
  paths = longestPaths m
  in paths `seq` findLongestLeastPath paths 0



ordByRootBond :: PGraph  -- ^ The first path to compare and its root index
              -> PGraph   -- ^ The second path to compare and its root index
              -> Int      -- ^ The index to comparea
              -> Ordering -- ^ The Ord result
ordByRootBond p1@PGraph {root=r1}
              p2@PGraph {root=r2} index = let
  mol = molecule p1
  output | index /= 0 = EQ
         | U.length r1 == 0 || U.length r2 == 0 = EQ
         | otherwise = compare (bondBetweenIndices mol index (U.head r1))
                               (bondBetweenIndices   mol index (U.head r2))
  in output



ordByOffPathBond :: PGraph  -- ^ The first path to compare and its root index
                 -> PGraph   -- ^ The second path to compare and its root index
                 -> Int      -- ^ The index to comparea
                 -> Ordering -- ^ The Ord result
ordByOffPathBond p1 p2 index = ordBondList (offPathBondTypeList p1 index)
                                           (offPathBondTypeList p2 index)


ordByBond :: PGraph  -- ^ The first path to compare and its root index
          -> PGraph   -- ^ The second path to compare and its root index
          -> Int      -- ^ The index to comparea
          -> Ordering -- ^ The Ord result
ordByBond p1 p2 index = ordBondList (bondTypeList p1 index)
                                    (bondTypeList p2 index)


ordBondList :: [NewBond] -> [NewBond] -> Ordering
ordBondList [] [] = EQ
ordBondList [] p = LT
ordBondList p [] = GT
ordBondList (bl1:bls1) (bl2:bls2) = let
  output | compare bl1 bl2 /= EQ = compare bl1 bl2
         | otherwise = ordBondList bls1 bls2
  in output


offPathBondTypeList :: PGraph -> Int -> [NewBond]
offPathBondTypeList path p_i = let
  index = pathIndex path p_i
  mol = molecule path
  bondTypeList = List.map (\a -> bondBetweenIndices mol index a) offPathList
  offPathList =  Set.toList $ Set.difference (bondIndexSet path p_i)
                                             (vToSet $ vertexList path)
  in List.sort bondTypeList

bondTypeList :: PGraph -> Int -> [NewBond]
bondTypeList path p_i = let
  index = pathIndex path p_i
  mol = molecule path
  bondTypeList = List.map (\a -> bondBetweenIndices mol index (bondsTo a))
                         $ Set.toList
                         $ atomBondSet
                         $ fromJust
                         $ getAtomAtIndex mol index
  in List.sort bondTypeList

-- | A comparison utility for ordAtom (below) that does the recursive path
-- comparison at a given index position
ordByPath :: PGraph  -- ^ The first path to compare and its root index
          -> PGraph   -- ^ The second path to compare and its root index
          -> Int      -- ^ The index to comparea
          -> Ordering -- ^ The Ord result
ordByPath p1 p2 index  = ordPathList (branchPaths p1 index)
                                     (branchPaths p2 index)

ordByNextBond :: PGraph  -- ^ The first path to compare and its root index
              -> PGraph   -- ^ The second path to compare and its root index
              -> Int      -- ^ The index to comparea
              -> Ordering -- ^ The Ord result
ordByNextBond p1 p2 index = let
  output | (index + 1) == (fromIntegral $ pathLength p1) = EQ
         | (index + 1) == (fromIntegral $ pathLength p2) = EQ
         | otherwise = compare (bondBetweenIndices mol m_i1 m_i1')
                               (bondBetweenIndices mol m_i2 m_i2')
  mol = molecule p1
  m_i1 = pathIndex p1 index
  m_i1' = pathIndex p1 (index + 1)
  m_i2 = pathIndex p2 index
  m_i2' = pathIndex p2 (index + 1)
  in output

bondIndexSet p p_i = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust $
                                 getAtomAtIndex (molecule p) (pathIndex p p_i)
pathIndexSet p = vToSet $ (root p) U.++ (vertexList p)
validIndexList p p_i = Set.toList $ Set.difference (bondIndexSet p p_i) (pathIndexSet p )
pLongest ps = V.filter (\p -> longest == pathLength p) ps
                  where longest = V.maximum $ (V.map pathLength ps)
branchPaths p p_i = V.fromList $ List.map (\a -> longestLeastAnchoredPath pNew a) (validIndexList p p_i)
  where pNew = p {root=(vertexList p)}




ordPathList :: V.Vector PGraph -> V.Vector PGraph -> Ordering
--ordPathList p1 p2 | (trace $ "OrdPath Compare-- A: " ++ (show p1) ++ (show p2))  False = undefined
ordPathList ps1 ps2 = let
  ll1 = findLongestLeastPath ps1 0
  ll2 = findLongestLeastPath ps2 0
  xp1 = V.filter (/=ll1) ps1
  xp2 = V.filter (/=ll2) ps2
  comp = comparePaths ll1 ll2
  output | V.length ps1 == 0 && V.length ps2 == 0 = EQ
         | V.length ps1 == 0 = LT
         | V.length ps2 == 0 = GT
         | comp /= EQ = comp
         | otherwise = ordPathList xp1 xp2
  in output


pathIndex :: PGraph -> Int -> Int
pathIndex p p_i  = (vertexList p) U.! p_i
{-# INLINE pathIndex #-}


{------------------------------------------------------------------------------}
-- | Orders atoms in a path to aid in path selection.  An ugly utility to affect
-- Ord for paths.  This function assumes the molecules in both paths are the same
-- It cannot check for this because such a check implcitly relies on this function.
-- Naturally, this is never used directly.
ordAtom :: PGraph   -- ^ The first path to compare
        -> PGraph   -- ^ The second path to compare
        -> Int      -- ^ The index to compare
        -> Ordering -- ^ The Ord result
ordAtom p1 p2 p_i = let
  m = molecule p1
  i1 = pathIndex p1 p_i
  i2 = pathIndex p2 p_i
  atom1 = fromJust $ getAtomAtIndex m i1
  atom2 = fromJust $ getAtomAtIndex m i2
  atoms = Map.size $ atomMap m
  byNumber = compare (atomicNumber atom1) (atomicNumber atom2)
  byIsotope = compare (neutronNumber atom1) (neutronNumber atom2)
  byVertex = compare (Set.size $ atomBondSet atom1) (Set.size $ atomBondSet atom2)
  byBranchBond = ordByBond p1 p2 p_i
  byRootBond = ordByRootBond p1 p2 p_i
  byOffPathBond = ordByOffPathBond p1 p2 p_i
  byNextBond = ordByNextBond p1 p2 p_i
  byPath = ordByPath p1 p2 p_i
  fAtom d = compared
    where test = compare (B.toLazyByteString $ atomBits_RECURSIVE d m atom1)
                         (B.toLazyByteString $ atomBits_RECURSIVE d m atom2)
          compared | d >= (Map.size $ atomMap m) = test
                   | test == EQ = fAtom (d+1)
                   | otherwise = test
  output = case atom1 of
            Element {} -> case atom2 of
              Element {} -> ordElements
              _  -> GT
            _  -> case atom2 of
              Element {} -> LT
              _ -> compare i1 i2

  ordElements

           | byNumber    /= EQ =  byNumber
           | byIsotope   /= EQ = byIsotope
           | byNextBond /= EQ = byNextBond
           | byRootBond /= EQ = byRootBond
           | byOffPathBond /= EQ  = byOffPathBond
           | byBranchBond    /= EQ =  byBranchBond
           | byVertex    /= EQ = byVertex
           | i1 == i2 = EQ
           | byPath      /= EQ = byPath

         -- If all of the above are EQ, then the atoms REALLY ARE chemically equivalent
         -- and cannot be distinguished.
          | otherwise         = EQ
  in output


{------------------------------------------------------------------------------}
-- | Find all paths starting from a given index, but excluding traversal through
-- the indices in the given exclusion set.
-- This is a utility function, not used directly.
findPathsExcluding :: Set Int  -- ^ The atom index set to exclude from paths
                   -> Int      -- ^ The maximum depth to recursively add to path
                   -> PGraph   -- ^ The path we are recursively adding to
                   -> Int      -- ^ The atom index to add to the growing path
                   -> V.Vector PGraph -- ^ The new paths created after terminal recursion
findPathsExcluding exclude depth path@PGraph {molecule=m, vertexList=l} index = let
  path' = path {vertexList=(l U.++ U.singleton index)}
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex m index
  pathIndexSet = Set.union exclude $ vToSet l
  validIndexSet = Set.difference bondIndexSet pathIndexSet
  accPath i p = p `seq` p V.++ (findPathsExcluding exclude depth path' i)
  paths | Set.size validIndexSet == 0      = V.singleton path'
        | U.length l > depth               = V.singleton path'
        | otherwise = Set.fold accPath V.empty validIndexSet
  in paths


{------------------------------------------------------------------------------}
-- | Find all possible paths starting from an atom index up to specified depth.
-- This is a utility function, not used directly.
findPaths :: Int      -- ^ The maximum depth
          -> PGraph   -- ^ The path we are building up
          -> Int      -- ^ The atom indices to add to the growing path
          -> V.Vector PGraph -- ^ The new paths created after terminal recursion
findPaths depth path@PGraph {molecule=m, vertexList=l} index = let
  path' = path {vertexList=(l U.++ U.singleton index)}
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust
                                           $ getAtomAtIndex m index
  pathIndexSet = vToSet l
  validIndexSet = Set.difference bondIndexSet pathIndexSet
  accPath i p = p `seq` p V.++ (findPaths depth path' i)
  paths | Set.size validIndexSet == 0   = V.singleton path'
        | U.length l > depth            = V.singleton path'
        | otherwise = Set.fold accPath V.empty validIndexSet
  in paths


{------------------------------------------------------------------------------}


