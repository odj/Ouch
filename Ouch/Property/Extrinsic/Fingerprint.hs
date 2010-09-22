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
    atomBits_OUCH
  , bondBits_OUCH
  , molBits_OUCH
  , molBits_N
  , molBits_ID
  , pathBits
  , findPaths
  , allPaths
  , longestPaths
  , longestLeastPath
  , longestLeastAnchoredPath
  , writeCanonicalPath
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
import Ouch.Structure.Molecule
import Ouch.Property.Ring
import Ouch.Data.Atom
import Ouch.Data.Bond
import Ouch.Input.Smiles
import Data.Maybe
import Data.Set as Set
import Data.List as List
import Data.Map as Map

data Fingerprint = SimpleFingerPrint

{------------------------------------------------------------------------------}
atomBits_OUCH :: Molecule -> Atom -> Builder
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
atomBits_RECURSIVE :: Int -> Molecule -> Atom -> Builder
atomBits_RECURSIVE depth  m a =  let
  i = fromJust $ getIndexForAtom a
  f_pb = pathBits atomBits_OUCH bondBits_OUCH
  paths = findPaths depth (PGraph m []) i
  in List.foldr (\p b -> f_pb p .||. b) B.empty paths



{------------------------------------------------------------------------------}
bondBits_OUCH :: Molecule -> Bond -> Builder
bondBits_OUCH m b = B.singleton (bit $ bondKey b::Word8)


{------------------------------------------------------------------------------}
molBits_N :: Int -> Molecule -> Builder
molBits_N depth m = B.append (sizeBits_OUCH m) (moleculeBits atomBits_OUCH bondBits_OUCH depth m)


{------------------------------------------------------------------------------}
molBits_ID :: Int -> Molecule -> Builder
molBits_ID depth m = B.append (sizeBits_OUCH m) (moleculeBits atomBits_R bondBits_OUCH depth m)
  where atomBits_R = atomBits_RECURSIVE depth



{------------------------------------------------------------------------------}
molBits_OUCH :: Molecule -> Builder
molBits_OUCH m = B.append (sizeBits_OUCH m) (moleculeBits atomBits_OUCH bondBits_OUCH 7 m)


{------------------------------------------------------------------------------}
sizeBits_OUCH :: Molecule -> Builder
sizeBits_OUCH m = B.singleton (n::Word8)
  where n = fromIntegral $ Map.size $ atomMap m


{------------------------------------------------------------------------------}
moleculeBits :: (Molecule -> Atom -> Builder) -> (Molecule -> Bond -> Builder) -> Int -> Molecule -> Builder
moleculeBits atomB bondB depth m = let
  pathBits_OUCH = pathBits atomB bondB
  paths = allPaths depth m
  in List.foldr (\p b -> b .||. pathBits_OUCH p) B.empty paths


{------------------------------------------------------------------------------}
pathBits :: (Molecule -> Atom -> Builder) -> (Molecule -> Bond -> Builder) -> PGraph -> Builder
pathBits atomB bondB p@PGraph {molecule=m, vertexList=[]} = B.empty
pathBits atomB bondB p@PGraph {molecule=m, vertexList=x:[]} = let
  atom = fromJust $ getAtomAtIndex m x
  atomBits = atomB m atom
  in atomBits
pathBits atomB bondB p@PGraph {molecule=m, vertexList=x:xs} = let
  atom = fromJust $ getAtomAtIndex m x
  bond = Set.findMax $ Set.filter (\b -> (bondsTo b) == List.head xs) $ atomBondSet atom
  atomBits = atomB m atom
  bondBits = bondB m bond
  bits = B.append atomBits bondBits
  in B.append bits $ pathBits atomB bondB p {vertexList=xs}


-- Logical OR where bytes expand to the length of the longest pair
{------------------------------------------------------------------------------}
(.||.) :: Builder -> Builder -> Builder
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
(.|||.) :: Builder -> Builder -> Builder
(.|||.) b1 b2 = let
  bytes1 = L.unpack $ B.toLazyByteString b1
  bytes2 = L.unpack $ B.toLazyByteString b2
  l1 = List.length bytes1
  l2 = List.length bytes2
  zipped = List.zip bytes1 bytes2
  logicalOrList = List.map (\(a1, a2) -> a1 .|. a2) zipped
  in  B.fromLazyByteString $ L.pack logicalOrList

-- Logical OR where list expands to the length of the RIGHT argument
(.||>.) :: Builder -> Builder -> Builder
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
allPaths :: Int -> Molecule -> [PGraph]
allPaths depth m = List.foldr (\i p -> p ++ findPaths depth (PGraph m []) i) [] indexList
  where indexList = Map.keys $ atomMap m

{------------------------------------------------------------------------------}
{-- longestPaths --}
-- Returns the longest paths found in a molecule.  Because paths can go
-- in either direction, this will always be an even numbered list.
longestPaths :: Molecule -> [PGraph]
longestPaths m = let
  depth = Map.size (atomMap m)
  paths = allPaths depth m
  maxLength = List.maximum $ List.map pathLength paths
  longest = List.filter ((==maxLength) . pathLength) paths
  in longest

{------------------------------------------------------------------------------}
-- findLongestLeastAnchoredPath
-- This takes a molecule and a starting position that is connected to the
-- first PGraph and finds the longest chain of connections in the molecule
-- choosing only atoms that are NOT in the first PGraph.  If more than
-- one exists, returns the LONGEST LEAST of these.
longestLeastAnchoredPath :: PGraph -> Int -> PGraph
longestLeastAnchoredPath exclude@PGraph{vertexList=l} anchor = let
  depth = fromIntegral $ pathLength exclude
  paths = findPathsExcluding (Set.fromList l) depth (exclude {vertexList=[]}) anchor
  nonExcludedPaths = List.filter (\a -> False == hasOverlap exclude a) paths
  output | List.length nonExcludedPaths == 0 = exclude {vertexList=[]}
         | otherwise = findLongestLeastPath nonExcludedPaths 0
  in output

{-- findLongestLeastPath --}
-- Takes a list of paths of the same length and from the same molecule and
-- returns the "least" path according to atom ordering rules.  Used in selecting a
-- path for canonicalization.
findLongestLeastPath :: [PGraph] -> Int -> PGraph
findLongestLeastPath [] i = PGraph emptyMolecule []
findLongestLeastPath gs i = let
  mol = molecule (gs!!0)
  ranks r acc | r  ==GT || acc==GT = GT
              | acc==EQ            = EQ
              | r  ==EQ            = EQ
              | r  ==LT            = LT
  foldRanks g = List.foldl (\acc a -> ranks (ordAtom mol ((vertexList g)!!i)
                                                         ((vertexList a)!!i)) acc ) LT gs
  mapRanks = List.map (\a -> foldRanks a) gs
  leastRank = List.minimum mapRanks
  gs' = List.filter ((==leastRank) . foldRanks) gs
  output | List.length gs == 1     = gs!!0
         | pathLength (gs!!0) == (fromIntegral i) = gs!!0
         | otherwise = findLongestLeastPath gs' (i+1)
  in output

{------------------------------------------------------------------------------}
{-- longestLeastPath --}
-- Finds the longest least path in a molecule.  Used for canonicalization.
longestLeastPath :: Molecule -> PGraph
longestLeastPath m = let
  paths = longestPaths m
  in findLongestLeastPath paths 0


{------------------------------------------------------------------------------}
{-- ordAtom --}
-- Orders atoms in a path to aid in path selection.
ordAtom :: Molecule -> Int -> Int -> Ordering
ordAtom m i1 i2 = let
  atom1 = fromJust $ getAtomAtIndex m i1
  atom2 = fromJust $ getAtomAtIndex m i2
  atoms = Map.size $ atomMap m
  byNumber = compare (atomicNumber atom1) (atomicNumber atom2)
  byIsotope = compare (neutronNumber atom1) (neutronNumber atom2)
  byVertex = compare (Set.size $ atomBondSet atom1) (Set.size $ atomBondSet atom2)
  byFinger = fAtom atoms
  fAtom i = compared
    where test = compare (B.toLazyByteString $ atomBits_RECURSIVE i m atom1)
                         (B.toLazyByteString $ atomBits_RECURSIVE i m atom2)
          compared | i >= (Map.size $ atomMap m) = test
                   | test == EQ = fAtom (i+1)
                   | otherwise = test
  output
         -- The atoms are the same index
         | i1 == i2          = EQ

         -- The atoms are the same element
         | byNumber    /= EQ = byNumber

         -- Atoms are the same isotope
         | byIsotope   /= EQ = byIsotope

         -- Atoms have the same number of connections
         | byVertex    /= EQ = byVertex

         -- Atoms have the same recursive fingerprint
         | byFinger    /= EQ = byFinger

         -- If all of the above are EQ, then the atoms REALLY ARE chemically equivalent
         -- and cannot be distinguished.
         | otherwise         = EQ
  in output


{------------------------------------------------------------------------------}
{-- findPathsExcluding --}
-- Find all paths starting from a given index, but excluding traversal through
-- the indices in the given exclusion set.
findPathsExcluding :: Set Int -> Int ->  PGraph -> Int -> [PGraph]
findPathsExcluding exclude depth path@PGraph {molecule=m, vertexList=l} index  = let
  path' = path {vertexList=(l ++ [index])}
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust $ getAtomAtIndex m index
  pathIndexSet = Set.union exclude $ Set.fromList l
  validIndexSet = Set.difference bondIndexSet pathIndexSet
  accPath i p = p ++ (findPathsExcluding exclude depth path' i)
  paths | Set.size validIndexSet == 0      = [path']
        | List.length l > depth            = [path']
        | otherwise = Set.fold accPath [] validIndexSet
  in paths

{------------------------------------------------------------------------------}
findPaths :: Int ->  PGraph -> Int -> [PGraph]
findPaths depth path@PGraph {molecule=m, vertexList=l} index  = let
  path' = path {vertexList=(l ++ [index])}
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust $ getAtomAtIndex m index
  pathIndexSet = Set.fromList l
  validIndexSet = Set.difference bondIndexSet pathIndexSet
  accPath i p = p ++ (findPaths depth path' i)
  paths | Set.size validIndexSet == 0      = [path']
        | List.length l > depth            = [path']
        | otherwise = Set.fold accPath [] validIndexSet
  in paths


{------------------------------------------------------------------------------}
{-- writeCanonicalPath --}
-- Writes the SMILES string for a given molecule
writeCanonicalPath :: Molecule -> String
writeCanonicalPath m = let
  backbone = longestLeastPath m
  in writePath [] backbone 0 False

{------------------------------------------------------------------------------}
{-- writePath --}
-- Writes the SMILES string for a given path
writePath :: [PGraph] -> PGraph -> Int -> Bool -> String
writePath gx g i subStructure = let
  mol = molecule g
  s = writeStep g i
  endOfPath = i == (fromInteger $ pathLength g)
  output | endOfPath && subStructure = ")"
         | endOfPath = ""
         | otherwise =  s ++ writeSubpath gx g i
                          ++ writePath gx g (i+1) subStructure
  in output

{------------------------------------------------------------------------------}
{-- writeSubpath --}
-- Writes the SMILES strings for all subpaths (if any exist) at position i in a
-- given path g, excluding travesal through any atoms in the paths gx
writeSubpath :: [PGraph] -> PGraph -> Int -> String
writeSubpath gx g i = let
  mol = molecule g
  vertices = vertexList g
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust $ getAtomAtIndex mol (vertices!!i)
  pathIndexSet = List.foldr (\a acc -> Set.union acc $ Set.fromList $ vertexList a) Set.empty (g:gx)
  validIndexList = Set.toList $ Set.difference bondIndexSet pathIndexSet
  pathIndexList = Set.toList pathIndexSet
  branchPaths = List.map (\a -> longestLeastAnchoredPath g {vertexList=pathIndexList} a) validIndexList
  nextBranch = findLongestLeastPath branchPaths 0
  s = writeStep g i
  endOfPath = i == (fromInteger $ pathLength g)
  output | (pathLength nextBranch) > 0 =  "(" ++ writePath (g:gx) nextBranch 0  True
                                              ++ writeSubpath (nextBranch:gx) g i
         | otherwise = ""
  in output

{------------------------------------------------------------------------------}
{-- writeStep --}
-- Writes atom and bond information from position i in a path
writeStep :: PGraph -> Int -> String
writeStep g i = let
  mol = molecule g
  atom = fromJust $ getAtomAtIndex mol ((vertexList g)!!i)
  in atomicSymbolForAtom atom
















