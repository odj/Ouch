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
allPaths :: Int -> Molecule -> [PGraph]
allPaths depth m = List.foldr (\i p -> p ++ findPaths depth (PGraph m []) i) [] indexList
  where indexList = Map.keys $ atomMap m

{------------------------------------------------------------------------------}
findPaths :: Int ->  PGraph -> Int -> [PGraph]
findPaths depth path@PGraph {molecule=m, vertexList=l} index  = let
  path' = path {vertexList=index:l}
  bondIndexSet = Set.map (\a -> bondsTo a) $ atomBondSet $ fromJust $ getAtomAtIndex m index
  pathIndexSet = Set.fromList l
  validIndexSet = Set.difference bondIndexSet pathIndexSet
  accPath i p = p ++ (findPaths depth path' i)
  paths | Set.size validIndexSet == 0      = [path']
        | List.length l > depth           = [path']
        | otherwise = Set.fold accPath [] validIndexSet
  in paths



















