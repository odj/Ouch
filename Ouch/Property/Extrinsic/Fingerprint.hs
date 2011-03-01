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
import Ouch.Property.Graph
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
  paths = findPaths depth [] (PGraph m U.empty U.empty) i
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


