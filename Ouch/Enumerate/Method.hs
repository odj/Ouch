{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Enumerate.Method
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


-- | Provides a framework for enumerating structures according to recursively
-- constructed 'Method's
module Ouch.Enumerate.Method (
    Method(..)
  , (>#>)
  , (>&&>)
  , addMethod
  , classSelector
  , openValenceSelector
  , halogens
  , alkyls
  , fingerprintFilterBuilder
  , addH
  , removeH
  , stripMol
  , removeE
  , elementSelector
  , addMols
  , makeUnique
  ) where

import Ouch.Structure.Atom
import {-# SOURCE #-} Ouch.Structure.Molecule
import Ouch.Structure.Marker
import {-# SOURCE #-} Ouch.Input.Smiles
import Ouch.Output.SDF
import Ouch.Output.Mol

import Data.Map as Map
import Data.List as List
import Data.Set as Set
import Data.Maybe as Maybe
import Control.Applicative
import {-# SOURCE #-} Ouch.Property.Extrinsic.Fingerprint
import Data.ByteString.Lazy as L
import Data.Binary.Builder as B


{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}



data Method   = NoMethod     { firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             }

              | AddMethod    { firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , addList      ::  ([(NewBond, Molecule)])
                             }

              | InsertMethod { firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , insertRule   ::  (Molecule -> [Int])
                             , insertList   ::  [(Molecule, [Int])]
                             }

              | ReplaceMethod { firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , replaceList  ::  [Molecule]
                             }

              | ReactMethod  { firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , reactList    ::  [(Molecule -> [Molecule])]
                             }

              | FilterMethod { firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , molFilter    ::  ([Molecule] -> [Molecule])
                             }




{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}


(>#>) :: [Molecule] -> (Maybe Method) -> [Molecule]
(>#>) ms mMethod = case mMethod of
  Nothing  -> ms
  Just method -> case method of
    NoMethod      {} -> ms >#> (firstApply method) >#> (lastApply method)
    AddMethod     {} -> addMethod ms method
    InsertMethod  {} -> ms >#> (firstApply method)
    ReplaceMethod {} -> replaceMethod ms method
    ReactMethod   {} -> ms >#> (firstApply method)
    FilterMethod  {} -> filterMethod ms method

(>##>) :: [Molecule] -> (Maybe Method) -> [Molecule]
(>##>) ms mMethod = ms ++ (ms >#> mMethod)


addMethod :: [Molecule] -> Method -> [Molecule]
addMethod ms method = let
  mols = ms >#> (firstApply method)
  atomList f m = Map.keys $ fst $ Map.partition (f m) (atomMap m)
  atomLists = List.map (atomList $ selector method) mols
  zipped = List.zip mols atomLists
  makeMol m1 i1 addItem = [connectMoleculesAtIndicesWithBond m1 i1 (snd addItem) 0 (fst addItem)]
  newMols m i = ((addList method) >>= (makeMol m i))
  molMapper m l = List.concat $ List.map (\i -> newMols m i) l
  output  = List.concat
          $ List.map (\(m, l) -> molMapper m l ) zipped
  in  output >#> (lastApply method)


replaceMethod :: [Molecule] -> Method -> [Molecule]
replaceMethod ms method = let
  mols = ms >#> (firstApply method)
  atomList f m = Map.keys $ fst $ Map.partition (f m) (atomMap m)
  atomLists = List.map (atomList $ selector method) mols
  zipped = List.zip mols atomLists
  output = undefined
  in ms

filterMethod :: [Molecule] -> Method -> [Molecule]
filterMethod ms method = let
   mols = ms >#> (firstApply method)
   output = (molFilter method ) mols
   in output >#> (lastApply method)

{------------------------------------------------------------------------------}
{-------------------------------Convenience Functions--------------------------}
{------------------------------------------------------------------------------}

fingerprintFilterBuilder :: (Ord a) => (Molecule -> a) -> ([Molecule] -> [Molecule])
fingerprintFilterBuilder fp = let
  newFilter ms = output
    where fingerprints = List.map (\m -> fp m `seq` fp m) ms
          filterMap = Map.fromList $ List.zip fingerprints ms
          output = Map.elems filterMap
  in newFilter

elementSelector :: String -> (Molecule -> Atom -> Bool)
elementSelector s = (\m a -> isElementType s a)

-- A convenience function to create a selector based on class number
classSelector :: Integer -> (Molecule -> Atom -> Bool)
classSelector i = mk
  where mk _ atom = case getMarker atom (Class 0) of
          Nothing -> False
          Just mk -> i == classNumber mk

openValenceSelector :: Molecule -> Atom -> Bool
openValenceSelector m a | freeValence > 0 = True | otherwise = False
  where freeValence = case getIndexForAtom a of
          Just i  -> freeValenceAtIndex m i
          Nothing -> 0

-- Convenience to create a new selector from two with AND logic
(>&&>) :: (Molecule -> Atom -> Bool) -> (Molecule -> Atom -> Bool) -> (Molecule -> Atom -> Bool)
(>&&>) sel1 sel2 = (\m a -> (sel1 m a) && (sel2 m a))

-- Convenience to create a new selector from two with OR logic
(>||>) :: (Molecule -> Atom -> Bool) -> (Molecule -> Atom -> Bool) -> (Molecule -> Atom -> Bool)
(>||>) sel1 sel2 = (\m a -> (sel1 m a) || (sel2 m a))

{------------------------------------------------------------------------------}
{-------------------------------List Generators--------------------------------}
{------------------------------------------------------------------------------}

halogens :: [(NewBond, Molecule)]
halogens = List.zip bond hal
  where hal = List.map makeScaffoldFromSmiles ["F", "Cl", "Br"]
        bond = List.replicate (List.length hal) Single

-- Create a list of all alkyls with up to 'i' carbon atoms.  Point of attachment is
-- atom number 0
alkyls :: Int -> [(NewBond, Molecule)]
alkyls i = let
  carbon = makeScaffoldFromSmiles "C"
  hydrogen = makeScaffoldFromSmiles "H"
  dummy = makeMoleculeFromAtom $ Open Set.empty Set.empty
  mth = Just $ AddMethod
    { firstApply=Nothing
    , lastApply=Nothing
    , selector=openValenceSelector
    , addList=[(Single, carbon)]
    }
  removeDummy = Just $ FilterMethod
    { firstApply=Nothing
    , lastApply=Nothing
    , molFilter=List.map (\m -> removeAtoms m isOpen )
    }

  alks = List.foldr (\enum mols -> enum mols) [dummy] $ List.replicate (i - 1) (>##> mth)
  uniqueAlks = alks >#> addH >#> makeUnique >#> removeDummy
  bond = List.replicate (List.length uniqueAlks) Single
  in List.zip bond uniqueAlks


addMols :: [Molecule] -> Maybe Method
addMols ms = Just $ AddMethod
  { firstApply=Nothing
  , lastApply=makeUnique
  , selector=openValenceSelector
  , addList=List.map (\m -> (Single, m)) ms
  }

addH = Just $ FilterMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , molFilter=List.map (\m -> fillMoleculeValence m )
  }

removeH = Just $ FilterMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , molFilter=List.map (\m -> removeAtoms m isHydrogen )
  }

removeLP = Just $ FilterMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , molFilter=List.map (\m -> removeAtoms m isLonePair )
  }

removeE = Just $ FilterMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , molFilter=List.map (\m -> removeAtoms m isElectron )
  }

stripMol = Just $ FilterMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , molFilter=(\ms -> ms >#> removeH >#> removeLP >#> removeE)
  }

makeUnique = Just $ FilterMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , molFilter=fingerprintFilterBuilder (\m -> writeCanonicalPath m)
  }

