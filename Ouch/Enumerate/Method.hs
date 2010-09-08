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

module Ouch.Enumerate.Method (
    Method(..)
  , (>#>)
  , addMethod
  , classSelector
  , openValenceSelector
  , halogens
  , alkyls
  ) where

import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Ouch.Structure.Marker
import Ouch.Input.Smiles

import Data.Map as Map
import Data.List as List
import Data.Set as Set
import Data.Maybe as Maybe
import Control.Applicative


{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}



data Method   = NoMethod      {firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method}
              | AddMethod     {firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , addList      ::  ([(NewBond, Molecule)])}
              | InsertMethod  {firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , insertRule   ::  (Molecule -> [Int])
                             , insertList   ::  [(Molecule, [Int])]}
              | ReplaceMethod {firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , replaceList  ::  [Molecule]}
              | ReactMethod   {firstApply   ::  Maybe Method
                             , lastApply    ::  Maybe Method
                             , selector     ::  (Molecule -> Atom -> Bool)
                             , reactList    ::  [(Molecule -> Molecule)]}




{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}


(>#>) :: [Molecule] -> (Maybe Method) -> [Molecule]
(>#>) ms mMethod = case mMethod of
  Nothing  -> ms
  Just method -> case method of
    NoMethod      {} -> ms >#> (firstApply method)
    AddMethod     {} -> addMethod ms method
    InsertMethod  {} -> ms >#> (firstApply method)
    ReplaceMethod {} -> replaceMethod ms method
    ReactMethod   {} -> ms >#> (firstApply method)

addMethod :: [Molecule] -> Method -> [Molecule]
addMethod ms method = let
  mols = ms >#> (firstApply method)
  atomList f m = Map.keys $ fst $ Map.partition (f m) (atomMap m)
  atomLists = List.map (atomList $ selector method) mols
  zipped = zip mols atomLists
  makeMol m1 i1 addItem = [connectMoleculesAtIndicesWithBond m1 i1 (snd addItem) 0 (fst addItem)]
  newMols m i = (addList method) >>= (makeMol m i)
  output = concat $ List.map (\(m, l) -> List.map (\i -> newMols m i) l ) zipped
  in concat output


replaceMethod :: [Molecule] -> Method -> [Molecule]
replaceMethod ms method = let
  mols = ms >#> (firstApply method)
  atomList f m = Map.keys $ fst $ Map.partition (f m) (atomMap m)
  atomLists = List.map (atomList $ selector method) mols
  zipped = zip mols atomLists
  output = undefined
  in ms

{------------------------------------------------------------------------------}
{-------------------------------Convenience Functions--------------------------}
{------------------------------------------------------------------------------}


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

-- Convenience to create a new selector from two with AND logic
(>||>) :: (Molecule -> Atom -> Bool) -> (Molecule -> Atom -> Bool) -> (Molecule -> Atom -> Bool)
(>||>) sel1 sel2 = (\m a -> (sel1 m a) || (sel2 m a))

{------------------------------------------------------------------------------}
{-------------------------------List Generators--------------------------------}
{------------------------------------------------------------------------------}

halogens :: [(NewBond, Molecule)]
halogens = zip bond hal
  where hal = List.map makeScaffoldFromSmiles ["F", "Cl", "Br"]
        bond = replicate (length hal) Single

-- Create a list of all alkyls with up to 'i' carbon atoms
alkyls :: Int -> [(NewBond, Molecule)]
alkyls i = let
  carbon = makeScaffoldFromSmiles "C"
  hydrogen = makeScaffoldFromSmiles "H"
  mth = Just $ AddMethod Nothing
                         Nothing
                         openValenceSelector
                         [(Single, carbon), (Single, hydrogen)]
  mths = replicate i (>#> mth)
  alks = List.foldr (\enum mols -> enum mols) [carbon] mths
  alks_noH = List.map (\m -> removeAtoms m isHydrogen ) alks
  bond = replicate (length alks_noH) Single
  in zip bond alks_noH



