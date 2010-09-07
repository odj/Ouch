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
) where

import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Ouch.Structure.Marker

import Data.Map as Map
import Data.List as List
import Data.Set as Set
import Data.Maybe as Maybe
import Control.Applicative


{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}



data Method   = NoMethod      {firstApply::Maybe Method}
              | AddMethod     {firstApply::Maybe Method
                             , selector::(Molecule -> Atom -> Bool)
                             , addList::([(NewBond, Molecule)])}
              | InsertMethod  {firstApply::Maybe Method}
              | ReplaceMethod {firstApply::Maybe Method}
              | ReactMethod   {firstApply::Maybe Method}




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
    ReplaceMethod {} -> ms >#> (firstApply method)
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



