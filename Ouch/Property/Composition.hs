{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Composition - a module to manage property data types

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

-------------------------------------------------------------------------------
------------------------------------------------------------------------------}

module Ouch.Property.Composition
  (
   molecularWeight
 , molecularFormula
 , exactMass
 , atomCount
 , heavyAtomCount
 , hBondAcceptorCount
 , hBondDonorCount
 , netCharge
  ) where



import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Molecule
import Ouch.Property.Property

import Data.Maybe
import Data.Set as Set
import Data.List as List
import Data.Map as Map


{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

molecularWeight :: Molecule -> Maybe Property
molecularWeight m = undefined


molecularFormula :: Molecule -> Maybe Property
molecularFormula m = undefined


exactMass :: Molecule -> Maybe Property
exactMass m = undefined


atomCount :: Molecule -> Maybe Property
atomCount m = undefined


heavyAtomCount :: Molecule -> MaybeProperty
heavyAtomCount m = undefined


hBondAcceptorCount :: Molecule -> Maybe Property
hBondAcceptorCount m = undefined


hBondDonorCount :: Molecule -> Maybe Property
hBondDonorCount m = undefined


netCharge :: Molecule -> Maybe Property
netCharge m = undefined


--molecularWeight
{------------------------------------------------------------------------------}
molecularWeight :: Molecule -> Either String Double
molecularWeight m = if (moleculeHasError m) then (Left "") else case m of
        Small {atomMap=atoms} -> Right $ mw atoms
        Markush  {}   -> Left "No MW for Markush"
        Polymer  {}   -> Left "No MW for Polymer"
        Biologic {}   -> Left "No MW for Biologic"
        where mw a = foldl (+) 0.0 $ List.map atomMW $ Map.fold (\a -> (++) [a]) [] a




--molecularFormula
{------------------------------------------------------------------------------}
molecularFormula :: Molecule -> Either String String
molecularFormula m = if (moleculeHasError m) then (Left "") else case m of
    Small {atomMap=at}    -> Right molFm
        where startMap = Map.empty
              endMap = List.foldr (updateMap) startMap $ List.map snd $ Map.toList at
              -- Use foldr to accumulate and count atoms
              updateMap a m | Map.notMember (atomicSymbolForAtom a) m = Map.insert (atomicSymbolForAtom a)  1 m
                            | otherwise                               = Map.adjust (+ 1) (atomicSymbolForAtom a) m
              -- Convert the map to a list of just the elements present, and in IUPAC order
              finalList = catMaybes $ List.map (\e -> lookupPair e endMap) molecularFormulaElements
              --  Build the final output string from the map
              molFm = List.foldr (\(e,n) -> if n>1 then ((e ++  (show n))++) else (e ++ ))  "" finalList
              -- simple little utility function which, strangely, is not already defined in Data.Map
              lookupPair k m = case v of
                  Just val -> Just (k, val)
                  Nothing -> Nothing
                  where v = Map.lookup k m
    Markush  {}   -> Left "No molecular formula defined for Markush"
    Polymer  {}   -> Left "No molecular formula defined for Polymer"
    Biologic {}   -> Left "No molecular formula defined for Biologic"
















