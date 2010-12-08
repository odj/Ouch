{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Enumerate.Formula
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


-- | Defines an enumerable molecular formula data type that can be used to create
-- enumerated sets that conform to particular formula.

module Ouch.Enumerate.Formula (

  )  where

import Ouch.Enumerate.Method
import Ouch.Structure.Molecule
import Ouch.Structure.Atom
import Ouch.Structure.Marker
import Ouch.Input.Smiles
import Data.Map as Map
import Data.List as List
import Data.Set as Set
import Data.Maybe as Maybe
import Control.Applicative hiding ( (<|>) )
import Text.Parsec
import Text.ParserCombinators.Parsec (GenParser)
import Text.Parsec.Language (haskellDef)
import qualified Text.Parsec.Token as P


data Formula = Formula { elements     :: Map Molecule Int
                       , unsaturation :: Int
                       } deriving (Show, Eq, Ord)


pHydrogen = skipMany (noneOf  "H") >>
            option 0 (try $ char 'H' >> (try pInt <|> return (0::Int) ))


pInt = do n <- option (1::Integer) natural
          return $ (fromInteger n::Int)

pAtom = do skipMany (try $ char 'H' >> try pInt)
           atom <- pAtomSymbol
           n    <- pInt
           return (fromSymbol 0 atom, n)

pFormula = do numberH <- lookAhead pHydrogen
              atoms <- manyTill pAtom (notFollowedBy pAtom <|> eof)
              let atomsMap = Map.fromListWith (+) atoms
                  unsat = degreeOfUnsaturation atomsMap numberH
              return $ Formula atomsMap unsat

degreeOfUnsaturation :: (Map Molecule Int) -> Int -> Int
degreeOfUnsaturation m h = (maxH - h) `div`  2
  where adds m = (fromInteger $ (freeValenceAtIndex m 0) - 2)::Int
        maxH = Map.foldrWithKey (\k n acc -> acc + n*(adds k)) 2 m

decompose :: Formula -> [FormulaMonad]
decompose f@Formula {elements=elem, unsaturation=u} = let
  foldFunc k _ b = (FormulaMonad ([k], decrement f k)):b
  in Map.foldWithKey (foldFunc) [] elem

decrement :: Formula -> Molecule -> Formula
decrement f@Formula {elements=elem, unsaturation=u} m = let
  val = Map.lookup m elem
  in case val of
    Nothing -> f
    Just i -> f {elements=dec}
      where dec | i > 1 = Map.insert  m (i-1) elem
                | i == 1 = Map.delete m elem
                | otherwise = Map.delete m elem

isEmpty :: Formula -> Bool
isEmpty f@Formula {elements=elem, unsaturation=u} = 0 == Map.size elem

instance Read Formula where
  readsPrec _ s =
    let pResult = parse pFormula "" s
        formula = case pResult of
          Right result -> result
          Left _       -> Formula Map.empty 0
    in [(formula, "")]

--instance Show Formula where
--  show = undefined



newtype FormulaMonad = FormulaMonad { formulaM :: ([Molecule], Formula) }
                                    deriving (Show)

compose :: FormulaMonad -> ([Molecule] -> Maybe Method) -> [Molecule]
compose fm@FormulaMonad {formulaM = (mols, f)} makeMethod = let
  dec = decompose f
  mapFunc fm'@FormulaMonad {formulaM = (mols', f')}
    = FormulaMonad ((mols >#> (makeMethod mols')), f')
  output | (isEmpty f) = mols
         | otherwise = List.foldr (\a acc -> (compose a makeMethod) ++ acc) [] (List.map mapFunc dec)
  in output

expand :: Formula -> [Molecule]
expand f = let
  fms = decompose f
  mols = List.map (\fm -> compose fm addMols) fms
  in (List.concat mols) >#> makeUnique











