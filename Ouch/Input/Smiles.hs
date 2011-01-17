{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Input.Smiles
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
{-# LANGUAGE NoMonomorphismRestriction #-} --Need this for ParserCombinators to work

-- | Handles parsing of a SMILES string into a 'Molecule' data type.
module Ouch.Input.Smiles (
     readSmi
   , makeScaffoldFromSmiles
   , natural
   , pAtomSymbol
   , fromSymbol
   ) where


import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Data.Atom
import Data.Maybe
import Data.Char
import Data.Set as Set
import Data.List as List
import Data.Map as Map
import Control.Applicative hiding ((<|>), optional, many)
import Text.Parsec
import Text.ParserCombinators.Parsec (GenParser)
import Text.Parsec.Language (haskellDef)
import qualified Text.Parsec.Token as P

{------------------------------------------------------------------------------}
{------------------------------------------------------------------------------}


-- | Return the scaffold 'Molecule' only with no implicit hydrogens added
makeScaffoldFromSmiles :: String -> Molecule
makeScaffoldFromSmiles s = case (parse pSmiles "" s) of
              Right (m, b) -> cyclizeMolecule m
              Left      er -> giveMoleculeError emptyMolecule (show er)

-- | Returns the complete 'Molecule' for the 'String' with implicit hydrogens
-- added as needed.  If parsing fails, the 'Molecule' is returned with
-- an internal error that gives the full text of the 'Parsec' error string
readSmi :: String -> Molecule
readSmi s = fillMoleculeValence $ makeScaffoldFromSmiles s

-- Define what we can from an established lexer (any one would do)
lexer = P.makeTokenParser haskellDef
natural = P.natural lexer
bracket = between (char '[') (char ']')

pSmiles = (emptyMolecule, Single) <$ char ')' <|>
          (emptyMolecule, Single) <$ eof <|>
            do bond <- pBond
               geometry <- optionMaybe pGeometry
               atom  <- pAtom
               subsmiles <- many (char '(' *> pSmiles)
               atoms <- pSmiles
               let branched = List.foldr addSub atom subsmiles
               return (addSub atoms branched, bond)

addSub (smi, bnd) mol = connectMoleculesAtIndicesWithBond mol 0 smi 0 bnd

pAtom = do atom    <- (try pSubsetAtom <|> try pBracket)
           closure <- many $ try pClosure
           return $ List.foldr
                    (\mk a -> addMarkerToAtomAtIndex a 0 mk) atom closure

pBracket = bracket $ do isotope    <- pIsotope
                        atomSymbol <- pAtomSymbol
                        chiral     <- optionMaybe pStereo
                        explicitH  <- optionMaybe pHydrogen
                        atomCharge <- optionMaybe pCharge
                        atomClass  <- optionMaybe pClass
                        return $ (fromSymbol isotope atomSymbol)
                                 >@> chiral >@> explicitH >@> atomCharge
                                 >@> atomClass >@> aromatic atomSymbol

pSubsetAtom = do atomSymbol <- (try pOrganicSubsetSymbols)
                 return $ fromSymbol 0 atomSymbol >@> (aromatic atomSymbol)

lowerFirst s = (toLower $ head s):(tail s)
upperFirst s = (toUpper $ head s):(tail s)
aromatic s | (isLower $ head s) = Just AromaticAtom | otherwise = Nothing

fromSymbol :: Integer -> String -> Molecule
fromSymbol i s = makeMoleculeFromAtom $ Element n (i - n) Set.empty Set.empty
                 where n = case Map.lookup (upperFirst s) atomicNumberFromSymbol of
                              Just atomicN -> atomicN
                              Nothing -> 0

organicSubsetList = ["Br", "Cl", "Si", "C", "N", "O", "H", "P", "S", "F", "B", "I"]
aromaticSubset = ["b", "c", "n", "o", "s", "p"]

pOrganicSubsetSymbols = choice $ (List.map (try. string) organicSubsetList)
                              ++ (List.map (try. string ) aromaticSubset)

pAtomSymbol = choice $ (List.map (try . string) smilesSymbols)
                    ++ (List.map (try . string . lowerFirst) smilesSymbols)

pClass = char ':' >> (Class <$> natural)

pClosure = do bond <- pBond
              geometry <- optionMaybe pGeometry  -- Not yet implemented semantically
              closure <- ((\s -> (read [s])::Integer) <$> digit) <|>
                         (char '%' >> natural)
              return $ Closure (fromIntegral closure) bond

pCharge = pPlus <|> pMinus

pPlus  = (try $ char '+' >> Charge <$> natural) <|>
         (try $ Charge 2 <$ string "++")        <|>
         (Charge 1 <$ char '+')

pMinus = (try $ char '-' >> Charge <$> natural) <|>
         (try $ Charge (-2) <$ string "--")     <|>
         (Charge (-1) <$ char '-')

pBond = option Single $
        NoBond <$ char '.' <|>
        Single <$ char '-' <|>
        Double <$ char '=' <|>
        Triple <$ char '#'

pStereo = (try $ (Chiral Smiles2) <$ string "@@") <|>
                 (Chiral Smiles1) <$ string "@"

pGeometry = char '\\' <|> char '/' -- Not yet implemented semantically

pHydrogen = option (ExplicitHydrogen 0) $
            (try $ char 'H' >> ExplicitHydrogen <$> natural) <|>
            (try $ ExplicitHydrogen 1 <$ char 'H')

pIsotope = option 0 natural




