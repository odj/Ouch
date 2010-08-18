{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Smiles - a module to parse smiles strings

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

module Ouch.Input.Smiles (
     nextChoppedSmile
   , makeMoleculeFromSmiles
   , makeScaffoldFromSmiles
   , makeAtomMoleculeFromChop
   , findNextSubSmile
   , parseSmiles
   , parseClosureAtomMarkers
   , makeAtomMoleculeFromBracketChop
   , chop
   ) where


import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Data.Atom
import Text.Regex.TDFA ((=~))
import Data.Maybe
import Data.Char
import Data.Set as Set
import Data.List as List
import Data.Map as Map
import Control.Applicative



{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}



data ChoppedSmile = Smile        {smile::String
                                , smiles::String
                                , newBond::NewBond
                                , mark::(Set AtomMarker)}
                  | SubSmile    {smile::String
                                , smiles::String
                                , newBond::NewBond
                                , mark::(Set AtomMarker)}
                                --Substructures don't technically
                                --need a marker set
                  | SmilesError {smile::String
                                , smiles::String
                                , newBond::NewBond
                                , mark::(Set AtomMarker)}
                  deriving (Show, Eq)



{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

-- makeMoleculeFromSmiles
-- Add all the "cleanup" functions here.
{------------------------------------------------------------------------------}
makeMoleculeFromSmiles::String -> Molecule
makeMoleculeFromSmiles smi = mol' {molMarkerSet=newMolInfo}
    where mol'  =   fillMoleculeValence
                  $ cyclizeMolecule
                  $ makeScaffoldFromSmiles smi
          info  = Set.singleton (Info $ "Produced from smile string: " ++ smi)
          closureWarning | hasHangingClosure mol'
                         = Set.singleton
                         $ Warning "Molecule has unmatched bond closures."
                         | otherwise
                         = Set.empty
          newMolInfo = foldr Set.union Set.empty
                       [(molMarkerSet mol'), info, closureWarning]




-- makeScaffoldFromSmiles
-- Right if successful
-- Left if error somewhere, with brief error description
-- This functions kicks off the recursive interpretations
{------------------------------------------------------------------------------}
makeScaffoldFromSmiles::String -> Molecule
makeScaffoldFromSmiles smi = case chop of
    Smile {}        -> growMoleculeAtIndexWithString newAtom 0 nextSmile
--    A well-formed smile should never start with a substructure, but if it does,
--    we connect the next atom or substructure to it's beginning (which may or
--    may not be permitted chemically)
    SubSmile {}     -> growMoleculeAtIndexWithString newSubstructure 0 nextSmile
    SmilesError {}  -> giveMoleculeError emptyMol
                       ("Error trying to parse the Smiles string: " ++ smi)
    where chop = nextChoppedSmile smi
          thisSmile = smile chop
          nextSmile = smiles chop
          newAtom = makeAtomMoleculeFromChop chop
          newSubstructure = makeScaffoldFromSmiles thisSmile
          emptyMol = Small Map.empty Set.empty




-- growPerhapsMoleculeWithString
-- Adds continuation of smiles string at the end of the molecule being built
{------------------------------------------------------------------------------}
growMoleculeAtIndexWithString :: Molecule -> Int -> String -> Molecule
growMoleculeAtIndexWithString m i smi
    | smi == ""  = m
    | otherwise  = if (moleculeHasError m) then m else case chop of
            Smile {}        -> growMoleculeAtIndexWithString newMolecule1
                               (newIndex) nextSmile
            SubSmile {}     -> growMoleculeAtIndexWithString newMolecule2
                               i nextSmile
            SmilesError {}  -> giveMoleculeError m
                               ("Error trying to grow from the Smiles string: " ++ smi)
            where chop = nextChoppedSmile smi
                  nextSmile = smiles chop
                  newAtom = makeAtomMoleculeFromChop chop
                  newSubStructure = makeScaffoldFromSmiles (smile chop)
                  newMolecule1 = connectMoleculesAtIndicesWithBond
                                 m i newAtom 0 (newBond chop)
                  newMolecule2 = connectMoleculesAtIndicesWithBond
                                 m i newSubStructure 0 (newBond chop)
                  newIndex = Map.size (atomMap newMolecule1) - 1



-- makeScaffoldFromSmiles
-- Right if successful
-- Left if error somewhere, with brief error description
{------------------------------------------------------------------------------}
makeAtomMoleculeFromChop :: ChoppedSmile -> Molecule
makeAtomMoleculeFromChop nb = case nb of
    SmilesError {} -> giveMoleculeError emptyMol
                      ("Could not make molecule from smile with string: "
                      ++ smile nb)
    Smile {smile=s} -> if ((head s) == '[')
                       then makeAtomMoleculeFromBracketChop nb
                       else makeAtomMolecule nb
    where markSetType = Set.singleton (if isLower $ head (smile nb)
                                       then AromaticAtom
                                       else Null)
          markSetClass = Set.singleton $ Class 0
          markSetAll = markSetType `Set.union`
                       (mark nb) `Set.union`
                       markSetClass
          a = [toUpper c | c <- smile nb]
          emptyMol = Small Map.empty Set.empty
          makeAtomMolecule nb
           | a == ""    =  giveMoleculeError emptyMol
                           "ERROR: Tried to make atom from empty string."
           | a == "C"   =  makeMoleculeFromAtom $ Element 6 0  Set.empty markSetAll
           | a == "N"   =  makeMoleculeFromAtom $ Element 7 0  Set.empty markSetAll
           | a == "O"   =  makeMoleculeFromAtom $ Element 8 0  Set.empty markSetAll
           | a == "H"   =  makeMoleculeFromAtom $ Element 1 0  Set.empty markSetAll
           | a == "P"   =  makeMoleculeFromAtom $ Element 15 0 Set.empty markSetAll
           | a == "S"   =  makeMoleculeFromAtom $ Element 16 0 Set.empty markSetAll
           | a == "F"   =  makeMoleculeFromAtom $ Element 9  0 Set.empty markSetAll
           | a == "B"   =  makeMoleculeFromAtom $ Element 5  0 Set.empty markSetAll
           | a == "BR"  =  makeMoleculeFromAtom $ Element 35 0 Set.empty markSetAll
           | a == "CL"  =  makeMoleculeFromAtom $ Element 17 0 Set.empty markSetAll
           | a == "I"   =  makeMoleculeFromAtom $ Element 53 0 Set.empty markSetAll
           | a == "*"   =  makeMoleculeFromAtom $ Unspecified  Set.empty markSetAll
           | otherwise  =  giveMoleculeError emptyMol
                           ("ERROR: Atom not recognized for symbol: " ++ a)



makeAtomMoleculeFromBracketChop::ChoppedSmile -> Molecule
makeAtomMoleculeFromBracketChop sb = mol
    where (s1, s2, s3)    = (smile sb)
                          =~ "(^(\\[)([+-@:a-zA-Z0-9])*(\\])([-=#\\.]{0,1}[/\\]{0,1}[0-9])*)"
                          ::(String, String, String)
          s =  init $ tail s2   -- Drop the brackets

          -- Get charge information
          (ch1, ch2, ch3)   = s
                            =~ "([+-][0-9]*)"
                            ::(String, String, String)
          (sign, int, _)    = ch2
                            =~ "([^+-][0-9]*)"
                            ::(String, String, String)
          charge | int == "" = 1 | otherwise = read int::Integer
          markCharge | ch2 == "" = Null
                     | sign == "-" = Charge (0 - charge)
                     | otherwise = Charge charge

          -- Get isotopic number
          (n1, n2, n3)   = s =~ "(^[0-9]*)"
                             ::(String, String, String)
          isotope | n2 == "" = 0 | otherwise = read n2::Integer

          -- Get element symbol
          (a1, a2, a3)   = s =~ "([A-Za-z]{1}[a-z]{0,1})"
                             ::(String, String, String)
          lookupString | isLower (head a2) && (length a2) == 1
                       = [toUpper $ a2!!0] | otherwise = a2
          markAromatic | isLower (head a2) && (length a2) == 1
                       = AromaticAtom | otherwise = Null
         -- Need to error check here!
          atomicNumber = case Map.lookup lookupString atomicNumberFromSymbol of
              Just n -> n
              Nothing -> 0  -- This needs to generate an error

          -- Get number of extra internal hydrogens (I think this notation is stupid)
          (h1, h2, h3)    = a3 =~ "([H]{1}[0-9]*)"
                               ::(String, String, String)
          (nh1, nh2, nh3) = h2 =~ "([0-9]+)"
                               ::(String, String, String)
          numberH | h2 == "" = 0 | nh2 == "" = 1 | otherwise = read nh2::Integer
          markH = ExplicitHydrogen numberH

          -- Get class information
          (c1, c2, c3)       = s  =~ "([:][0-9]+)"
                                  ::(String, String, String)
          (cn1, cn2, cn3)    = c2 =~ "([0-9]+)"
                                  ::(String, String, String)
          classNumber | c2 == "" = 0
                      | cn2 == "" =0
                      | otherwise = read cn2::Integer
          markClass = Class classNumber

          -- Get Stereochemical Information
          scCount = foldr (\c n -> if (c == '@') then n+1 else n) 0 s
          markStereo | scCount == 1 = Null -- (Chiral Levo)
                     | scCount == 2 = Null -- (Chiral Dextro)
                     | otherwise = Null

          -- Now make the molecule
          mol = makeMoleculeFromAtom $ Element atomicNumber (isotope - atomicNumber) Set.empty markSetAll
          markSetAll = Set.union (mark sb)
                     $ Set.fromList ([markAromatic, markClass, markH, markStereo, markCharge]
                     ++ (parseClosureAtomMarkers s3 []))

-- nextSmilesSubstring
-- Lots, lots, lots!!!!! more to fill in here
{------------------------------------------------------------------------------}
parseSmiles :: String -> (String, String, String)
parseSmiles s = s =~ pat::(String, String, String)
    -- Need a more general format to recognize two-letter elements (or maybe just enumerate them?)
    where pat = foldr (++) "" patList
          patList = [ "("
                    -- Search for first atom + bond/marker
                    , "^[-=#\\./]{0,1}"
                   ++ "[\\]{0,1}([A-Za-z]|Br|Cl|Si|Sn|Li|Na|Cs){1}"
                   ++ "([-=#\\.]{0,1}[/\\]{0,1}[0-9])*[@]*"
                    -- return the whole string for anything
                    -- that STARTS with open parens
                    , "|^[(].*|(^[/][(].*)|(^[\\][(].*)"
                    -- return next atom segment within square
                    -- brackets plus closure ID's afterwards
                    , "|(^([-=#\\.][\\[]|[/][\\[]|[\\\\[])"
                   ++ "([+-@:a-zA-Z0-9])*(\\])([-=#\\.]{0,1}"
                   ++ "[/\\]{0,1}[0-9])*)"
                    , ")"
                    ]


-- nextSmilesSubstring
-- Function strips atom/bond info into a new type.
{------------------------------------------------------------------------------}
nextChoppedSmile :: String -> ChoppedSmile
nextChoppedSmile s
   | s2 == "" || s1 /= "" = SmilesError {smile=s1, smiles=s, newBond=NoBond, mark=(Set.empty)}
   | isSubSmile           = SubSmile    {smile=b', smiles=ss3, newBond=nb', mark=markerSet}
   | otherwise            = Smile       {smile=smi, smiles=s3, newBond=nb, mark=markerSet}

   where (s1, s2, s3)       = parseSmiles s   -- Get initial parse

         -- Strip off double-bond geometry info first!
         (g1, g2, g3)       = s2 =~ "(^[/\\])"::(String, String, String)
         markGeo            | g2 == "/"   = Set.singleton (GeoIsomer ProCis)
                            | g2 == "\\"  = Set.singleton (GeoIsomer ProTrans)
                            | g2 == ""    = Set.singleton Null
                            | otherwise   = Set.singleton Null
        -- Use the remainder that contains information for all operations below
         g                   | length g1 /= 0 = g1 | length g3 /= 0 = g3

        -- What kind of initial parse we get determines how to handle it.
         isSubSmile | length g /= 0 = head g == '(' | otherwise = False
         (ss1, ss2, ss3)  | length g /= 0 = findNextSubSmile g 1
                          | otherwise = ("", "", "")

         -- Get bond info for single atoms (or bracketed atoms)
         (b1, b2, b3)       = g =~ "(^[-=#\\.])"::(String, String, String)
         b                   | length b1 /= 0 = b1
                             | length b3 /= 0 = b3
         nb   | b2 == "."              = NoBond
              | b2 `elem` ["", "-"]    = Single
              | b2 == "="              = Double
              | b2 == "#"              = Triple
              | otherwise              = Single
        -- Get bond info for subsmiles
         (b1', b2', b3')    = ss2 =~ "(^[-=#\\.])"::(String, String, String)
         b'                   | length b1' /= 0 = b1'
                              | length b3' /= 0 = b3'
                              -- | otherwise = ""
         nb'   | b2' == "."            = NoBond
               | b2' `elem` ["", "-"]  = Single
               | b2' == "="            = Double
               | b2' == "#"            = Triple
               | otherwise             = Single

         -- Get atom closure info and parse it,
         -- then merge with above geoAtomMarker
         markerSet       = Set.union markGeo (Set.fromList
                         $ parseClosureAtomMarkers lb2 [])
         (tr1, tr2, tr3) = b   =~ "(.*\\])|()"
                               ::(String, String, String)
         (lb1, lb2, lb3) = tr3 =~ "(([-=#\\.]{0,1}[%]{0,1}[/\\]{0,1}[0-9])+)"
                               ::(String, String, String)
         smi | tr2 == ""  = lb1
             | otherwise  = tr2

--parseClosureAtomMarkers
parseClosureAtomMarkers :: String -> [AtomMarker] -> [AtomMarker]
parseClosureAtomMarkers [] ml = ml
parseClosureAtomMarkers s ml = parseClosureAtomMarkers s3
                               (Closure {labelNumber=closureNumber, bondType=nb}:ml)
    where (_, s2, s3) = s  =~ "([-=#\\.]{0,1}[/\\]{0,1}[%]{0,1}[0-9]{0,1}){0,1}"
                           ::(String, String, String)
          (_, n, _)   = s  =~ "([0-9])"
                           ::(String, String, String)
          (_, nbs, _) = s2 =~ "^[-=#\\.]"
                           ::(String, String, String)
          nb | nbs == "-" = Single
             | nbs == "=" = Double
             | nbs == "#" = Triple
             | otherwise = Single
          closureNumber = read n::Int


-- findNextSubSmile
{------------------------------------------------------------------------------}
findNextSubSmile::String -> Int -> (String, String, String)
findNextSubSmile s i | length s >= i && (open == closed)  = ("", s2', s3)
                     | length s > i                       = findNextSubSmile s (i + 1)
                     | length s <= i                      = ("", "", s)
                     | otherwise                          = ("", "", s)
                     where s2 = take i s
                           s3 = drop i s
                           s2' | i > 2 = tail $ init s2 | otherwise = ""
                           open = sum [1 | c <- s2, c=='(']
                           closed = sum [1 | c <- s2, c==')']



--  debugging function
{------------------------------------------------------------------------------}
chop :: String -> String
chop [] = ""
chop a = show nb ++ "\n" ++ chop (smiles nb)
   where nb = nextChoppedSmile a

