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
   , parseClosureMarkers
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
import qualified Data.List as List
import Data.Map as Map
import Control.Applicative 



{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}



data ChoppedSmile = Smile {smile::String, smiles::String, newBond::NewBond, mark::(Set Marker)}
                  | SubSmile {smile::String, smiles::String, newBond::NewBond, mark::(Set Marker)}  --Substructures don't need a marker set
                  | SmilesError {smile::String, smiles::String, newBond::NewBond, mark::(Set Marker)}  --
                  deriving (Show, Eq)



{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

-- makeMoleculeFromSmiles
-- Add all the "cleanup" functions here.
{------------------------------------------------------------------------------}
makeMoleculeFromSmiles::String -> PerhapsMolecule
makeMoleculeFromSmiles smi = fillMoleculeValence $ makeScaffoldFromSmiles smi
    
    

-- makeScaffoldFromSmiles
-- Right if successful
-- Left if error somewhere, with brief error description
-- This functions kicks off the recursive interpretations
{------------------------------------------------------------------------------}
makeScaffoldFromSmiles::String -> PerhapsMolecule
makeScaffoldFromSmiles smi = case chop of 
    Smile {}        -> growPerhapsMoleculeAtIndexWithString newAtom 0 nextSmile 
    -- A well-formed smile should never start with a substructure, but if it does, we connect the next 
    -- atom or substructure to it's beginning (which may or may not be permitted chemically)
    SubSmile {}     -> growPerhapsMoleculeAtIndexWithString newSubstructure 0 nextSmile  
    SmilesError {}  -> (Left $ "Error trying to parse the Smiles string: " ++ smi)
    where chop = nextChoppedSmile smi
          thisSmile = smile chop
          nextSmile = smiles chop
          newAtom = makeAtomMoleculeFromChop chop
          newSubstructure = makeScaffoldFromSmiles thisSmile
          

                            

-- growPerhapsMoleculeWithString
-- Adds continuation of smiles string at the end of the molecule being built
{------------------------------------------------------------------------------}
growPerhapsMoleculeAtIndexWithString :: PerhapsMolecule -> Int -> String -> PerhapsMolecule
growPerhapsMoleculeAtIndexWithString pm i smi 
    | smi == ""  = pm
    | otherwise  = case pm of
    Right {} -> case chop of
            Smile {}        -> growPerhapsMoleculeAtIndexWithString newMolecule1 (newIndex) nextSmile 
            SubSmile {}     -> growPerhapsMoleculeAtIndexWithString newMolecule2 i nextSmile
            SmilesError {}  -> (Left $ "Error trying to grow from the Smiles string: " ++ smi)
            where chop = nextChoppedSmile smi
                  nextSmile = smiles chop
                  newAtom = makeAtomMoleculeFromChop chop
                  newSubStructure = makeScaffoldFromSmiles (smile chop)
                  newMolecule1 = connectPerhapsMoleculesAtIndicesWithBond pm i newAtom 0 (newBond chop)
                  newMolecule2 = connectPerhapsMoleculesAtIndicesWithBond pm i newSubStructure 0 (newBond chop)
                  -- We just made this thing, so there shouldn't be any errors, right?
                  newIndex = case newMolecule1 of Right m1 -> Map.size (atomMap m1) - 1          
    Left {} -> pm


-- makeScaffoldFromSmiles
-- Right if successful
-- Left if error somewhere, with brief error description
{------------------------------------------------------------------------------}
makeAtomMoleculeFromChop::ChoppedSmile -> PerhapsMolecule
makeAtomMoleculeFromChop nb = case nb of 
    SmilesError {} -> Left $ "Could not make molecule from smile with string: " ++ smile nb
    Smile {smile=s} -> if ((head s) == '[') then makeAtomMoleculeFromBracketChop nb else makeAtomMolecule nb
        where makeAtomMolecule nb  | a == ""        =  Left "ERROR: Tried to make atom from empty string."
                                   | a == "C"       =  Right $ Small $ Map.singleton 0 $ Element 6 0 [] markSetAll
                                   | a == "N"       =  Right $ Small $ Map.singleton 0 $ Element 7 0 [] markSetAll
                                   | a == "O"       =  Right $ Small $ Map.singleton 0 $ Element 8 0 [] markSetAll
                                   | a == "H"       =  Right $ Small $ Map.singleton 0 $ Element 1 0 [] markSetAll

                                   | a == "P"       =  Right $ Small $ Map.singleton 0 $ Element 15 0 [] markSetAll
                                   | a == "S"       =  Right $ Small $ Map.singleton 0 $ Element 16 0 [] markSetAll
                                   | a == "F"       =  Right $ Small $ Map.singleton 0 $ Element 9 0 [] markSetAll
                                   | a == "B"       =  Right $ Small $ Map.singleton 0 $ Element 5 0 [] markSetAll
                                   | a == "BR"      =  Right $ Small $ Map.singleton 0 $ Element 35 0 [] markSetAll
                                   | a == "CL"      =  Right $ Small $ Map.singleton 0 $ Element 17 0 [] markSetAll
                                   | a == "I"       =  Right $ Small $ Map.singleton 0 $ Element 53 0 [] markSetAll
                                   | a == "*"       =  Right $ Small $ Map.singleton 0 $ Unspecified [] markSetAll -- Wildcard Atom
                                   | otherwise      =  Left  $ "ERROR: Atom not recognized for symbol: " ++ a
                                   where markSetType = Set.singleton (if isLower $ head (smile nb) then AromaticAtom else Null)
                                         markSetClass = Set.empty 
                                         markSetAll = markSetType `Set.union` mark nb
                                         a = [toUpper c | c <- smile nb] 


makeAtomMoleculeFromBracketChop::ChoppedSmile -> PerhapsMolecule
makeAtomMoleculeFromBracketChop sb = mol
    where (s1, s2, s3)    = (smile sb) =~ "(^[\\[].*(\\]))"::(String, String, String)  -- Contents of the brackets
          s =  init $ tail s2   -- Drop the brackets
          
          -- Get charge information
          (ch1, ch2, ch3)   = s     =~ "([+-][0-9]*)"::(String, String, String)
          (sign, int, _)    = ch2   =~ "([^+-][0-9]*)"::(String, String, String)
          charge | int == "" = 1 | otherwise = read int::Integer
          markCharge | ch2 == "" = Null | sign == "-" = Charge (0 - charge) | otherwise = Charge charge 
          
          -- Get isotopic number
          (n1, n2, n3)   = s     =~ "(^[0-9]*)"::(String, String, String)
          isotope | n2 == "" = 0 | otherwise = read n2::Integer
          
          -- Get element symbol
          (a1, a2, a3)   = s     =~ "([A-Za-z]{1}[a-z]{0,1})"::(String, String, String)
          lookupString | isLower (head a2) && (length a2) == 1 = [toUpper $ a2!!0] | otherwise = a2
          markAromatic | isLower (head a2) && (length a2) == 1 = AromaticAtom | otherwise = Null
          atomicNumber = case Map.lookup lookupString atomicNumberFromSymbol of  -- Need to error check here!
              Just n -> n
              Nothing -> 0  -- This needs to generate an error
          
          -- Get number of extra internal hydrogens (I think this notation is stupid)
          (h1, h2, h3)    = a3     =~ "([H]{1}[0-9]*)"::(String, String, String)
          (nh1, nh2, nh3) = h2     =~ "([0-9]+)"::(String, String, String)
          numberH | h2 == "" = 0 | nh2 == "" = 1 | otherwise = read nh2::Integer
          markH = ExplicitHydrogen numberH
          
          -- Get class information
          (c1, c2, c3)       = s     =~ "([:][0-9]+)"::(String, String, String)
          (cn1, cn2, cn3)    = c2    =~ "([0-9]+)"::(String, String, String)
          classNumber | c2 == "" = 0 | cn2 == "" =0 | otherwise = read cn2::Integer
          markClass = Class classNumber
          
          -- Get Stereochemical Information
          scCount = List.foldr (\c n -> if (c == '@') then n+1 else n) 0 s
          markStereo | scCount == 1 = Null -- (Chiral Levo)
                     | scCount == 2 = Null -- (Chiral Dextro)
                     | otherwise = Null

          -- Now make the molecule
          mol = Right $ Small $ Map.singleton 0 $ Element atomicNumber (isotope-atomicNumber) [] markSetAll
         -- hydrogen = Right $ Small $ Map.singleton 0 $ Element 1 0 [] Set.empty
          -- hydrogens = take (fromIntegral numberH) $ repeat hydrogen
          -- Need to fill the rest with radicals to keep the bracket designation explicit.
          -- molH = List.foldr (\a mol -> connectPerhapsMoleculesAtIndicesWithBond mol 0 a 0 Single) mol hydrogens
          
         
          markSetAll = Set.union (mark sb) $ Set.fromList ([markAromatic, markClass, markH, markStereo, markCharge] ++ (parseClosureMarkers s3 [])) 

-- nextSmilesSubstring
-- Lots, lots, lots!!!!! more to fill in here
{------------------------------------------------------------------------------}
parseSmiles :: String -> (String, String, String)
parseSmiles s = s =~ pat::(String, String, String)
    -- Need a more general format to recognize two-letter elements (or maybe just enumerate them?) 
    where pat = List.foldr (++) "" patList 
          patList = [ "("
                    , "^[-=#\\./]{0,1}[\\]{0,1}([A-Za-z]|Br|Cl|Si|Sn|Li|Na|Cs){1}([-=#\\.]{0,1}[/\\]{0,1}[0-9])*[@]*"  -- Search for first atom + bond/marker
                    , "|^[(].*|(^[/][(].*)|(^[\\][(].*)" -- return the whole string for anything that STARTS with open parens
                    , "|(^[\\[]([+-@:a-zA-Z]|[0-9])*(\\])([-=#\\.]{0,1}[/\\]{0,1}[0-9])*)"  -- return next atom segment within square brackets plus closure ID's afterwards
                    , "|(^[-=#\\.]{0,1}[\\[]([+-@:a-zA-Z]|[0-9])*(\\])([-=#\\.]{0,1}[/\\]{0,1}[0-9])*)"
                    , "|(^[/][\\[]([+-@:a-zA-Z]|[0-9]|])*(\\])([-=#\\.]{0,1}[/\\]{0,1}[0-9])*)"
                    , "|(^[\\][\\[]([+-@:a-zA-Z]|[0-9]|])*(\\])([-=#\\.]{0,1}[/\\]{0,1}[0-9])*)"
                    , ")"
                    ]

                          
-- nextSmilesSubstring
-- Function strips atom/bond info into a new type.  
{------------------------------------------------------------------------------}
nextChoppedSmile :: String -> ChoppedSmile 
nextChoppedSmile s
   | s2 == "" || s1 /= ""   = SmilesError {smile=s1, smiles=s, newBond=NoBond, mark=Set.singleton $ Comment s2}
   | (take 2 s2) == "/["    = Smile {smile=(drop 1 s2), smiles=s3, newBond=nb, mark=(Set.fromList [Comment s2, GeoIsomer ProCis])}
   | (take 2 s2) == "\\["   = Smile {smile=(drop 1 s2), smiles=s3, newBond=nb, mark=(Set.fromList [Comment s2, GeoIsomer ProTrans])}
   | (take 2 s2) == "/("    = SubSmile {smile=bb2', smiles=ss3', newBond=nb2', mark=(Set.fromList [Comment s2, GeoIsomer ProCis])}
   | (take 2 s2) == "\\("   = SubSmile {smile=bb2', smiles=ss3', newBond=nb2', mark=(Set.fromList [Comment s2, GeoIsomer ProTrans])}
   | head s2 == '['         = Smile {smile=s2, smiles=s3, newBond=nb, mark=(Set.singleton $ Comment s2)}
   | head s2 == '('         = SubSmile {smile=bb2, smiles=ss3, newBond=nb2, mark=Set.singleton $ Comment ss2}
   | head s2 == '/'         = Smile {smile=a2, smiles=s3, newBond=nb, mark=(Set.fromList [Comment s2, GeoIsomer ProCis])}
   | head s2 == '\\'        = Smile {smile=a2, smiles=s3, newBond=nb, mark=(Set.fromList [Comment s2, GeoIsomer ProTrans])}
   | head (take 2 s2) == '['= Smile {smile=(drop 1 s2), smiles=s3, newBond=nb, mark=(Set.singleton $ Comment s2)}
   | otherwise              = Smile {smile=a2, smiles=s3, newBond=nb, mark=markerSet `Set.union` (Set.singleton $ Comment s2)}
   
   where (s1, s2, s3)       = parseSmiles s                                 -- Get initial parse
         (b1, b2, b3)       = s2 =~ "(^[-=#\\.])"::(String, String, String) -- Get bond info for single atoms
         (lb1, lb2, lb3)    = s2 =~ "([-=#\\.]{0,1}[%]{0,1}[/\\]{0,1}[0-9])+"::(String, String, String) -- Get atom closure bond substring
         (a1, a2, a3)       = s2 =~ "([A-Za-z]+)"::(String, String, String) -- Atom only, remove bond comments, etc

         (ss1, ss2, ss3)    = findNextSubSmile s 1
         (ss1', ss2', ss3') = findNextSubSmile (drop 1 s) 1    

         (bb1, bb2, bb3)    = ss2 =~ "([\\[A-Za-z]+.*)"::(String, String, String) -- Get bond info for subsmiles
         (bb1', bb2', bb3') = ss2' =~ "([\\[A-Za-z]+.*)"::(String, String, String) -- Get bond info for subsmiles with geometry marker     
         
         nb  | b2 == "."              = NoBond
             | b2 `elem` ["", "-"]    = Single
             | b2 == "="              = Double
             | b2 == "#"              = Triple
             | otherwise              = Single
             
         -- need to process a different string to determine bond type if parse gives a subsmile
         nb2 | bb1 == "."             = NoBond
             | bb1 == "" || bb2 == "-"= Single
             | bb1 == "="             = Double
             | bb1 == "#"             = Triple
             | otherwise              = Single
             
         -- need to process a different string to determine bond type if parse gives a subsmile
         nb2'| bb1' == "."             = NoBond
             | bb1' == "" || bb2' == "-"= Single
             | bb1' == "="             = Double
             | bb1' == "#"             = Triple
             | otherwise              = Single
         -- This is a very simple parse of ring closure markers.  
         -- Does not accomodate "%" notation (yet)
         markerSet = Set.fromList $ parseClosureMarkers lb2 []





parseClosureMarkers :: String -> [Marker] -> [Marker]
parseClosureMarkers [] ml = ml
parseClosureMarkers s ml = parseClosureMarkers s3 (Closure {labelNumber=closureNumber, bondType=nb}:ml)
    where (_, s2, s3) = s =~ "([-=#\\.]{0,1}[/\\]{0,1}[%]{0,1}[0-9]{0,1}){0,1}"::(String, String, String)
          (_, n, _)    = s =~ "([0-9])"::(String, String, String)
          (_, nbs, _)  = s2 =~ "^[-=#\\.]"::(String, String, String)
          nb | nbs == "-" = Single
             | nbs == "=" = Double
             | nbs == "#" = Triple
             | otherwise = Single
          closureNumber = read n::Integer

-- findNextSubSmile
{------------------------------------------------------------------------------}
findNextSubSmile::String -> Int -> (String, String, String)
findNextSubSmile s i | length s >= i && (open == closed)  = ("", s2', s3)
                     | length s > i                       = findNextSubSmile s (i + 1)
                     | length s <= i                      = ("", "", s)  -- Will trigger smiles error message
                     | otherwise                          = ("", "", s)  -- Will trigger smiles error message
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

