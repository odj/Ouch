module Ouch.Input.Smiles (
     nextChoppedSmile
   , makeMoleculeFromSmiles
   , makeAtomMoleculeFromChop
   , findNextSubSmile
   , parseSmiles
   , chop
   ) where
       
import Ouch.Structure.Atom
import Ouch.Structure.Molecule
import Text.Regex.TDFA ((=~))
import Data.Maybe
import Data.Char
import Data.Set as Set
import Data.List as List

-- NewBond
-- We package each smile-parsing step as a "NewBond" data type to make
-- it convenient for different functions to determine what they should do without
-- needing to re-parse parts of the string every single time.  In Haskell,
-- this probably doesn't make things more efficient, but it helps for DRY style compliance.
data ChoppedSmile = Smile {smile::String, smiles::String, newBond::NewBond, mark::(Set Marker)}
                  | SubSmile {smile::String, smiles::String, newBond::NewBond, mark::(Set Marker)}  --Substructures don't need a marker set
                  | SmilesError {smile::String, smiles::String, newBond::NewBond, mark::(Set Marker)}  --
                  deriving (Show, Eq)

-- We do not include aromatic bonds here because they are not specified by 
-- labels, but by atom case or pure inference.


-- makeMoleculeFromSmiles
-- Right if successful
-- Left if error somewhere, with brief error description
-- This functions kicks off the recursive interpretations
makeMoleculeFromSmiles::String -> PartialMolecule
makeMoleculeFromSmiles smi = case chop of 
    Smile {}        -> growPartialMoleculeAtIndexWithString newAtom 0 nextSmile 
    -- A well-formed smile should never start with a substructure, but if it does, we connect the next 
    -- atom or substructure to it's beginning (which may or may not be permitted chemically)
    SubSmile {}     -> growPartialMoleculeAtIndexWithString newSubstructure 0 nextSmile  
    SmilesError {}  -> (Left $ "Error trying to parse the Smiles string: " ++ (smi))
    where chop = nextChoppedSmile smi
          thisSmile = smile chop
          nextSmile = smiles chop
          newAtom = makeAtomMoleculeFromChop chop
          newSubstructure = makeMoleculeFromSmiles thisSmile
          
connectPartialMoleculesAtIndicesWithBond::PartialMolecule -> Int -> PartialMolecule -> Int -> NewBond -> PartialMolecule
connectPartialMoleculesAtIndicesWithBond pm1 i1 pm2 i2 b =
    case pm1 of 
        Left {} -> pm1
        Right m1 -> case pm2 of
            Left {} -> pm2
            Right m2  | errorTest == False    -> (Left ("Could not connect molecules at index: " ++ (show i1) ++ " " ++ (show i2)))
                      | otherwise             -> connectMoleculesAtIndicesWithBond m1 i1 m2 i2 b
                      where errorTest = (test1 && test2)
                            n1 = numberOfAtoms m1
                            n2 = numberOfAtoms m2
                            test1 = case n1 of
                                Just n -> (fromIntegral n > i1) && (i1 >= 0)
                                Nothing -> True
                            test2 = case n2 of
                                Just n -> (fromIntegral n > i2) && (i2 >= 0)
                                Nothing -> False

-- growPartialMoleculeWithString
-- Adds continuation of smiles string at the end of the molecule being built
growPartialMoleculeAtIndexWithString :: PartialMolecule -> Int -> String -> PartialMolecule
growPartialMoleculeAtIndexWithString pm i smi 
    | smi == ""  = pm
    | otherwise  = case pm of
    Right {} -> case chop of
            Smile {}        -> growPartialMoleculeAtIndexWithString newMolecule1 newIndex nextSmile 
            SubSmile {}     -> growPartialMoleculeAtIndexWithString newMolecule2 i nextSmile
            SmilesError {}  -> (Left $ "Error trying to grow from the Smiles string: " ++ smi)
            where chop = nextChoppedSmile smi
                  nextSmile = smiles chop
                  newAtom = makeAtomMoleculeFromChop chop
                  newSubStructure = makeMoleculeFromSmiles (smile chop)
                  newMolecule1 = connectPartialMoleculesAtIndicesWithBond pm i newAtom 0 (newBond chop)
                  newMolecule2 = connectPartialMoleculesAtIndicesWithBond pm i newSubStructure 0 (newBond chop)
                  -- We just made this thing, so there shouldn't be any errors, right?
                  newIndex = case newMolecule1 of Right m1 -> (fromIntegral $ fromJust $ numberOfAtoms m1)-1          
    Left {} -> pm    




-- makeMoleculeFromSmiles
-- Right if successful
-- Left if error somewhere, with brief error description
makeAtomMoleculeFromChop::ChoppedSmile -> PartialMolecule
makeAtomMoleculeFromChop nb    | a == ""        =  (Left "ERROR: Tried to make atom from empty string.")
                               | a == "C"       =  (Right $ Small [Element 6 0 [] markSet1])
                               | a == "N"       =  (Right $ Small [Element 7 0 [] markSet1])
                               | a == "O"       =  (Right $ Small [Element 8 0 [] markSet1])
                               | a == "H"       =  (Right $ Small [Element 1 0 [] markSet1])
                               
                               | a == "P"       =  (Right $ Small [Element 15 0 [] markSet1])
                               | a == "S"       =  (Right $ Small [Element 16 0 [] markSet1])
                               | a == "F"       =  (Right $ Small [Element 9 0 [] markSet1])
                               | a == "B"       =  (Right $ Small [Element 5 0 [] markSet1])
                               | a == "BR"       =  (Right $ Small [Element 35 0 [] markSet1])
                               | a == "CL"       =  (Right $ Small [Element 17 0 [] markSet1])
                               | a == "I"       =  (Right $ Small [Element 53 0 [] markSet1])
                               | a == "*"       =  (Right $ Small [Unspecified [] markSet1]) -- Wildcard Atom
                               | otherwise      =  (Left  $ "ERROR: Atom not recognized for symbol: " ++ a)
                               where markSet1 = Set.singleton (if (isUpper $ head (smile nb)) then AromaticAtom else Null)
                                     a = [toUpper c | c <- (smile nb)] 



-- nextSmilesSubstring
-- Lots, lots, lots!!!!! more to fill in here
parseSmiles :: String -> (String, String, String)
parseSmiles s = s =~ pat::(String, String, String)
    -- Need a more general format to recognize two-letter elements (or maybe just enumerate them?) 
    where pat = List.foldr (++) "" patList 
          patList = [ "("
                    , "^[-=#\\.]*([A-Za-z]|Br|Cl|Si|Sn|Li|Na|Cs){1}([-=#\\.]{0,1}[0-9])*[@]*"  -- Search for first atom + bond/marker
                    , "|^[\\(].*" -- return the whole string for anything that STARTS with open parens
                    , "|(^[\\[].*(\\]))"  -- return next atom segment within square brackets
                    , ")"
                    ]
                          

-- nextSmilesSubstring
-- Function strips atom/bond info into a new type. 
-- String residual in unmodified (smiles), while atom string (smile) is suitable
-- for parsing by  
nextChoppedSmile :: String -> ChoppedSmile 
nextChoppedSmile s
   | s2 == "" || s1 /= ""   = SmilesError {smile=s1, smiles=s, newBond=NoBond, mark=(Set.singleton $ Comment s2)}
   | (head s2) == '('       = SubSmile {smile=bb2, smiles=ss3, newBond=nb2, mark=(Set.singleton $ Comment ss2)}
   | otherwise              = Smile {smile=a2, smiles=s3, newBond=nb, mark=Set.union markerSet (Set.singleton $ Comment s2)}
   where (s1, s2, s3)       = parseSmiles s                                 -- Get initial parse
         (b1, b2, b3)       = s2 =~ "(^[-=#\\.])"::(String, String, String) -- Get bond info for single atoms
         (l1, l2, l3)       = s2 =~ "([0-9]+)"::(String, String, String)    -- Get atom marker info
         (a1, a2, a3)       = s2 =~ "([A-Za-z]+)"::(String, String, String) -- Atom only, remove bond comments, etc
         (ss1, ss2, ss3)    = findNextSubSmile s 1
         (bb1, bb2, bb3)    = ss2 =~ "([A-Za-z]+.*)"::(String, String, String) -- Get bond info for subsmiles
         nb  | b2 == "."              = NoBond
             | b2 == "" || b2 == "-"  = Single
             | b2 == "="              = Double
             | b2 == "#"              = Triple
             | otherwise              = Single
         -- need to process a different string to determine bond type if parse gives a subsmile
         nb2 | bb1 == "."             = NoBond
             | bb1 == "" || bb2 == "-"= Single
             | bb1 == "="             = Double
             | bb1 == "#"             = Triple
             | otherwise              = Single
         -- This is a very simple parse of ring closure markers.  
         -- Does not accomodate "%" notation (yet)
         markerNumbers = List.map (\a -> read a::Integer) $ foldr (\acc x -> [acc] : x) [] l2
         markerSet = Set.fromList $ List.map Label markerNumbers



findNextSubSmile::String -> Int -> (String, String, String)
findNextSubSmile s i | (length s) >= i && (open == closed)  = ("", s2', s3)
                     | (length s) > i                       = findNextSubSmile s (i + 1)
                     | (length s) <= i                      = ("", "", s)  -- Will trigger smiles error message
                     | otherwise                            = ("", "", s)  -- Will trigger smiles error message
                     where s2 = take i s
                           s3 = drop i s
                           s2' | (i > 2) = tail $ init s2 | otherwise = ""
                           open = sum [1 | c <- s2, c=='(']
                           closed = sum [1 | c <- s2, c==')']



-- GHCI debugging function
chop :: String -> String
chop [] = ""
chop a = (show nb) ++ "\n" ++ chop (smiles nb)
   where nb = nextChoppedSmile a



-----------
-- Some old stuff from the graveyard
{-

parseSmiles :: String -> (String, String, String)
parseSmiles s = s =~ (a ++ b ++ c ++ d ++ e)::(String, String, String)
    where (a, b, c, d, e) = findRecursionLevelPattern startPattern s 1
          startPattern =  ( "([=#]*.[0-9]*[@]*|^[\\(][^\\)^\\(]*"
                          , "([\\(][^\\)]"
                          , "*"
                          , "[^\\(]*[\\)])"
                          , "[^\\)]*[\\)])"
                          )




findRecursionLevelPattern::(String, String, String, String, String) -> String -> Int -> (String, String, String, String, String) 
findRecursionLevelPattern (fr, inf, c, inb, b) s i
        | open == closed      = (fr, inf, c, inb, b)
        | i > (length s)      = (fr, inf, c, inb, b)   -- Not a well formed SMILES.  Should throw error.
        | otherwise           = findRecursionLevelPattern (fr ++ inf, innerFront, center, innerBack, inb ++ b) s (i+1)
            where innerFront = "([\\(][^\\)]"
                  center     = "*"
                  innerBack  = "[^\\)]*[\\)])*"
                  (s1, s2, s3) = s =~ (fr ++ inf ++ c ++ inb ++ b)::(String, String, String)
                  open = sum [1 | c <- s2, c=='(']
                  closed = sum [1 | c <- s2, c==')']

-}
