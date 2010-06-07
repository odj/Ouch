-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-- Project !Ouch
-- No license selected yet-- project still under development

-- Orion D. Jankowski
-- 2010-May-24
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

module Ouch.Structure.Molecule 
    (
       Molecule(..)
     , PerhapsMolecule(..)
     , addAtom
     , addMolecule
     , numberOfAtoms
     , numberOfHeavyAtoms
     , fillMoleculeValence
     , molecularWeight
     , exactMass
     , numberOfHydrogenBondDonors
     , numberOfHydrogenBondAcceptors
     , numberOfRings
     , numberOfRotatableBonds
     , molecularFormula
     , connectMoleculesAtIndicesWithBond
     ) where

import Ouch.Structure.Atom
import Ouch.Data.Atom
import Data.Either
import Data.Maybe
import Data.Map as Map
import Data.List as List

data Molecule = Small {atomList::[Atom]} 
                | Markush {}
                | Polymer {}
                | Biologic {}

-- Use this data structure when a molecule is still being constructed and
-- might encounter an error later on.  When an error is encountered, stop construction
-- and propogate error message to a function that actually cares.  Used mostly in
-- INPUT module
type PerhapsMolecule =  (Either String Molecule) 

-- addAtom
-- Default sigma-bonds to top atom on the list
addAtom :: Molecule -> Atom -> Molecule
addAtom m a = case m of
    Small {atomList=(atom:atoms)} -> Small (atoms' ++ [atom'] ++ atoms)
         where (atom', atoms') = sigmaBondToAtom atom a
    Markush            -> m
    Polymer            -> m
    Biologic           -> m

connectMoleculesAtIndicesWithBond :: Molecule -> Int -> Molecule -> Int -> NewBond -> PerhapsMolecule
connectMoleculesAtIndicesWithBond m1 i1 m2 i2 b = Right $ Small {atomList=(a1b ++ [newAtom1] ++ a1e ++ a2b ++ [newAtom2] ++ a2e)}
    where atom1 = (atomList m1) !! i1
          atom2 = (atomList m2) !! i2
          (a1b, a1e) = (take i1 (atomList m1), drop (i1+1) (atomList m1))
          (a2b, a2e) = (take i2 (atomList m2), drop (i2+1) (atomList m2))
          (newAtom1, newAtom2) = connectAtomsWithBond atom1 atom2 b 
          
-- addMolecule
-- Combines two molecules and connects bond-markers if required
-- Otherwise, adds as disconnected structure
addMolecule :: Molecule -> Molecule -> Molecule
addMolecule m1 m2 = undefined
   
makeMoleculeFromAtom:: Atom -> Molecule
makeMoleculeFromAtom a = undefined       

numberOfAtoms :: Molecule -> Maybe Integer
numberOfAtoms m = case m of
    Small {atomList=atoms} -> Just $ num atoms
    Markush                -> Nothing
    Polymer                -> Nothing
    Biologic               -> Nothing
    where elements a = List.filter isElement a
          num a = fromIntegral $ length $ elements a

numberOfHeavyAtoms :: Molecule -> Maybe Integer
numberOfHeavyAtoms m = case m of
    Small {atomList=atoms} -> Just $ num atoms
    Markush     -> Nothing
    Polymer     -> Nothing
    Biologic    -> Nothing
    where heavy a = List.filter isHeavyAtom a
          num a = fromIntegral $ length $ heavy a

-- Fills valence with hydrogens and lone-pairs to give neutral species
-- If valence is already complete, returns molecule unchanged 
-- If cannot fill because of an error, returns error string.
fillMoleculeValence :: PerhapsMolecule -> PerhapsMolecule
fillMoleculeValence pm = case pm of
    Left {} -> pm
    Right m -> case m of
        Small {atomList=atoms} -> Right $ newMolecule atoms
        Markush     -> Left "Can't fill valence on a Markush."
        Polymer     -> Left "Can't fill valence on a Polymer."
        Biologic    -> Left "Can't fill valence on a Biologic."
        where newAtomTuple a = List.map (\a -> fillValence a []) a
              newAtomsList a = (List.map fst $ newAtomTuple a) ++ (foldl (++) [] (List.map snd $ newAtomTuple a))
              newMolecule a = Small $ newAtomsList a


molecularWeight :: PerhapsMolecule -> Either String Double
molecularWeight pm = case pm of
    Left m  -> Left m
    Right m -> case m of
        Small {atomList=atoms} -> Right $ mw atoms
        Markush     -> Left "No MW for Markush"
        Polymer     -> Left "No MW for Polymer"
        Biologic    -> Left "No MW for Biologic"
        where mw a = foldl (+) 0.0 $ List.map atomMW a

-- exactMass
-- This is a distribution of masses according to natural isotopic abundance, normalized to 100
-- This is the normal way in which mass-spec people report mass distribution fractions
exactMass :: PerhapsMolecule -> Either String [(Integer, Double)]
exactMass m = undefined

numberOfHydrogenBondDonors :: PerhapsMolecule -> Either String Integer
numberOfHydrogenBondDonors m = Right (0::Integer)

numberOfHydrogenBondAcceptors :: PerhapsMolecule -> Either String Integer
numberOfHydrogenBondAcceptors m = Right (0::Integer)

numberOfRings :: PerhapsMolecule -> Either String Integer
numberOfRings m = Right (0::Integer)

numberOfRotatableBonds :: PerhapsMolecule -> Either String Integer
numberOfRotatableBonds m = Right (0::Integer)

molecularFormula :: PerhapsMolecule -> Either String String
molecularFormula pm = case pm of
    Left s      -> Left s
    Right mol   -> case mol of 
        Small at    -> Right molFm 
            where startMap = Map.empty
                  endMap = List.foldr (updateMap) startMap at
                  -- Use foldr to accumulate and count atoms
                  updateMap a m | Map.notMember (atomicSymbolForAtom a) m   = Map.insert (atomicSymbolForAtom a)  1 m
                             | otherwise                                    = Map.adjust (+ 1) (atomicSymbolForAtom a) m           
                  -- Convert the map to a list of just the elements present, and in IUPAC order
                  finalList = catMaybes $ List.map (\e -> lookupPair e endMap) molecularFormulaElements
                  --  Build the final output string from the map
                  molFm = List.foldr (\(e,n) -> if n>1 then ((e ++  (show n))++) else (e ++ ))  "" finalList
                  -- simple little utility function which, strangely, is not already defined in Data.Map
                  lookupPair k m = case v of
                      Just val -> Just (k, val)
                      Nothing -> Nothing
                      where v = Map.lookup k m
        Markush     -> Left "No molecular formula defined for Markush"
        Polymer     -> Left "No molecular formula defined for Polymer"
        Biologic    -> Left "No molecular formula defined for Biologic"





instance Show Molecule where
    show m = case m of
       Small {atomList=atoms} -> "Is a small molecule with atoms: " ++ show atoms
       Markush     -> "Is a markush"
       Polymer     -> "Is a polymer"
       Biologic    -> "Is a Biologic"


instance Eq Molecule where
    a == b = True --Obviously that needs fixing!

