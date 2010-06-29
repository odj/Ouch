{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Molecules - a module to manage molecule data types

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


module Ouch.Structure.Molecule
    (
       Molecule(..)
     , addAtom
     , addMarkerToAtomAtIndex
     , addMolecule
     , numberOfAtoms
     , numberOfHeavyAtoms
     , fillMoleculeValence
     , hasHangingClosure
     , molecularWeight
     , exactMass
     , giveMoleculeError
     , moleculeHasError
     , numberOfHydrogenBondDonors
     , numberOfHydrogenBondAcceptors
     , numberOfRings
     , numberOfRotatableBonds
     , updateAtomLabelMarkers
     , molecularFormula
     , connectMoleculesAtIndicesWithBond
     ) where


import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Data.Atom
import Data.Either
import Data.Maybe
import Data.Map as Map
import qualified Data.List as List
import Data.Set as Set
import Control.Applicative



{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}

data Molecule = Small {atomMap::(Map Int Atom)
                     , molMarkerSet::(Set MoleculeMarker)}
                | Markush {molMarkerSet::(Set MoleculeMarker)}
                | Polymer {molMarkerSet::(Set MoleculeMarker)}
                | Biologic {molMarkerSet::(Set MoleculeMarker)}


{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

-- addAtom
-- Default sigma-bonds to top atom on the list
{------------------------------------------------------------------------------}
addAtom :: Molecule -> Atom -> Molecule
addAtom m a = if (moleculeHasError m) then m else case m of
    Small {atomMap=atomM} -> m {atomMap=newAtomMap}
         where (atom, newAtom) = connectAtomsWithBond topAtom a Single
               (k, topAtom) = Map.findMax atomM
               newAtomMap = Map.insert (k+1) newAtom $ Map.insert k atom atomM
    Markush  {}           -> giveMoleculeError m "Cannot add atom to Markush"
    Polymer  {}           -> giveMoleculeError m "Cannot add atom to Markush"
    Biologic {}           -> giveMoleculeError m "Cannot add atom to Markush"




-- findConnectedAtoms
{------------------------------------------------------------------------------}
findConnectedAtoms :: Atom -> [Atom]
findConnectedAtoms a = [a]



-- cyclizePerhapsMolecule
-- Find all matching closure instances and cyclize on matched pairs.
{------------------------------------------------------------------------------}
cyclizeMolecule :: Molecule -> Molecule
cyclizeMolecule m = if (moleculeHasError m) then m else case m of
    Small {} -> case tpl of
                Nothing       -> m
                Just (a1, a2) -> cyclizeMoleculeAtIndexesWithBond m a1 a2 bond
                    where  bond = getMatchingClosureBondType atom1 atom2
                           atom1 = (\(Just a) -> a) $ Map.lookup a1 (atomMap m)
                           atom2 = (\(Just a) -> a) $ Map.lookup a2 (atomMap m)
    Markush  {}    -> giveMoleculeError m "Can't cyclize on a Markush."
    Polymer  {}    -> giveMoleculeError m "Can't cyclize on a Polymer."
    Biologic {}    -> giveMoleculeError m "Can't cyclize on a Biologic."
    where markers = List.map atomMarkerSet $ List.map snd $ Map.toList (atomMap m)
          isClosure mk  = case mk of Closure {} -> True ; _ -> False
          splitMk = List.map fst $ List.map (Set.partition isClosure) markers
          hasPair ms1 ms2 = List.elem True $ (==) <$> labelSet1 <*> labelSet2
          where labelSet1 = List.map labelNumber (Set.toList ms1)
                labelSet2 = List.map labelNumber (Set.toList ms2)
          firstClosure = List.findIndex (/=Set.empty) splitMk
          tpl = case firstClosure of
              Nothing -> Nothing
              Just atom1 -> case secondClosure of
                  Nothing -> Nothing
                  Just atom2 -> Just (atom1, atom2)
                  where secondClosure = List.findIndex (hasPair (splitMk !! atom1)) splitMk2
                        splitMk2 = (take atom1 splitMk) ++ [Set.empty] ++ (drop (atom1+1) splitMk)
--hasHangingClosure
{------------------------------------------------------------------------------}
hasHangingClosure :: Molecule -> Bool
hasHangingClosure m = if (moleculeHasError m) then False else case m of
        Small {}       -> output
        Markush  {}    -> False
        Polymer  {}    -> False
        Biologic {}    -> False
        where markers       = List.map atomMarkerSet $ List.map snd $ Map.toList (atomMap m)
              isClosure mk  = case mk of Closure {} -> True ; _ -> False
              splitMk = List.map fst $ List.map (Set.partition isClosure) markers
              firstClosure = List.findIndex (/=Set.empty) splitMk
              output = case firstClosure of
                  Nothing -> False
                  Just atom1 -> True

-- cyclizePerhapsMoleculeAtIndexesWithBond
{------------------------------------------------------------------------------}

cyclizeMoleculeAtIndexesWithBond :: Molecule -> Int -> Int -> NewBond -> Molecule
cyclizeMoleculeAtIndexesWithBond m i1 i2 b = if (moleculeHasError m) then m else case a1 of
        Nothing     -> (giveMoleculeError m ("Could not connect molecules at index: "
                                     ++ (show i1) ++ " " ++ (show i2)))
        Just atom1  -> case a2 of
            Nothing     -> (giveMoleculeError m ("Could not connect molecules at index: "
                                         ++ (show i1) ++ " " ++ (show i2)))
            Just atom2 -> cyclizeMoleculeAtAtomsWithBond m atom1 atom2 b
                where cyclizeMoleculeAtAtomsWithBond m a1 a2 b
                        | errorTest = cyclizeMolecule (m {atomMap=newMap}) -- !!!Check this!!!
                        | otherwise = giveMoleculeError m "Could not cyclize molecule"
                        where markerLabel = getMatchingClosureNumber atom1 atom2
                              errorTest = case markerLabel of
                                  Nothing -> False
                                  Just label -> True
                               -- This should not evaluate if 'Nothing' because of the guards
                              label = (\(Just l) -> l) markerLabel
                              (newAtom1, newAtom2) = connectAtomsWithBond (removeClosureAtomMarker atom1 label)
                                                     (removeClosureAtomMarker atom2 label) b
                              newMap = Map.insert i1 newAtom1 $ Map.insert i2 newAtom2 (atomMap m)
    where a1 = Map.lookup i1 (atomMap m)
          a2 = Map.lookup i2 (atomMap m)






-- connectPerhapsMoleculesAtIndicesWithBond
-- Takes two 'PerhapsMolecules' and connects them with a 'NewBond' at their
-- respective indices.  Return an error if indices or PerhapsMolecules are
-- invalid.
{------------------------------------------------------------------------------}
connectMoleculesAtIndicesWithBond::Molecule -> Int -> Molecule -> Int -> NewBond -> Molecule
connectMoleculesAtIndicesWithBond m1 i1 m2 i2 b =
    if (moleculeHasError m1) then m1 else
    if (moleculeHasError m2) then m2 else m12
        where m12 | hasClosure = cyclizeMolecule $ connectMolecules m1 i1 m2 i2 b
                  | otherwise  = connectMolecules m1 i1 m2 i2 b
              markers       = List.foldr ((++) . Set.toList . atomMarkerSet) [] $ List.map snd $ Map.toList (atomMap m2)
              isClosure mk  = case mk of Closure {} -> True ; _ -> False
              hasClosure    = List.elem True $ List.map (isClosure) markers
              connectMolecules m1 i1 m2 i2 b
                  | errorTest = Small {atomMap=newMap, molMarkerSet=(Set.union (molMarkerSet m1) (molMarkerSet m2))}
                  | otherwise = giveMoleculeError m1 ("Could not connect molecules at index: "
                                    ++  (show i1) ++ " " ++ (show i2))
                  where a1 = Map.lookup i1 (atomMap m1)
                        a2 = Map.lookup i2 (atomMap m2)
                        -- An error gives a 'False' value
                        errorTest = case a1 of
                            Nothing -> False
                            Just atom1 -> case a2 of
                                Nothing -> False
                                Just atom2 -> True
                        atom1 = (\(Just a) -> a) a1
                        atom2 = (\(Just a) -> a) a2
                        (newAtom1, newAtom2) = connectAtomsWithBond atom1 atom2 b
                        newMap1 = Map.insert i1 newAtom1 (atomMap m1)
                        newMap2 = Map.mapKeysMonotonic (+ orginalLength) $ Map.insert i2 newAtom2 (atomMap m2)
                        orginalLength = Map.size newMap1
                        newMap = Map.union newMap1 newMap2


-- addMolecule
-- Combines two molecules and connects bond-markers if required
-- Otherwise, adds as disconnected structure.  This is not really meant to be accessed
-- directly.
{------------------------------------------------------------------------------}
addMolecule :: Molecule -> Molecule -> Molecule
addMolecule m1 m2 = if (moleculeHasError m1) then m1 else
                    if (moleculeHasError m1) then m2 else m12
    where m12 = cyclizeMolecule (Small {atomMap=newAtomMap, molMarkerSet=newMarkerSet})
          newatomMap = Map.union atomMap1 atomMap2
          atomMap1 = atomMap m1
          atomMap2 = Map.mapKeysMonotonic (+startIndex) $ atomMap m2
          startIndex = Map.size atomMap1
          newAtomMap = Map.union atomMap1 atomMap2
          newMarkerSet = Set.union (molMarkerSet m1) (molMarkerSet m2)


--
{------------------------------------------------------------------------------}
makeMoleculeFromAtom:: Atom -> Molecule
makeMoleculeFromAtom a = Small {atomMap = (Map.singleton 0 a), molMarkerSet=Set.empty}


-- numberOfAtoms
{------------------------------------------------------------------------------}
numberOfAtoms :: Molecule -> Maybe Integer
numberOfAtoms m = if (moleculeHasError m) then Nothing else case m of
    Small {atomMap=atoms}     -> Just $ num atoms
    Markush  {}               -> Nothing
    Polymer  {}               -> Nothing
    Biologic {}               -> Nothing
    where num a = fromIntegral $ Map.size $ Map.filter isElement a


-- numberOfHeavyAtoms
{------------------------------------------------------------------------------}
numberOfHeavyAtoms :: Molecule -> Maybe Integer
numberOfHeavyAtoms m = case m of
    Small {atomMap=atoms} -> Just $ num atoms
    Markush  {}    -> Nothing
    Polymer  {}    -> Nothing
    Biologic {}    -> Nothing
    where heavy a = Map.filter isHeavyAtom a
          num a = fromIntegral $ Map.size $ heavy a


-- fillMoleculeValence
-- Fills valence with hydrogens and lone-pairs to give neutral species
-- If valence is already complete, returns molecule unchanged
-- If cannot fill because of an error, adds error to molecule and gives it back.
{------------------------------------------------------------------------------}
fillMoleculeValence :: Molecule -> Molecule
fillMoleculeValence m = if (moleculeHasError m) then m else case m of
        Small {atomMap=atoms} -> m {atomMap=newAtomMap}
            where newAtomTuple = Map.fold (\a -> (++) [fillValence a []]) [] atoms
                  addedAtomList =  List.foldr (\a -> (++) (snd a)) [] newAtomTuple
                  numberOriginalAtoms = Map.size atoms
                  numberAtomsAdded = length addedAtomList
                  newAtomMap1 = Map.fromList $ zip [0..(numberOriginalAtoms-1)] $ List.map fst newAtomTuple
                  newAtomMap2 = Map.fromList $ zip [numberOriginalAtoms..(numberOriginalAtoms+numberAtomsAdded-1)] addedAtomList
                  newAtomMap  = Map.union newAtomMap1 newAtomMap2
        Markush  {}   -> giveMoleculeError m "Can't fill valence on a Markush."
        Polymer  {}   -> giveMoleculeError m "Can't fill valence on a Polymer."
        Biologic {}   -> giveMoleculeError m "Can't fill valence on a Biologic."

 -- moleculeHasError
 -- Yes if atomMarkerSet contains any 'MolError' value
{------------------------------------------------------------------------------}
moleculeHasError :: Molecule -> Bool
moleculeHasError = Set.member (MError "") . molMarkerSet

 -- giveMoleculeError
 -- Adds an error marker containing a string to the molecule.
{------------------------------------------------------------------------------}
giveMoleculeError :: Molecule -> String -> Molecule
giveMoleculeError m err = markMolecule m $ MError err

 -- getMoleculeError
 -- Adds an error marker containing a string to the molecule.
{------------------------------------------------------------------------------}
getMoleculeError :: Molecule -> Maybe String
getMoleculeError m  = let err =  Set.filter (==(MError "")) $  molMarkerSet m
                      in  if Set.size err == 0 then Nothing
                          else Just $ molMarker $ Set.findMax err



{------------------------------------------------------------------------------}
markMolecule :: Molecule -> MoleculeMarker -> Molecule
markMolecule m mm = newMol
    where newSet = Set.insert mm (molMarkerSet m)
          newMol = m {molMarkerSet=newSet}

-- atomAtIndex
-- Gives you back the atom at index specified.
{------------------------------------------------------------------------------}
atomAtIndex :: Molecule -> Int -> Maybe Atom
atomAtIndex m i = Map.lookup i $ atomMap m

-- addMarkerToAtomAtIndex
-- Adds marker to the specified atom.  If the atom doesn't exist, original
-- molecule is returned with a warnng string added that this operation failed
{------------------------------------------------------------------------------}
addMarkerToAtomAtIndex :: Molecule -> Int -> AtomMarker -> Molecule
addMarkerToAtomAtIndex m i am = if (moleculeHasError m) then m else case atom of
    Nothing -> markMolecule m warning
    Just a  -> m {atomMap = Map.insert i (markAtom a am) (atomMap m)}
    where atom = atomAtIndex m i
          warning = Warning $ "Unable to add marker to atom at position " ++ show i



-- updateAtomLabelMarkers
-- Adds atom label corresponding to the atom map key.  This lets the atom type "know"
-- where it is in the overall Molecule data structure
{------------------------------------------------------------------------------}
updateAtomLabelMarkers :: Molecule -> Molecule
updateAtomLabelMarkers m = if (moleculeHasError m) then m
    else m {atomMap=newMap}
        where jst = (\(Just a) -> a)
              foldFunc = (\k m -> Map.insert k (markAtom (jst $ Map.lookup k m) (Label $ fromIntegral k)) m)
              newMap = List.foldr foldFunc (atomMap m) [0..(Map.size $ atomMap m)-1]

{------------------------------------------------------------------------------}
molecularWeight :: Molecule -> Either String Double
molecularWeight m = if (moleculeHasError m) then (Left "") else case m of
        Small {atomMap=atoms} -> Right $ mw atoms
        Markush  {}   -> Left "No MW for Markush"
        Polymer  {}   -> Left "No MW for Polymer"
        Biologic {}   -> Left "No MW for Biologic"
        where mw a = foldl (+) 0.0 $ List.map atomMW $ Map.fold (\a -> (++) [a]) [] a



-- exactMass
-- This is a distribution of masses according to natural isotopic abundance, normalized to 100
-- This is the normal way in which mass-spec people report mass distribution fractions
{------------------------------------------------------------------------------}
exactMass :: Molecule -> Maybe [(Integer, Double)]
exactMass m = undefined

--
{------------------------------------------------------------------------------}
numberOfHydrogenBondDonors :: Molecule -> Maybe Integer
numberOfHydrogenBondDonors m = undefined

--
{------------------------------------------------------------------------------}
numberOfHydrogenBondAcceptors :: Molecule -> Maybe Integer
numberOfHydrogenBondAcceptors m = undefined

--
{------------------------------------------------------------------------------}
numberOfRings :: Molecule -> Maybe Integer
numberOfRings m = Just (0::Integer)

--
{------------------------------------------------------------------------------}
numberOfRotatableBonds :: Molecule -> Maybe Integer
numberOfRotatableBonds m = Just (0::Integer)

--
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






{------------------------------------------------------------------------------}
{-------------------------------Typeclass Intances-----------------------------}
{------------------------------------------------------------------------------}




instance Show Molecule where
    show m = if (moleculeHasError m) then ("Molecule has error.") else case m of
       Small {atomMap=atoms, molMarkerSet=mm} -> "\nIs a small molecule with formula: "
                ++ (\(Right a) -> a) (molecularFormula $ m) ++ "\n"
                ++ (List.foldr (\b ->  (++) ((show $ fst b) ++ " -- "
                    ++ (show $ snd b))) "" (Map.toList atoms))
                ++ (List.foldr (\b ->  (++) (show b)) "" (Set.toList mm)) ++ "\n"
       Markush  {}   -> "Is a markush"
       Polymer  {}   -> "Is a polymer"
       Biologic {}   -> "Is a Biologic"


{------------------------------------------------------------------------------}
instance Eq Molecule where
    a == b = True --Obviously that needs fixing!

