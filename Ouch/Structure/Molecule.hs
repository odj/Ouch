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
     , addBond
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
import Data.List as List
import Data.Set as Set
import Data.Maybe as Maybe
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

-- checkMolFunSmall
-- Performs check to make sure molecule is valid, otherwise performs the function
-- specified. If check fails, because of type, adds specified string.
{------------------------------------------------------------------------------}
checkMolFunSmall :: Molecule -> Molecule -> String -> Molecule
checkMolFunSmall mIn mOut s = if (moleculeHasError mIn) then mIn else case mIn of
    Small {} -> mOut
    Markush  {}    -> giveMoleculeError mIn "Can't perform on Markush: " ++ s
    Polymer  {}    -> giveMoleculeError mIn "Can't perform on Polymer: " ++ s
    Biologic {}    -> giveMoleculeError mIn "Can't perform on Biologic: " ++ s

checkMolFun :: Molecule -> Molecule -> Molecule
checkMolFun mIn mOut s = if (moleculeHasError mIn) then mIn else mOut



--addAtom
-- Adds atom to top of the atom list with no bonds to the molecule.
{------------------------------------------------------------------------------}
addAtom :: Molecule -> Atom -> Molecule
addAtom m a = checkMolFunSmall m mOut "addAtom"
    where mOut =  m {atomMap = Map.insert (1 + Map.size $ atomMap m) a}

-- numberBondsAtIndex
-- Returns number of covalent connections to other atoms in the molecule
-- graph (i.e. one sigma and two pi bonds count as a 'one' bond)
{------------------------------------------------------------------------------}
numberBondsAtIndex :: Molecule -> Int -> Integer
numberBondsAtIndex a = fromIntegral $ Map.size $ atomBondSet a



-- numberBondsToAtomsAtIndex
{------------------------------------------------------------------------------}
numberBondsToAtomsAtIndex :: Molecule -> Int -> Integer
numberBondsToAtomsAtIndex a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = Set.size b


-- numberBondsToRadicalsAtIndex
{------------------------------------------------------------------------------}
numberBondsToRadicalsAtIndex :: Molecule -> Int -> Integer
numberBondsToRadicalsAtIndex a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = Set.size b

-- numberBondsToHydrogensAtIndex
{------------------------------------------------------------------------------}
numberBondsToHydrogensAtIndex ::  Molecule -> Int -> Integer
numberBondsToHydrogensAtIndex a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = Set.size b


-- numberAromaticBondsToAtomsAtIndex !!! Check this - not right !!!!
{------------------------------------------------------------------------------}
numberAromaticBondsToAtomsAtIndex ::  Molecule -> Int -> Integer
numberAromaticBondsToAtomsAtIndex a = case a of
    Element z n b _ -> nt b
    LonePair b m -> nt b
    Electron b m -> nt b
    Unfilled b m -> nt b
    where nt b = Set.size b

-- numberBondsToHeavyAtomsAtIndices
{------------------------------------------------------------------------------}
numberBondsToHeavyAtomsAtIndices :: Atom -> [Int]
numberBondsToHeavyAtomsAtIndices atom = undefined

-- addNewBond
-- Connects two atom positions with a new bond
{------------------------------------------------------------------------------}
addBond :: Molecule -> Int -> Int -> NewBond -> Molecule
addBond m i1 i2 b  = checkMolFunSmall m mOut "addBond"
    where mOut = if (isNothing a1 || isNothing a2)
                 then giveMoleculeError m "Cannot connect atoms."
                 else m {atomMap=newAtomMap}
          a1 = (\(Just a) -> a) $ Map.lookup i1 $ atomMap m
          a2 = (\(Just a) -> a) $ Map.lookup i2 $ atomMap m
          (s1, s2)  = (atomBondSet a1, atomBondSet a2)
          (ns1, ns2) = case b of
            Single -> (Set.insert (Sigma i2) s1, Set.insert (Sigma i1) s2)
            Double -> (Set.insert (Pi i2) s1, Set.insert (Pi i1) s2)
            Triple -> (Set.insert (PiPi i2) s1, Set.insert (PiPi i1) s2)
            NoBond -> (s1, s2)
          (na1, na2) = (a1 {atomBondSet = ns1}, a2 {atomBondSet = ns2})
          newAtomMap = Map.insert i1 na1 $ Map.insert i2 na2

-- addLonePairAtIndex
-- Create a new lone-pair centered on the atom.
-- Return atom and list containing new lone-pair.
{------------------------------------------------------------------------------}
addLonePairAtIndex :: Molecule -> Int -> Molecule
addLonePairAtIndex m i  = checkMolFunSmall m mOut "addLonePairAtIndex"
    where mOut = addAtom m newAtom
          newAtom = LonePair (Set.singleton $ Single i) Set.empty

--addHydrogenAtIndex
{------------------------------------------------------------------------------}
addHydrogenAtIndex :: Molecule -> Int -> Molecule
addHydrogenAtIndex m i  = checkMolFunSmall m mOut "addHydrogenAtIndex"
    where mOut = addAtom m newAtom
          newAtom = Element 1 0 (Set.singleton $ Single i) Set.empty

--addElectronAtIndex
{------------------------------------------------------------------------------}
addElectronAtIndex :: Molecule -> Int -> Molecule
addElectronAtIndex m i  = checkMolFunSmall m mOut "addElectronAtIndex"
    where mOut = addAtom m newAtom
          newAtom = Electron (Set.singleton $ Single i) Set.empty


-- addUnfilled
-- Create a new unfilled orbital centered on the atom.
-- Return atom and list containing new unfilled orbital.
{------------------------------------------------------------------------------}
addUnfilled :: Molecule -> Int -> Molecule
addUnfilled a = undefined


-- checkValence
-- Verify valence rules are met.  True is what you want.
{------------------------------------------------------------------------------}
checkValenceAtIndex :: Molecule -> Int -> Bool
checkValenceAtIndex m i = True




-- fillValence
-- Populate free valences with hydrogens/lone-pairs.  Return new atom plus an
-- atom list containing all the hydrogens added (adding to second arg).
-- Lone pairs (i.e. for Nitrgen atoms) and empty orbitals (i.e. on Boron)
-- are also added and incuded in the list.
{------------------------------------------------------------------------------}
fillValenceAtIndex :: Molecule -> Int -> Molecule
fillValence m i =
    let val          = valence a
        hBool        = Set.member (ExplicitHydrogen 0) (atomMarkerSet a)
        h            | hBool = numberH $ Set.findMax $ Set.filter (== (ExplicitHydrogen 0)) (atomMarkerSet a)
                     | otherwise = 0
        nba          = (numberOfBondsToAtoms a) + (numberOfBondsToRadicals a)
        nb           = numberOfBonds a
        nbh          = numberOfBondsToHydrogens a
        (aH, asH)    = addHydrogen a as
        (aLP, asLP)  = addLonePair a as
        (aEL, asEL)  = addElectron a as
        nbrB         = (numberOfBondsToRadicals a) == 0
        outputXH   | nbh < h                                       = fillValence aH asH
                   | nb >= ((fst val) + (abs(snd val)))             = (a, as)
                     -- Fill marked aromatics with a radical
                   | nbrB && Set.member AromaticAtom (atomMarkerSet a) = fillValence aEL asEL

                     -- Fill empty valences with radical
                   | nba < fst val                                 = fillValence aEL asEL

                     -- Then, fill lone-pairs if needed
                   | nb < ((fst val) + (abs(snd val)))             = fillValence aLP asLP

                     -- Pattern completion
                   | otherwise = (a, as)
        output     | nb >= ((fst val) + (abs(snd val)))            = (a, as)
                     -- Fill marked aromatics with a radical
                   | nbrB && Set.member AromaticAtom (atomMarkerSet a) = fillValence aEL asEL

                     -- Fill empty valences with hydrogen
                   | nba < fst val                                 = fillValence aH asH

                     -- Then, fill lone-pairs if needed
                   | nb < ((fst val) + (abs(snd val)))             = fillValence aLP asLP

                     -- Pattern completion
                   | otherwise = (a, as)
    in if hBool then outputXH else output


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
          hasPair ms1 ms2 = List.elem True $ (==) <$> labelSet ms1 <*> labelSet ms2
          labelSet ms = List.map labelNumber (Set.toList ms)
          firstClosure = List.findIndex (/=Set.empty) splitMk
          tpl = case firstClosure of
              Nothing -> Nothing
              Just atom1 -> case secondClosure of
                  Nothing -> Nothing
                  Just atom2 -> Just (atom1, atom2)
                  where secondClosure = List.findIndex
                                       (hasPair (splitMk !! atom1)) splitMk2
                        splitMk2 = (take atom1 splitMk)
                                 ++ [Set.empty]
                                 ++ (drop (atom1+1) splitMk)


--hasHangingClosure
{------------------------------------------------------------------------------}
hasHangingClosure :: Molecule -> Bool
hasHangingClosure m = if (moleculeHasError m) then False else case m of
    Small {}  -> output
    _         -> False
    where markers = List.map atomMarkerSet
                  $ List.map snd
                  $ Map.toList (atomMap m)
          isClosure mk  = case mk of Closure {} -> True ; _ -> False
          splitMk = List.map fst $ List.map (Set.partition isClosure) markers
          firstClosure = List.findIndex (/=Set.empty) splitMk
          output = case firstClosure of
              Nothing    -> False
              Just atom1 -> True


--cyclizeMoleculeAtIndexesWithBond
{------------------------------------------------------------------------------}
cyclizeMoleculeAtIndexesWithBond :: Molecule -> Int -> Int -> NewBond -> Molecule
cyclizeMoleculeAtIndexesWithBond m i1 i2 b =
    let a1 = Map.lookup i1 (atomMap m)
        a2 = Map.lookup i2 (atomMap m) in
    if (moleculeHasError m) then m else case a1 of
        Nothing     -> (giveMoleculeError m
                       ("Could not connect molecules at index: "
                    ++ (show i1) ++ " " ++ (show i2)))
        Just atom1  -> case a2 of
            Nothing -> (giveMoleculeError m
                       ("Could not connect molecules at index: "
                    ++ (show i1) ++ " " ++ (show i2)))
            Just atom2 -> cyclizeMoleculeAtAtomsWithBond m atom1 atom2 b
                where cyclizeMoleculeAtAtomsWithBond m a1 a2 b
                        | errorTest = cyclizeMolecule (m {atomMap=newMap})
                        | otherwise = giveMoleculeError m
                          "Could not cyclize molecule"
                        where markerLabel = getMatchingClosureNumber atom1 atom2
                              errorTest = case markerLabel of
                                  Nothing -> False
                                  Just label -> True
                              label = (\(Just l) -> l) markerLabel
                              (newAtom1, newAtom2) = connectAtomsWithBond (removeClosureAtomMarker atom1 label)
                                                     (removeClosureAtomMarker atom2 label) b
                              newMap = Map.insert i1 newAtom1 $ Map.insert i2 newAtom2 (atomMap m)


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



--makeMoleculeFromAtom
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

--markMolecule
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


--molecularWeight
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


--Function Stubs

--numberOfHydrogenBondDonors
{------------------------------------------------------------------------------}
numberOfHydrogenBondDonors :: Molecule -> Maybe Integer
numberOfHydrogenBondDonors m = undefined

--numberOfHydrogenBondAcceptors
{------------------------------------------------------------------------------}
numberOfHydrogenBondAcceptors :: Molecule -> Maybe Integer
numberOfHydrogenBondAcceptors m = undefined

--numberOfRings
{------------------------------------------------------------------------------}
numberOfRings :: Molecule -> Maybe Integer
numberOfRings m = Just (0::Integer)

--numberOfRotatableBonds
{------------------------------------------------------------------------------}
numberOfRotatableBonds :: Molecule -> Maybe Integer
numberOfRotatableBonds m = Just (0::Integer)

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

