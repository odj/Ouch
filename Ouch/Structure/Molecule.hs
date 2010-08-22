{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Molecule - a module to manage molecule data types

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
     , atomAtIndex
     , addHydrogenAtIndex
     , addLonePairAtIndex
     , addElectronAtIndex
     , addMolMarker
     , addProperty
     , addPropertyFromFunction
     , bondTargetSetForIndex
     , setAtom
     , getAtomAtIndex
     , addMarkerToAtomAtIndex
     , addMolecule
     , numberOfAtoms
     , makeMoleculeFromAtom
     , numberOfHeavyAtoms
     , fillMoleculeValence
     , fillValenceAtIndex
     , hasHangingClosure
     , cyclizeMolecule
     , molecularWeight
     , exactMass
     , giveMoleculeError
     , moleculeHasError
     , numberBondsAtIndex
     , occupiedValenceAtIndex
     , numberBondsToAtomsAtIndex
     , numberBondsToRadicalsAtIndex
     , numberBondsToHydrogensAtIndex
     , molecularFormula
     , connectMoleculesAtIndicesWithBond
     , incrementAtomMap
     ) where


import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Data.Atom
import Ouch.Property.Property

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

data Molecule =   Molecule {atomMap::(Map Int Atom)
                          , molMarkerSet::(Set MoleculeMarker)
                          , molPropertyMap::(Map String Property)}


{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}


checkMolFun :: Molecule -> Molecule -> Molecule
checkMolFun mIn mOut = if (moleculeHasError mIn) then mIn else mOut

-- getAtomAtIndex
{------------------------------------------------------------------------------}
getAtomAtIndex :: Molecule -> Int -> Maybe Atom
getAtomAtIndex m i = Map.lookup i $ atomMap m

-- setAtom
-- Inserts atom into given molecule at whatever index the atom thinks it should
-- be at.  If the index is invalid or non-existant, then it adds to the next
-- atom index.
{------------------------------------------------------------------------------}
setAtom :: Atom -> Molecule -> Molecule
setAtom a m  = checkMolFun m mOut
  where mOut | (isJust atomIndex) /= True  =  addAtom a m
             | fromJust atomIndex > (Map.size $ atomMap m) = withWarning
             | otherwise = output
        atomIndex = getIndexForAtom a
        withWarning = addAtom a (giveMoleculeWarning m "Tried to set atom with invalid index")
        output = case atomIndex of
          Just i  -> m {atomMap = Map.insert i a $ atomMap m}
          Nothing -> addAtom a m


--addAtom
-- Adds atom to top of the atom list with no bonds to the molecule.
{------------------------------------------------------------------------------}
addAtom :: Atom -> Molecule -> Molecule
addAtom a m = checkMolFun m mOut
    where mOut =  m {atomMap = Map.insert atomNumber newAtom $ atomMap m}
          atomNumber = Map.size $ atomMap m
          newAtom = a {atomMarkerSet=(Set.insert (Label atomNumber) $ atomMarkerSet a) }


addMolMarker :: Molecule -> MoleculeMarker -> Molecule
addMolMarker m mm = m {molMarkerSet = Set.insert mm $ molMarkerSet m}

addProperty :: Molecule -> Property -> Molecule
addProperty m p = m {molPropertyMap = Map.insert (propertyKey p) p $ molPropertyMap m}

addPropertyFromFunction :: Molecule -> (Molecule -> Property) -> Molecule
addPropertyFromFunction m f = addProperty m $ f m


occupiedValenceAtIndex :: Molecule -> Int -> Integer
occupiedValenceAtIndex m i = case atomAtIndex m i of
  Just a  -> occupiedValence a
  Nothing -> 0


bondTargetSetForIndex :: Molecule -> Int -> (Set Atom)
bondTargetSetForIndex m i = let atom = atomAtIndex m i in
  case atom of
  Nothing -> Set.empty
  Just a -> targets
    where bonds = atomBondSet a
          labels = Set.map bondsTo bonds
          targets = Set.map fromJust $ Set.filter isJust $ Set.map (atomAtIndex m) labels


numberOfBondToAtomWithProperty :: Molecule -> Int -> (Atom -> Bool) -> Integer
numberOfBondToAtomWithProperty m i f = let atoms = bondTargetSetForIndex m i in
  fromIntegral $ Set.size $ Set.filter f atoms



-- numberBondsAtIndex
-- Returns number of covalent connections to other atoms in the molecule
-- graph (i.e. one sigma and two pi bonds count as a 'one' bond)
{------------------------------------------------------------------------------}
numberBondsAtIndex :: Molecule -> Int -> Integer
numberBondsAtIndex m i = fromIntegral $ Set.size $ bondTargetSetForIndex m i


-- numberBondsToAtomsAtIndex
{------------------------------------------------------------------------------}
numberBondsToAtomsAtIndex :: Molecule -> Int -> Integer
numberBondsToAtomsAtIndex m i = numberOfBondToAtomWithProperty m i isElement


-- numberBondsToRadicalsAtIndex
{------------------------------------------------------------------------------}
numberBondsToRadicalsAtIndex :: Molecule -> Int -> Integer
numberBondsToRadicalsAtIndex m i = numberOfBondToAtomWithProperty m i isElectron


-- numberBondsToHydrogensAtIndex
{------------------------------------------------------------------------------}
numberBondsToHydrogensAtIndex ::  Molecule -> Int -> Integer
numberBondsToHydrogensAtIndex m i = numberOfBondToAtomWithProperty m i isHydrogen


-- numberAromaticBondsToAtomsAtIndex !!! Check this - not right !!!!
{------------------------------------------------------------------------------}
numberAromaticBondsToAtomsAtIndex ::  Molecule -> Int -> Integer
numberAromaticBondsToAtomsAtIndex m i = numberOfBondToAtomWithProperty m i
                                        (\a -> hasMarker a AromaticAtom)

-- numberBondsToHeavyAtomsAtIndices
{------------------------------------------------------------------------------}
numberBondsToHeavyAtomsAtIndex :: Molecule -> Int -> Integer
numberBondsToHeavyAtomsAtIndex m i  = numberOfBondToAtomWithProperty m i isHeavyAtom


-- numberBondsToLonePairAtIndices
{------------------------------------------------------------------------------}
numberBondsToLonePairAtIndex :: Molecule -> Int -> Integer
numberBondsToLonePairAtIndex m i  = numberOfBondToAtomWithProperty m i isLonePair



-- addBond
-- Connects two atom positions with a new bond
{------------------------------------------------------------------------------}
addBond :: Molecule -> Int -> Int -> NewBond -> Molecule
addBond m i1 i2 b  = checkMolFun m mOut
    where mOut | (isNothing a1 || isNothing a2) =
                    giveMoleculeError m $ "Cannot connect atoms at positions: "
                    ++ (show i1) ++ " " ++ (show i2)
               | otherwise =  setAtom na1 $ setAtom na2 m
          a1 = getAtomAtIndex m i1
          a2 = getAtomAtIndex m i2
          atom1 = fromJust a1
          atom2 = fromJust a2
          (bs1, bs2)  = (atomBondSet atom1, atomBondSet atom2)
          (nbs1, nbs2) = case b of
            Single -> (Set.insert (Sigma i2) bs1, Set.insert (Sigma i1) bs2)
            Double -> (Set.insert (Pi i2) bs1, Set.insert (Pi i1) bs2)
            Triple -> (Set.insert (PiPi i2) bs1, Set.insert (PiPi i1) bs2)
            NoBond -> (bs1, bs2)
          (na1, na2) = (atom1 {atomBondSet = nbs1}, atom2 {atomBondSet = nbs2})
          --newAtomMap = Map.insert i1 na1 $ Map.insert i2 na2 $ atomMap m

-- addLonePairAtIndex
-- Create a new lone-pair centered on the atom.
{------------------------------------------------------------------------------}
addLonePairAtIndex :: Molecule -> Int -> Molecule
addLonePairAtIndex m i  = checkMolFun m mOut "addLonePairAtIndex"
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = LonePair Set.empty Set.empty

--addHydrogenAtIndex
{------------------------------------------------------------------------------}
addHydrogenAtIndex :: Molecule -> Int -> Molecule
addHydrogenAtIndex m i  = checkMolFun m mOut "addHydrogenAtIndex"
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = (Element 1 0) Set.empty Set.empty

--addElectronAtIndex
{------------------------------------------------------------------------------}
addElectronAtIndex :: Molecule -> Int -> Molecule
addElectronAtIndex m i  = checkMolFun m mOut "addElectronAtIndex"
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = Electron Set.empty Set.empty


-- addUnfilled
-- Create a new unfilled orbital centered on the atom.
-- Return atom and list containing new unfilled orbital.
{------------------------------------------------------------------------------}
addUnfilled :: Molecule -> Int -> Molecule
addUnfilled m i  = checkMolFun m mOut "addUnfilled"
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = Unfilled Set.empty Set.empty


-- checkValence
-- Verify valence rules are met.  True is what you want.
{------------------------------------------------------------------------------}
checkValenceAtIndex :: Molecule -> Int -> Bool
checkValenceAtIndex m i = True




-- fillValenceAtIndex
-- Populate free valences with hydrogens/lone-pairs.  Return new atom plus an
-- atom list containing all the hydrogens added (adding to second arg).
-- Lone pairs (i.e. for Nitrgen atoms) and empty orbitals (i.e. on Boron)
-- are also added and incuded in the list.
{------------------------------------------------------------------------------}
fillValenceAtIndex :: Molecule -> Int -> Molecule
fillValenceAtIndex m i = checkMolFun m mOut "fillValenceAtIndex"
  where atom  = getAtomAtIndex m i
        mOut  = case atom of
          Just _  -> output
          Nothing -> giveMoleculeError m $ "Cannot fill valence at position: " ++ (show i)
        val   = valence a
        a     = fromJust atom
        hBool = Set.member (ExplicitHydrogen 0) (atomMarkerSet a)
        h | hBool = numberH $ Set.findMax $ Set.filter (== (ExplicitHydrogen 0)) (atomMarkerSet a)
          | otherwise = 0
        nba   = (occupiedValenceAtIndex m i) - (numberBondsToLonePairAtIndex m i)
        nb    = occupiedValenceAtIndex m i
        nbh   = numberBondsToHydrogensAtIndex m i
        nbrB  = (numberBondsToRadicalsAtIndex m i) == 0
        mH    = addHydrogenAtIndex m i
        mLP   = addLonePairAtIndex m i
        mEL   = addElectronAtIndex m i

        outputXH   | nbh < h                            = fillValenceAtIndex mH i
                   | nb >= ((fst val) + (abs(snd val))) = m
                     -- Fill marked aromatics with a radical
                   | nbrB && Set.member AromaticAtom (atomMarkerSet a) = fillValenceAtIndex mEL i

                     -- Fill empty valences with radical
                   | nba < fst val                      = fillValenceAtIndex mEL i

                     -- Then, fill lone-pairs if needed
                   | nb < ((fst val) + (abs(snd val)))  = fillValenceAtIndex mLP i

                     -- Pattern completion
                   | otherwise                          = m

        outputH    | nb >= ((fst val) + (abs(snd val))) = m
                     -- Fill marked aromatics with a radical
                   | nbrB && Set.member AromaticAtom (atomMarkerSet a) = fillValenceAtIndex mEL i

                     -- Fill empty valences with hydrogen
                   | nba < fst val                      = fillValenceAtIndex mH i

                     -- Then, fill lone-pairs if needed
                   | nb < ((fst val) + (abs(snd val)))  = fillValenceAtIndex mLP i

                     -- Pattern completion
                   | otherwise = m

        output     | hBool      = outputXH
                   | otherwise  = outputH


-- cyclizeMolecule
-- Find all matching closure instances and cyclize on matched pairs.
{------------------------------------------------------------------------------}
cyclizeMolecule :: Molecule -> Molecule
cyclizeMolecule m = checkMolFunSmall m mOut "cyclizeMolecule"
  where mOut = case tpl of
                Nothing       -> m
                Just (a1, a2) -> cyclizeMoleculeAtIndexesWithBond m a1 a2
                      $ getMatchingClosureBondTypeAtIndices m a1 a2
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
                          where secondClosure = List.findIndex (hasPair (splitMk !! atom1)) splitMk2
                                splitMk2 = (take atom1 splitMk) ++ [Set.empty] ++ (drop (atom1+1) splitMk)



getMatchingClosureBondTypeAtIndices :: Molecule -> Int -> Int -> NewBond
getMatchingClosureBondTypeAtIndices m i1 i2 = getMatchingClosureBondType a1 a2
  where a1 = fromJust $ getAtomAtIndex m i1
        a2 = fromJust $ getAtomAtIndex m i2




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
            Just atom2  | hasClosure mOut -> cyclizeMolecule mOut
                        | otherwise       -> mOut
                where mOut
                        | errorTest = addBond (m {atomMap=newMap}) i1 i2 b
                        | otherwise = giveMoleculeError m
                          "Could not cyclize molecule"
                        where markerLabel = getMatchingClosureNumber atom1 atom2
                              errorTest = case markerLabel of
                                  Nothing -> False
                                  Just label -> True
                              label = fromJust markerLabel
                              newAtom1 = removeClosureAtomMarker atom1 label
                              newAtom2 = removeClosureAtomMarker atom2 label
                              newMap = Map.insert i1 newAtom1 $ Map.insert i2 newAtom2 (atomMap m)


-- connectMoleculesAtIndicesWithBond
-- Takes two molecule and connects them with a 'NewBond' at their
-- respective indices.  Return an error if indices or molecules are
-- invalid.
{------------------------------------------------------------------------------}
connectMoleculesAtIndicesWithBond::Molecule -> Int -> Molecule -> Int -> NewBond -> Molecule
connectMoleculesAtIndicesWithBond m1 i1 m2 i2 b = checkMolFunSmall m1 mOut "connectMoleculesAtIndicesWithBond"
  where mOut = connectMolecules m1 i1 m2 i2 b
        connectMolecules m1 i1 m2 i2 b | errorTest = giveMoleculeError m1 "Could not connect molecules, invalid index"
                                       | otherwise = output
        a1 = getAtomAtIndex m1 i1
        a2 = getAtomAtIndex m2 i2
        errorTest = (isJust a1 && isJust a1) /= True
        m1Length = Map.size $ atomMap m1
        output = addBond (addMolecule m1 m2) i1 (i2 + m1Length)  b


{------------------------------------------------------------------------------}
incrementAtomMap :: (Map Int Atom) -> Int -> (Map Int Atom)
incrementAtomMap map i = Map.mapKeys (+i) $ Map.map (\a -> incrementAtom a i) map




-- hasClosure
-- True if molecule contains any atom with a colosure marker
{------------------------------------------------------------------------------}
hasClosure :: Molecule -> Bool
hasClosure m = List.elem True $ List.map (isClosure) markers
    where markers = List.foldr ((++) . Set.toList . atomMarkerSet) [] $
                    List.map snd $ Map.toList (atomMap m)
          isClosure mk  = case mk of Closure {} -> True ; _ -> False


-- addMolecule
-- Combines two molecules and connects bond-markers if required
-- Otherwise, adds as disconnected structure.  This is not really meant to be accessed
-- directly.
{------------------------------------------------------------------------------}
addMolecule :: Molecule -> Molecule -> Molecule
addMolecule m1 m2 = if (moleculeHasError m1) then m1 else
                    if (moleculeHasError m1) then m2 else m12
    where m12 = Small {atomMap=newMap, molMarkerSet=newMarkerSet, molPropertyMap=Map.empty}
          newMap = Map.union (atomMap m1) incMap
          newMarkerSet = Set.union (molMarkerSet m1) (molMarkerSet m2)
          m1Length = Map.size $ atomMap m1
          incMap = incrementAtomMap (atomMap m2) m1Length




--makeMoleculeFromAtom
{------------------------------------------------------------------------------}
makeMoleculeFromAtom:: Atom -> Molecule
makeMoleculeFromAtom a = Small {atomMap=(Map.singleton 0 (markAtom a $ Label 0))
                              , molMarkerSet=Set.empty
                              , molPropertyMap=Map.empty}


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
fillMoleculeValence m = Map.foldWithKey foldMol m $ atomMap m
  where foldMol k atomMap' molecule = fillValenceAtIndex molecule k


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

 -- giveMoleculeWaarning
 -- Adds an error marker containing a string to the molecule.
{------------------------------------------------------------------------------}
giveMoleculeWarning :: Molecule -> String -> Molecule
giveMoleculeWarning m warning = markMolecule m $ Warning warning

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






{------------------------------------------------------------------------------}
{-------------------------------Typeclass Intances-----------------------------}
{------------------------------------------------------------------------------}

instance Show Molecule where
    show m = if (moleculeHasError m) then ("Molecule has error.") else case m of
       Small {atomMap=atoms, molMarkerSet=mm, molPropertyMap=mp} -> "\nIs a small molecule with formula: "
                ++ (\(Right a) -> a) (molecularFormula $ m) ++ "\n"
                ++ (List.foldr (\b ->  (++) ((show $ fst b) ++ " -- "
                    ++ (show $ snd b))) "" (Map.toList atoms))
                ++ (List.foldr (\b ->  (++) (show b)) "" (Set.toList mm)) ++ "\n"
                ++ (Map.fold (\b ->  (++) (show b)) "" mp) ++ "\n"
       Markush  {}   -> "Is a markush"
       Polymer  {}   -> "Is a polymer"
       Biologic {}   -> "Is a Biologic"


{------------------------------------------------------------------------------}
instance Eq Molecule where
    a == b = True --Obviously that needs fixing!

