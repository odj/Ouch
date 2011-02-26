{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Structure.Molecule
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

{-# LANGUAGE BangPatterns #-}

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
     , bondBetweenIndices
     , debugShow
     , setAtom
     , getAtomAtIndex
     , getBondMap
     , getName
     , getPropertyForKey
     , getMoleculeError
     , addMarkerToAtomAtIndex
     , addMolecule
     , makeMoleculeFromAtom
     , fillMoleculeValence
     , fillValenceAtIndex
     , freeValenceAtIndex
     , hasHangingClosure
     , cyclizeMolecule
     , giveMoleculeError
     , moleculeHasError
     , numberBondsAtIndex
     , occupiedValenceAtIndex
     , numberBondsToAtomsAtIndex
     , numberBondsToHeavyAtomsAtIndex
     , numberBondsToRadicalsAtIndex
     , numberBondsToHydrogensAtIndex
     , connectMoleculesAtIndicesWithBond
     , replaceAtomAtIndexWithBond
     , incrementAtomMap
     , removeAtoms
     , updateBondSet
     , emptyMolecule
     , (>@>)
     ) where


import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker


import {-# SOURCE #-} Ouch.Property.Builder
import {-# SOURCE #-} Ouch.Property.Extrinsic.Fingerprint
import {-# SOURCE #-} Ouch.Input.Smiles
import Ouch.Output.Smiles
import Ouch.Enumerate.Method

import Data.Either
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

-- | Infix function to filter out molecules with errors when performing
-- some operation on a molecule of dubious origin.
(>>>) :: Molecule -> Molecule -> Molecule
(>>>) !mIn !mOut = if (moleculeHasError mIn) then mIn else mOut

emptyMolecule = Molecule Map.empty Set.empty Map.empty


{------------------------------------------------------------------------------}
-- | Returns the atom at the specified index, if it exists.
getAtomAtIndex :: Molecule -> Int -> Maybe Atom
getAtomAtIndex m i = Map.lookup i $ atomMap m

{------------------------------------------------------------------------------}
-- | Inserts atom into given molecule at whatever index the atom thinks it should
-- be at.  If the index is invalid or non-existant, then it adds to the next
-- atom index.  If the molecule is flagged with an error, then returns the (corrupt)
-- molecule unchanged.
setAtom :: Atom -> Molecule -> Molecule
setAtom !a !m  = m >>> mOut
  where mOut | (isJust atomIndex) /= True  =  addAtom a m
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
addAtom !a !m = m >>> mOut
    where mOut =  m {atomMap = Map.insert atomNumber newAtom $ atomMap m}
          atomNumber = Map.size $ atomMap m
          !newAtom = a {atomMarkerSet=(Set.insert (Label atomNumber) $ atomMarkerSet a) }


addMolMarker :: Molecule -> MoleculeMarker -> Molecule
addMolMarker m mm = m {molMarkerSet = Set.insert mm $ molMarkerSet m}

addProperty :: Molecule -> Property -> Molecule
addProperty m p = m {molPropertyMap = Map.insert (propertyKey p) p $ molPropertyMap m}

addPropertyFromFunction :: Molecule -> (Molecule -> Maybe Property) -> Molecule
addPropertyFromFunction m f = let
  prop = f m
  output = case prop of
    Just p ->  addProperty m p
    Nothing -> m
  in output

getPropertyForKey :: Molecule -> String -> Maybe Property
getPropertyForKey m k = Map.lookup k $ molPropertyMap m


occupiedValenceAtIndex :: Molecule -> Int -> Integer
occupiedValenceAtIndex m i = case atomAtIndex m i of
  Just a  -> occupiedValence a
  Nothing -> 0

freeValenceAtIndex :: Molecule -> Int -> Integer
freeValenceAtIndex m i = case atomAtIndex m i of
  Just a  -> (fst $ valence a) - (occupiedValence a)
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


-- | Find the bond type between two atom indices
bondBetweenIndices :: Molecule -> Int -> Int -> NewBond
bondBetweenIndices m i1 i2 = let
  a1 = getAtomAtIndex m i1
  a2 = getAtomAtIndex m i2
  newBond = case a1 of
              Nothing -> NoBond
              Just atom1 -> case a2 of
                Nothing -> NoBond
                Just atom2 -> bondBetweenAtoms atom1 atom2
  bondBetweenAtoms a b | (Set.size filteredSet1) /= 1 = NoBond
                       | (Set.size filteredSet2) /= 1 = NoBond
                       | newBondForBond (Set.findMax filteredSet1) ==
                         newBondForBond (Set.findMax filteredSet2) = newBondForBond (Set.findMax filteredSet2)
                       | otherwise = NoBond
  filteredSet1 = Set.filter (hasBondTo i2) (atomBondSet $ fromJust a1)
  filteredSet2 = Set.filter (hasBondTo i1) (atomBondSet $ fromJust a2)
  hasBondTo i b = (bondsTo b) == i
  in newBond



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
addBond !m !i1 !i2 !b  =  m >>> mOut
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

getBondMap :: Molecule -> Map (Int, Int) NewBond
getBondMap m = cleanMap
  where atoms = atomMap m
        cleanMap = Map.mapKeys (\a -> sortTuple a) atomsToMap
        atomsToMap = Map.fold (\a m -> Map.union m $ bondSetToMap a) Map.empty atoms
        bondSetToMap a = Set.fold (\bond m -> Map.union m $ bondToMapElem bond $ ind a)
                         Map.empty $ atomBondSet a
        bondToMapElem b i = Map.singleton (i, bondsTo b) $ showNewBond b
        ind a = fromJust $ getIndexForAtom a
        sortTuple a@(a1, a2) | a1 > a2 = (a2, a1) | a1 < a2 = a | otherwise = a
        showNewBond b = case b of
          Sigma {}    -> Single
          Pi {}       -> Double
          PiPi {}     -> Triple
          Aromatic {} -> AromaticOnly


-- addLonePairAtIndex
-- Create a new lone-pair centered on the atom.
{------------------------------------------------------------------------------}
addLonePairAtIndex :: Molecule -> Int -> Molecule
addLonePairAtIndex m i  = m >>> mOut
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = LonePair Set.empty Set.empty

--addHydrogenAtIndex
{------------------------------------------------------------------------------}
addHydrogenAtIndex :: Molecule -> Int -> Molecule
addHydrogenAtIndex m i  = m >>> mOut
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = (Element 1 0) Set.empty Set.empty

--addElectronAtIndex
{------------------------------------------------------------------------------}
addElectronAtIndex :: Molecule -> Int -> Molecule
addElectronAtIndex m i  = m >>> mOut
    where mOut = addBond (addAtom newAtom m) i (Map.size $ atomMap m) Single
          newAtom = Electron Set.empty Set.empty


-- addUnfilled
-- Create a new unfilled orbital centered on the atom.
-- Return atom and list containing new unfilled orbital.
{------------------------------------------------------------------------------}
addUnfilled :: Molecule -> Int -> Molecule
addUnfilled m i  = m >>> mOut
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
fillValenceAtIndex !m !i = m >>> mOut
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
cyclizeMolecule m = m >>> mOut
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
hasHangingClosure m = if (moleculeHasError m) then False else output
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
cyclizeMoleculeAtIndexesWithBond !m !i1 !i2 !b =
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
connectMoleculesAtIndicesWithBond !m1 !i1 !m2 !i2 !b = m1 >>> mOut
  where mOut = connectMolecules m1 i1 m2 i2 b
        connectMolecules m1 i1 m2 i2 b | m2isEmpty = m1
                                       | errorTest = giveMoleculeError m1 "Could not connect molecules, invalid index"
                                       | otherwise = output
        a1 = getAtomAtIndex m1 i1
        a2 = getAtomAtIndex m2 i2
        errorTest = (isJust a1 && isJust a1) /= True
        m1Length = Map.size $ atomMap m1
        m2isEmpty = (Map.size $ atomMap m2) == 0
        output = addBond (addMolecule m1 m2) i1 (i2 + m1Length)  b


{------------------------------------------------------------------------------}
incrementAtomMap :: (Map Int Atom) -> Int -> (Map Int Atom)
incrementAtomMap !map !i = Map.mapKeys (+i) $ Map.map (\a -> incrementAtom a i) map



removeAtoms :: Molecule -> (Atom -> Bool) -> Molecule
removeAtoms !m !f = incrementAtoms
  where filterMol  = m {atomMap = Map.filter ((/=True) . f) $ atomMap m}
        filterBond = filterMol `seq` Map.fold (\a acc -> setAtom (updateBondSet acc a) acc) filterMol $ atomMap filterMol
        incrementAtoms = filterBond `seq` Map.foldlWithKey f filterBond $ atomMap filterBond
            where f !acc !i1 !a = let
                    i2 = (Map.findIndex (fromJust $ getIndexForAtom a) $ atomMap acc)
                    in mapIndex acc i1 i2
        mapIndex !m !i1 !i2 = newAtom `seq` mapBondxs (setAtom newAtom m') i1 i2
            where newAtom = incrementAtomMarker (fromJust $ getAtomAtIndex m i1) (i2 - i1)
                  m' = m {atomMap = Map.delete i1 $ atomMap m}
                  mapBondxs !m !i1 !i2 = m {atomMap = Map.map (\a -> mapBonds a i1 i2) $ atomMap m}
                  mapBonds !a !i1 !i2 = a {atomBondSet = Set.map (\b -> mapBond b i1 i2) $ atomBondSet a}
                  mapBond !b !i1 !i2 | bondsTo b == i1 = b {bondsTo = i2} | otherwise = b
                  incrementAtomMarker !a !i = markAtom a (Label $ (fromJust $ getIndexForAtom a) + i)

updateBondSet :: Molecule -> Atom -> Atom
updateBondSet !m !a = a {atomBondSet = valid}
    where keyIsValid !m !b = case Map.lookup (bondsTo b) $ atomMap m of
                           Just _  -> True
                           Nothing -> False
          (valid, _) = Set.partition (keyIsValid m) $ atomBondSet a


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
                    if (moleculeHasError m2) then m1 else m12
    where m12 = Molecule {atomMap=newMap, molMarkerSet=newMarkerSet, molPropertyMap=Map.empty}
          m1err = giveMoleculeError m1 "Tried to added a molecule with error marker."
          newMap = Map.union (atomMap m1) incMap
          newMarkerSet = Set.union (molMarkerSet m1) (molMarkerSet m2)
          m1Length = Map.size $ atomMap m1
          incMap = incrementAtomMap (atomMap m2) m1Length


-- replaceAtomAtIndexWithBond
-- Combines the two molecules then substitutes the atom at i1 with the
-- atom at i2 (from the old molecule's numbering system).  Most commonly, m2
-- will be a single atom
{------------------------------------------------------------------------------}
replaceAtomAtIndexWithBond :: Molecule -> Int -> Molecule -> Int -> Molecule
replaceAtomAtIndexWithBond m1 i1 m2 i2 = m1

--makeMoleculeFromAtom
{------------------------------------------------------------------------------}
makeMoleculeFromAtom:: Atom -> Molecule
makeMoleculeFromAtom a = Molecule {atomMap=(Map.singleton 0 (markAtom a $ Label 0))
                              , molMarkerSet=Set.empty
                              , molPropertyMap=Map.empty}



-- fillMoleculeValence
-- Fills valence with hydrogens and lone-pairs to give neutral species
-- If valence is already complete, returns molecule unchanged
-- If cannot fill because of an error, adds error to molecule and gives it back.
{------------------------------------------------------------------------------}
fillMoleculeValence :: Molecule -> Molecule
fillMoleculeValence !m = Map.foldrWithKey foldMol m $ atomMap m
  where foldMol !k !atomMap' !molecule = fillValenceAtIndex molecule k


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
getMoleculeError m  = let err =  Set.filter (\a -> EQ == compare a (MError "")) $  molMarkerSet m
                      in  if Set.size err == 0 then Nothing
                          else Just $ molMarker $ Set.findMax err

--markMolecule
{------------------------------------------------------------------------------}
markMolecule :: Molecule -> MoleculeMarker -> Molecule
markMolecule m mm = newMol
    where newSet = Set.insert mm (molMarkerSet m)
          newMol = m {molMarkerSet=newSet}

getName :: Molecule -> Maybe String
getName m | Set.size nameSet == 0 = Nothing
          | otherwise = Just $ molMarker $ Set.findMax nameSet
  where nameSet = Set.filter (\a ->  EQ == compare a (Name "")) $ molMarkerSet m


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
addMarkerToAtomAtIndex !m !i !am = if (moleculeHasError m) then m else case atom of
    Nothing -> markMolecule m warning
    Just a  -> m {atomMap = Map.insert i (markAtom a am) (atomMap m)}
    where atom = atomAtIndex m i
          warning = Warning $ "Unable to add marker to atom at position " ++ show i

(>@>) :: Molecule -> Maybe AtomMarker -> Molecule
(>@>) a m = case m of Nothing -> a
                      Just m' -> addMarkerToAtomAtIndex a 0 m'




{------------------------------------------------------------------------------}
{-------------------------------Typeclass Intances-----------------------------}
{------------------------------------------------------------------------------}

instance Show Molecule where
  show m = writeSmiles m


instance Read Molecule where
  readsPrec _ s = [(readSmi s, "")]


debugShow m = if (moleculeHasError m) then ("Molecule has error.") else case m of
   Molecule {atomMap=atoms, molMarkerSet=mm, molPropertyMap=mp} ->
           --"\nIs a small molecule with formula: "
           -- ++ (\(Right a) -> a) (molecularFormula $ m) ++ "\n" ++
               (List.foldr (\b ->  (++) ((show $ fst b) ++ " -- "
                ++ (show $ snd b))) "" (Map.toList atoms))
            ++ (List.foldr (\b ->  (++) (show b)) "" (Set.toList mm)) ++ "\n"
            ++ (Map.fold (\b ->  (++) (show b)) "" mp) ++ "\n"




{------------------------------------------------------------------------------}
instance Eq Molecule where
    (==) a b = if (Map.size $ atomMap a) /= (Map.size $ atomMap b) then False
               else if (show a) /= (show b) then False
               else True
instance Ord Molecule where
   compare a b = compare (show a) (show b)




