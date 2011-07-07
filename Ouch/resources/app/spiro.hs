#!/usr/bin/env runghc
-- build with 'ghc --make spiro.hs -O2 -o spiro'
import Ouch
import Ouch.Output.SDF
import System.IO
import System.Exit
import System.Environment
import Data.Map as M
import Data.List as L
import Data.Set as S

-- A command-line utility for creating a series of spiro-cycloalkane
-- dendrimers.  Takes two integers from the command line: a ring size to
-- template with and a generation number.  Returns a single structure as 
-- a SMILES string.


cyclopropane = read "C1CC1" :: Molecule
ethenyl = read "CC" :: Molecule

carbon = Element 6 6 S.empty S.empty

spirolize branch (m:_) = let
    openIndices = M.keys $ M.filter (openValenceSelector m) (atomMap m)
    addBranch m i = connectMoleculesAtIndicesWithBond m i branch 0 Single
    cyclize i m = addBond m i i1 Single
      where i1 = M.size (atomMap m) - 1
    in (L.foldr (\i acc -> cyclize i $ addBranch acc i) m openIndices):[]

exitWithUsage = do
    putStrLn "Usage: spiro <m> <n>"
    exitSuccess

growWith atom molecule 
    | size == 0  = addAtom atom molecule
    | otherwise  = addBond (addAtom atom molecule) (size -1) size Single
    where size = M.size $ atomMap molecule


makeCycle molecule =
    addBond molecule 0 (size - 1) Single
    where size = M.size $ atomMap molecule


makeChain atom size molecule
    | currentSize /= size = makeChain atom size $ growWith atom molecule
    | otherwise           = molecule
    where currentSize = M.size $ atomMap molecule

makeRing atom size molecule = makeCycle $ makeChain atom size molecule

-- Prints a list as SMILES
ppList [] = return ()
ppList (x:xs) = do
    putStrLn $ show x
    ppList xs

main = do
    args <- getArgs
    if (length args /= 2)
        then exitWithUsage
        else run args

run args = do
    let m = read (args!!0) :: Int
        n = read (args!!1) :: Int
        core = makeRing carbon m emptyMolecule
        branch = makeChain carbon (m - 1) emptyMolecule
        spiroMethod = Just $ FilterMethod { firstApply = Nothing
                                           , lastApply  = Nothing
                                           , molFilter  = spirolize branch
                                           }
        mols = L.foldl (\acc i -> acc >#> spiroMethod) [core] [1..n]
    ppList mols


