import Ouch
import Ouch.Output.SDF
import System.IO
import System.Environment
import Data.Map as M
import Data.List as L

-- A command-line utility for creating a series of spiro-cyclopropane
-- dendrimers.  Takes a single inteeger [n] from the command line
-- and writes [n] structrues to std-out


cyclopropane = read "C1CC1" :: Molecule
ethenyl = read "CC" :: Molecule

spiroMethod = Just $ FilterMethod { firstApply = Nothing
                                   , lastApply  = Nothing
                                   , molFilter  = spirolize
                                   }
spirolize (m:_) = let
  openIndices = M.keys $ M.filter (openValenceSelector m) (atomMap m)
  addEthyl m i = connectMoleculesAtIndicesWithBond m i ethenyl 0 Single
  cyclize i m = addBond m i i1 Single 
    where i1 = M.size (atomMap m) - 1
  in (L.foldr (\i acc -> cyclize i $ addEthyl acc i) m openIndices):[]



-- Prints a list as SMILES
ppList [] = return ()
ppList (x:xs) = do
  putStrLn $ show x
  ppList xs

main = do
  arg:_ <- getArgs
  let n = read arg :: Int
  let mols = L.foldl (\acc i -> acc >#> spiroMethod) [cyclopropane] [1..n]
  {-putStrLn $ Ouch.Output.SDF.sdf mols-}
  ppList mols
