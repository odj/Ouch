import Ouch.Structure.Molecule
import System.IO
import System.Environment
import Ouch.Enumerate.Method

main = do
  arg:_ <- getArgs
  contents <- readFile arg
  let molecules = map (\a -> read a :: Molecule) $ lines contents
      uniqueMols = molecules >#> makeUnique
  ppList uniqueMols

-- Prints a list as SMILES     
ppList [] = return ()
ppList (x:xs) = do
  putStrLn $ show x
  ppList xs  
