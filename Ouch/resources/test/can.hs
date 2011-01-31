import System.IO
import System.Environment
import Ouch.Structure.Molecule

-- A command-line utility for enumerating
-- a molecular formula.
main = do 
  arg:_ <- getArgs
  contents <- readFile arg
  let l = lines contents 
  ppSmi l 

-- Prints a list as SMILES     
ppSmi [] = return ()
ppSmi (x:xs) = do
  putStrLn $ show $ (read x :: Molecule)
  ppSmi xs
