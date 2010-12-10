import Ouch.Structure.Molecule
import Ouch.Input.Smiles
import Ouch.Enumerate.Formula
import System.IO
import Data.List
import System.Environment
import Control.Parallel.Strategies
import Control.Parallel



main = do 
  arg:_ <- getArgs
  let f = read arg::Formula
      l = expand f
  l `seq` ppList l
  putStrLn $ show $ length l

ppList [] = return ()
ppList (x:xs) = do
  putStrLn $ show x
  ppList xs
