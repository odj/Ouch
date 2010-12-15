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
  let n = read arg::Int
      l = inputs n
      results = l `pseq` parMap rseq (\input -> (ppList (expand (read input::Formula)))) l
      zipped = results `pseq` zip results l
  zipped `pseq` mapM (\z -> writeFile (((snd z) ++ ".out")::FilePath) (fst z)) zipped


ppList [] = ""
ppList (x:xs) = (show x) ++ "\n"
                         ++ ppList xs


inputs n = let
  range = [1..(n-1)]
  elements = ["Si", "N", "O", "Cl"]
  combinedFormulas = [ "C" ++ (show (n-num))
                           ++ elem
                           ++ (show num) | num <- range, elem <- elements]
  singleElem = ["C" ++ (show n)] ++
               [elem ++ (show n) | elem <- elements]
  in singleElem ++ combinedFormulas             
