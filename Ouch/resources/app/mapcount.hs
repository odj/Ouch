import Ouch
import System.IO
import Control.Monad
import System.Environment
import Data.Map as M
import Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Maybe



main = do 
  arg:_ <- getArgs
  input <- readFile arg
  let readMol s = read s :: Molecule
  sequence (L.map (putStrLn . show . readMol) (lines input))
