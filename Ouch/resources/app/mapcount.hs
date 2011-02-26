import Ouch
import System.IO
import Control.Monad
import System.Environment
import Data.Map as M
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Maybe

  
main = do 
  args <- getArgs
  let m = read (Prelude.head args)::Molecule
      {-ll = longestLeastPath m-}
      p = growPath m (U.singleton 0)
  putStrLn version
  putStrLn $ "Map Length: " ++ (show $ p)
  {-putStrLn $ debugShow m-}


