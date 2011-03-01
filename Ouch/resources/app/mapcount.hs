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
      {-p = growPath m (U.singleton 0)-}
  putStrLn $ show m
  {-putStrLn version-}
  {-putStrLn $ Prelude.head args-}
  {-putStrLn $ "Number of paths from index 0: " ++ (show $ V.length p) ++ "\n"-}
  {-sequence_ $ V.foldr (\a acc -> (do-}
      {-putStrLn $ show $ U.length a):acc) [] p-}
  {-putStrLn $ show $ pathLength $ longestLeastPath m-}
  {-putStrLn $ show $ V.length $ longestPaths m-}
  {-sequence_ $ M.foldWithKey (\k a acc -> (do-}
        {-putStr $ (show k) ++ ": " -}
        {-putStrLn $ show $ V.length $ growPath m (U.singleton k)):acc) [] $ atomMap m-}
  {-putStrLn $ "Map Length: " ++ (show $ p)-}
  {-putStrLn $ debugShow m-}


