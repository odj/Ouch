import Ouch.Enumerate.Formula
import System.IO
import System.Environment

-- A command-line utility for enumerating
-- a molecular formula.
main = do 
  args <- getArgs
  let f = read (head args)::Formula
      l = expand f
  if args == [] then putStrLn "Enter a molecular formula."
     else l `seq` ppList l


-- Prints a list as SMILES     
ppList [] = return ()
ppList (x:xs) = do
  putStrLn $ show x
  ppList xs
