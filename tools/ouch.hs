--
-- So far simple SMILES file CLI printer without options.
--
-- Usage: ouch <SMILES_file>
--

import Control.Monad(when)
import Data.Functor((<$>))
import Ouch.Input.Smiles(readSmi)
import Ouch.Property.Composition(molecularFormula, molecularWeight)
import Ouch.Structure.Molecule
import System.Environment
import System.Exit(exitFailure)
import Text.Printf(printf)
import qualified Ouch.Property.Builder as OPB

showMolecule :: String -- ^ SMILES string
             -> String -- ^ description
showMolecule smile = smile ++ "\n  --> " ++ showProperty molecularFormula mol ++ ", " ++ showMolWeight mol
  where
    mol = readSmi smile

showProperty :: OPB.Property -> Molecule -> String
showProperty p m = show $ getPropertyValue p m

showMolWeight :: Molecule -> String
showMolWeight mol = let (OPB.DoubleValue d) = getPropertyValue molecularWeight mol
                    in printf "%.2f" d

getPropertyValue :: OPB.Property -> Molecule -> OPB.Value
getPropertyValue prop m = case OPB.value prop of
                            Left  v -> v
                            Right f -> f m

main :: IO ()
main = do
  args <- getArgs
  when (null args) $ do
    putStrLn "ouch: SMILES file operand required"
    exitFailure
  smiles <- lines <$> readFile (head args)
  mapM_ (putStrLn . showMolecule) smiles

