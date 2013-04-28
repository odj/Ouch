{-# LANGUAGE RecordWildCards, DeriveDataTypeable #-}
--
-- CLI OUCH tool.
--
-- Common flags:
--   -f --file=FILE          List compounds from SMILES file
--   -e --enumerate=FORMULA  Enumerate molecular formula
--   -? --help               Display help message
--   -V --version            Print version information
--

import Control.Monad(unless)
import Data.Functor((<$>))
import Ouch.Enumerate.Formula(expand, Formula)
import Ouch.Input.Smiles(readSmi)
import Ouch.Property.Composition(molecularFormula, molecularWeight)
import Ouch.Structure.Molecule
import System.Console.CmdArgs.Implicit
import System.Environment(getArgs, withArgs)
import Text.Printf(printf)
import qualified Ouch.Property.Builder as OPB

ouchVersion, ouchProgram, ouchSummary :: String
ouchVersion   = "0.0.1"
ouchProgram   = "ouch"
ouchSummary   = ouchProgram ++ " v" ++ ouchVersion

data CmdOptions = CmdOptions {
  file      :: String,
  enumerate :: String
} deriving (Show, Data, Typeable)

defOpts :: CmdOptions
defOpts = CmdOptions
  { file      = "" &= typ "FILE" &= help "List compounds from SMILES file"
  , enumerate = "" &= typ "FORMULA" &= help "Enumerate molecular formula"
  } &= summary ouchSummary
    &= program ouchProgram
    &= help "OUCH CLI tool"


showMolecule :: Molecule -> String
showMolecule mol = showProperty molecularFormula mol ++ ", " ++ showMolWeight mol

showSmile :: String -- ^ SMILES string
          -> String -- ^ description
showSmile smile = smile ++ "\n  --> " ++ showMolecule (readSmi smile)

showProperty :: OPB.Property -> Molecule -> String
showProperty p m = show $ getPropertyValue p m

showMolWeight :: Molecule -> String
showMolWeight mol = let (OPB.DoubleValue d) = getPropertyValue molecularWeight mol
                    in printf "%.2f" d

getPropertyValue :: OPB.Property -> Molecule -> OPB.Value
getPropertyValue prop m = case OPB.value prop of
                            Left  v -> v
                            Right f -> f m

-- | If no arguments are passed to a program, displays help message and exits,
-- returns parsed arguments otherwise.
cmdArgsOrHelp :: Data a => a -> IO a
cmdArgsOrHelp opts = do
  mainArgs <- getArgs
  (if null mainArgs then withArgs ["-?"] else id) (cmdArgs opts)

main :: IO ()
main = do
  CmdOptions{..} <- cmdArgsOrHelp defOpts
  unless (null file) $ do
    smiles <- lines <$> readFile file
    mapM_ (putStrLn . showSmile) smiles
  unless (null enumerate) $
    mapM_ print (expand (read enumerate::Formula))


