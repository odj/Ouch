import Ouch.Test.Methods
import Ouch.Structure.Molecule
import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Input.Smiles
import Ouch.Data.Atom

import Data.Either
import Data.Maybe


-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
main = do
    putStrLn $ performTests tests
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
