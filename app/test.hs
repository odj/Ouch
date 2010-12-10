import System.Environment
import Data.Maybe
import Data.Either
import Ouch.Structure.Molecule
import Ouch.Structure.Marker
import Ouch.Structure.Atom
-- import Ouch.Input.Smiles
import Ouch.Enumerate.Method
import Data.Map as Map
import Data.List as List
import Data.Set as Set
import Data.ByteString.Lazy as L
import Ouch.Property.Extrinsic.Fingerprint
import Ouch.Property.Composition
import Ouch.Input.Smiles
import Data.Binary.Builder as B
import System.IO

carbon = makeMoleculeFromAtom $ Element 6 0 Set.empty Set.empty

mth = Just $ AddMethod
  { firstApply=Nothing
  , lastApply=makeUnique
  , selector=openValenceSelector
  , addList=[(Single, carbon)]
  }

--makeUnique = Just $ FilterMethod
--  { firstApply=Nothing
--  , lastApply=Nothing
--  , molFilter=fingerprintFilterBuilder (\m -> writeCanonicalPath m)
--  }

alk :: Int -> [Molecule]
alk i = List.foldr (\enum mols -> enum mols) [carbon]
                 $ List.replicate (i-1) (>#>  mth)

main = do
  arg:_ <- getArgs
  let ns = [1..read arg::Int]
      lengths = List.map (List.length . alk) ns
  Prelude.putStrLn $ show $ List.zip ns lengths


mth' = Just $ AddMethod
  { firstApply=Nothing
  , lastApply=Nothing
  , selector=(openValenceSelector >&&> elementSelector "O")
  , addList=[(Single, carbon)]
  }
