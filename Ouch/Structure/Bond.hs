{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Bond - a module to define bond data type
    
    Copyright (c) 2010 Orion D. Jankowski
    
    This file is part of Ouch, a chemical informatics toolkit
    written entirely in Haskell.

    Ouch is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Ouch is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Ouch.  If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------
------------------------------------------------------------------------------}
{-# LANGUAGE ForeignFunctionInterface, CPP, Generics #-}

module Ouch.Structure.Bond (
      Bond(..)
    ) where
    
-- This line terminates recursive import sequences        
import {-# SOURCE #-} Ouch.Structure.Atom
import Ouch.Structure.Marker
import Data.Set as Set




{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}

data Bond = Sigma {bondsTo::Atom}
          | Pi {bondsTo::Atom} 
          | Aromatic {bondsTo::Atom}
          | Delta {bondsTo::Atom}
          | Hbond {bondsTo::Atom}
          | Ionic {bondsTo::Atom}
          | Antibond {bondsTo::Atom}
          | Any {bondsTo::Atom}



      
{------------------------------------------------------------------------------}
{-------------------------------Typeclass Intances-----------------------------}
{------------------------------------------------------------------------------}    
instance Show Bond where
  show b = case b of
      Sigma {} -> " -" ++ desc
      Pi {} -> " =" ++ desc
      Aromatic {} -> " ~" ++ desc
      Delta {} -> " delta-" ++ desc
      Hbond {} -> " H." ++ desc
      Ionic {} -> " +/-" ++ desc
      Antibond {} -> " !" ++ desc
      where atom = bondsTo b
            desc = case atom of
                Element {} -> atomicSymbolForAtom atom
                LonePair {} -> "LP"
                Electron {} -> "EL"
                Unfilled {} -> "UF"
  
 
instance Eq Bond where
  a == b = True
          
