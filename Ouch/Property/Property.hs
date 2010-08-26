{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Property - a module to manage property data types

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

module Ouch.Property.Property (
   Property(..)
 , Value(..)
   ) where


import {-# SOURCE #-} Ouch.Structure.Molecule
import Data.Maybe
import Data.Set as Set
import Data.List as List
import Data.Map as Map


{------------------------------------------------------------------------------}

data Property = Property {propertyKey::String
                        , value::Value
                        , func::Maybe (Molecule -> Maybe Property)}


instance Ord Property where
    compare a b | (propertyKey a) >  (propertyKey b) = GT
                | (propertyKey a) <  (propertyKey b) = LT
                | (propertyKey a) == (propertyKey b) = EQ

instance Eq Property where
    a == b = (propertyKey a) == (propertyKey b)


instance Show Property where
  show a = (show $ propertyKey a) ++ ": " ++ (show $ value a)


data Value =
    IntegerValue    Integer
  | DoubleValue     Double
  | StringValue     String
  | BoolValue       Bool
  | TupleArrayValue [(Double, Double)]
  deriving (Eq, Ord, Show)


