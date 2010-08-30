{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
    String -

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

--------------------------------------------------------------------------------
-------------------------------------------------------------------------------}

module Ouch.Text.String (
    Justify (..)
  , padString
  ) where



import Data.Maybe
import Data.Char
import Data.List


{------------------------------------------------------------------------------}
{-------------------------------Date Types-------------------------------------}
{------------------------------------------------------------------------------}

data Justify = LeftJustify | RightJustify | CenterJustify


{------------------------------------------------------------------------------}
{-------------------------------Functions--------------------------------------}
{------------------------------------------------------------------------------}

padString :: Justify -> Int -> Char -> String -> String
padString jst width c s = output
  where str = take width s
        padWidth = (width - (length str))
        padStr = replicate padWidth c
        lPadStr | odd padWidth = replicate ((padWidth - 1) `div` 2) c
                | otherwise = replicate (padWidth `div` 2) c
        rPadStr | odd padWidth = replicate (1 + (padWidth - 1) `div` 2) c
                | otherwise = replicate (padWidth `div` 2) c
        output = case jst of
          LeftJustify   -> str ++ padStr
          RightJustify  -> padStr ++ str
          CenterJustify -> lPadStr ++ str ++ rPadStr
