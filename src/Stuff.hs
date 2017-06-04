{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE ExplicitForAll #-}
--{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
--{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RankNTypes #-}
--{-# LANGUAGE RebindableSyntax #-}
--{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE OverloadedLists #-}
--{-# LANGUAGE NamedFieldPuns #-}
module Stuff ( (<>), Semiring (..), Ring (..)
  , (<<<), uidx, uJust,debug
  , Number, Array
  , Tuple (..)
  , List (..)
  , to0fix, to1fix, to2fix, to3fix, to5fix
  , toNumber, intFromString, numberFromString
  , sqr, div', mod', divMod'
  )  where

import Prelude
import qualified Data.Vector.Unboxed as A
  ( Vector, length, fromList, toList, unsafeIndex, create, replicate
  , singleton, map, foldl, zipWith )
import Data.Maybe (  Maybe (..), fromJust )
import Text.Printf ( printf )
import Text.Read ( readMaybe )

--------------------------------------------------------------
-- adapting for PureScript
import Data.Semigroup ( (<>) )
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

type Number = Double
type Array a = A.Vector a
uidx :: Array Number -> Int -> Number
uidx = A.unsafeIndex
uJust :: forall a. Maybe a -> a
uJust = fromJust
(<<<) :: (b -> c) -> (a -> b) -> a -> c
(<<<) = (.)
infixr 9 <<<

-- | A simple product type for wrapping a pair of component values.
-- type Tuple a b = (a, b)
data Tuple a b = Tuple a b

-- List, PureScript does not provide sugar
type List a = [a]

-- | Allows `Tuple`s to be rendered as a string with `show` whenever there are
-- | `Show` instances for both component types.
instance (Show a, Show b) => Show (Tuple a b) where
  show (Tuple a b) = "(Tuple " <> show a <> " " <> show b <> ")"

class Semiring a where
  add  :: a -> a -> a
  zero :: a
  mul  :: a -> a -> a
  one  :: a
class Semiring a => Ring a where
  sub :: a -> a -> a
--------------------------------------------------------------

toNumber :: Int -> Number
toNumber = fromIntegral

to0fix :: Number -> String
to0fix = printf "%4.0f"
to1fix :: Number -> String
to1fix = printf "%6.1f"
to2fix :: Number -> String
to2fix = printf "%7.2f"
to3fix :: Number -> String
to3fix = printf "%8.3f"
to5fix :: Number -> String
to5fix = printf "%10.5f"

sqr :: Number -> Number
sqr x = x*x

-- | generalisation of 'div' to any instance of Real
div' :: Number -> Number -> Int
div' n d = floor ( n /  d)

-- | generalisation of 'divMod' to any instance of Real
divMod' :: Number -> Number -> (Tuple Int Number)
divMod' n d = (Tuple f (n - (toNumber f) * d)) where
    f = div' n d

-- | generalisation of 'mod' to any instance of Real
mod' :: Number -> Number -> Number
mod' n d = n - (toNumber f) * d where
    f = div' n d

-- | convert from String to Number and Int
numberFromString :: String -> Maybe Number
numberFromString s = readDouble s where
  readDouble = readMaybe :: String -> Maybe Double
intFromString :: String -> Maybe Int
intFromString s = readInt s where
  readInt = readMaybe :: String -> Maybe Int

