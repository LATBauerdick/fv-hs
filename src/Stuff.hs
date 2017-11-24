-- {-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE ExplicitForAll #-}
--{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
-- {-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
--{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RankNTypes #-}
--{-# LANGUAGE RebindableSyntax #-}
--{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE OverloadedLists #-}
--{-# LANGUAGE NamedFieldPuns #-}
module Stuff ( (<>), Semiring (..), Ring (..)
  , (<<<), uidx, uJust, debug, trace, unsafePartial
  , Number, Array
  , Tuple (..)
  , List, range, fromList
  , to0fix, to1fix, to2fix, to3fix, to5fix
  , iflt, irem
  , toNumber, intFromString, numberFromString
  , sqr, div', mod', divMod'
  , normals
  )  where

import Prelude
import qualified Data.Vector.Unboxed as A
  ( Vector, Unbox, length, fromList, toList, unsafeIndex
    , create, replicate
    , singleton, map, foldl, zipWith
    , replicateM, concat, take )
import Data.Maybe (  Maybe (..), fromJust )
import Text.Printf ( printf )
import Text.Read ( readMaybe )
import System.Random ( RandomGen, random )
import qualified Data.List as L ( take )

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
unsafePartial :: forall a. a -> a
unsafePartial x = x

-- filter list of objects given list of indices in [a]
-- return list with only those b that have  indices that  are in rng [a]
iflt :: ( Eq a, Enum a, Num a ) => [a] -> [b] -> [b]
iflt rng hl =
  [h | (h, i) <- zip hl [0..], i `elem` rng ]

irem :: (Eq a, Enum a, Num a) => a -> [b] -> [b]
irem indx hl = [ h | (h,i) <- zip hl [0..], i /= indx ]

-- | A simple product type for wrapping a pair of component values.
-- type Tuple a b = (a, b)
data Tuple a b = Tuple a b

-- List, PureScript does not provide sugar
type List a = [a]
range :: Int -> Int -> [Int]
range f t = [ f .. t ]
fromList :: [Number] ->  Array Number
fromList = A.fromList

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
divMod' :: Number -> Number -> Tuple Int Number
divMod' n d = Tuple f (n - toNumber f * d) where
    f = div' n d

-- | generalisation of 'mod' to any instance of Real
mod' :: Number -> Number -> Number
mod' n d = n - toNumber f * d where
    f = div' n d

-- | convert from String to Number and Int
numberFromString :: String -> Maybe Number
numberFromString  = readDouble where
  readDouble = readMaybe :: String -> Maybe Double
intFromString :: String -> Maybe Int
intFromString = readInt where
  readInt = readMaybe :: String -> Maybe Int


-- | generate a list of n normally distributed random values
-- | usinging the Box-Muller method and the random function
boxMuller :: RandomGen g => g -> (Number, Number, g)
boxMuller g = let
                (u1, g1) = random g
                (u2, g2) = random g1
                r = sqrt (-2.0 * log u1)
                t = 2.0 * pi * u2
                b1 = r * cos t
                b2 = r * sin t
             in ( b1, b2, g2 )
normals :: RandomGen g => Int -> g -> (Array Number, g)
normals n g = (ns, g') where
  nto :: Int
  nto = n `div` 2
  is :: List Int
  is = [0 .. n `div` 2]
  doNormals :: RandomGen g => (g, List Number) -> Int -> (g, List Number)
  doNormals (g, ns) _ = (g', n1:n2:ns) where
                        (n1,n2,g') = boxMuller g
  (g', ls) = foldl doNormals (g, []) is
  ns = A.fromList <<< L.take n $ ls


