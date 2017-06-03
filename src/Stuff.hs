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
  , to0fix, to1fix, to2fix, to3fix, to5fix
  , sqr
  )  where

import Prelude
import qualified Data.Vector.Unboxed as A
  ( Vector, length, fromList, toList, unsafeIndex, create, replicate
  , singleton, map, foldl, zipWith )
import Data.Maybe (  Maybe (..), fromJust )
import Text.Printf ( printf )

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
data Tuple a b = Tuple a b

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

