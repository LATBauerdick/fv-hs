module Stuff where

import Prelude
import Math ( sqrt )
import Data.Ord (signum)
import Data.String ( takeWhile, dropWhile, toCharArray, fromCharArray, split, Pattern (..) )
import Data.Char (toCharCode)
import Data.List ( List(..), (:))
import Data.Array ( unsafeIndex, range, length, take, concat ) as A
import Data.Unfoldable ( replicateA )
import Data.Tuple ( Tuple(..), fst, snd )
import Data.Maybe ( Maybe(..), fromMaybe', fromJust )
import Data.Foldable ( class Foldable, foldr, sum )
import Partial.Unsafe (unsafePartial, unsafePartialBecause, unsafeCrashWith)
import Unsafe.Coerce ( unsafeCoerce ) as Unsafe.Coerce
import Data.List ( fromFoldable )
import Data.Int ( round, toNumber, floor )
import Text.Format ( format, precision, width )
import Data.Enum ( class Enum )
import Control.MonadZero ( guard )
import Control.Monad.Eff.Unsafe (unsafePerformEff)
import Math ( log, sqrt, pi, sin, cos ) as Math

-- | Returns `True` for any Unicode space character, and the control
-- | characters `\t`, `\n`, `\r`, `\f`, `\v`.
-- |
-- | `isSpace` includes non-breaking space.
-- avoiding to include unicode which increases make time
isSpace :: Char -> Boolean
-- The magic 0x377 used in the code below isn't really that magical. As of
-- 2014, all the codepoints at or below 0x377 have been assigned, so we
-- shouldn't have to worry about any new spaces appearing below there.
isSpace c = if uc <= 0x337
               then uc == 32 || (uc >= 9 && uc <= 13) || uc == 0xa0
               else false
  where
    uc :: Int
    uc = toCharCode c

-- | square of a number
sqr :: Number -> Number
sqr a = a*a

-- | unsafe index to Array
uidx :: forall a. Array a -> Int -> a
uidx = unsafePartial A.unsafeIndex

uJust :: forall a. Maybe a -> a
uJust = unsafePartial $ fromJust

-- | filter list of objects given list of indices in [a]
-- | return list with only those b that have  indices that  are in rng [a]
iflt :: forall a. Array Int -> Array a  -> Array a
iflt rng hl = do
  i <- rng
  pure $ uidx hl i

-- | remove element at index
irem :: forall a. Int -> Array a -> Array a
irem indx hl = do
  i <- A.range 0 ((A.length hl)-1)
  guard $ i /= indx
  pure $ uidx hl i


-- | round to 3 decimal
roundDec :: Number -> Number
roundDec x = (toNumber (round ( 1000.0 * x )))/1000.0

to0fix :: Number -> String
to0fix = format (width 4 <> precision 0)
to1fix :: Number -> String
to1fix = format (width 6 <> precision 1)
to2fix :: Number -> String
to2fix = format (width 7 <> precision 2)
to3fix :: Number -> String
to3fix = format (width 8 <> precision 3)
to5fix :: Number -> String
to5fix = format (width 10 <> precision 5)

-- | simultaneous 'quot' and 'rem'
quotRem :: Int -> Int -> (Tuple Int Int)
--quotRem             :: a -> a -> (a,a)
quotRem n d = if signum r == - signum d
                 then (Tuple (q+1) (r-d))
                 else qr
  where qr = divMod n d
        q = fst qr
        r = snd qr

-- | generate a list of n normally distributed random values
-- | usinging the Box-Muller method and the random function
boxMuller :: forall e. Eff (random :: RANDOM | e) (Array Number)
boxMuller = do
              u1 <- random
              u2 <- random
              let r = Math.sqrt (-2.0 * Math.log u1)
                  t = 2.0 * Math.pi * u2
                  b1 = r * Math.cos t
                  b2 = r * Math.sin t
              pure $ [ b1, b2 ]
normals :: forall e. Int -> Eff (random :: RANDOM | e) (Array Number)
normals n = do
  ls <- replicateA ((n+1)/2) $ boxMuller
  pure $ A.take n $ A.concat ls

-- | Calculate mean and standard deviation
stats :: Array Number -> Tuple Number Number
stats xs = Tuple mean stddev where
  n = toNumber (A.length xs)
  mean = sum xs / n
  stddev = sqrt $ sum (map (\v -> sqr (v-mean)) xs) / n


