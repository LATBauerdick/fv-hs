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

module Test.Random ( testRandom ) where

import Prelude
-- import Control.Monad.Eff (Eff)
-- import Control.Monad.Eff.Console ( CONSOLE, log )
-- import Control.Monad.Eff.Random ( RANDOM )
-- import Data.List.Lazy ( replicateM )
import Data.Vector.Unboxed as A
  ( zipWith, replicateM, fromList )
-- import Data.Traversable ( for )
-- import Data.Tuple ( Tuple (..) )
import Statistics.Sample ( meanVariance )
import Statistics.Sample ( meanVariance )
import System.Random ( RandomGen, newStdGen )

import Data.Cov ( chol, fromArray, toArray, (*.), Vec5, Cov5 )
import FV.Fit ( fit )
import FV.Types ( VHMeas(..), HMeas(..), MMeas(..)
  , invMass, fromQMeas, fitMomenta )
import Stuff

{-- import qualified Graphics.Gnuplot.Frame.OptionSet as Opts --}
{-- import Graphics.Histogram --}

-- | Randomizable TypeClass to provide randomize method
-- | for MC smearing of a measurement
class Randomizable a where
  randomize :: RandomGen g => a -> g -> (a, g)
-- | randomize a single helix parameters measurement, based on the cov matrix
-- | return randomized helix
instance Randomizable HMeas where
  randomize (HMeas h hh w0) g = (HMeas h' hh w0, g') where
    (rs, g') = normals 5 g
    r5 :: Vec5
    r5 = fromArray rs
    h' = fromArray $ A.zipWith (+) (toArray h) (toArray (chol hh *. r5))

-- | randomize a vertex measurement by randomizing each helix parameter measurement
-- | leaving the initial vertex untouched
instance Randomizable VHMeas where
  randomize (VHMeas { vertex= v, helices= hl}) g = 
    (VHMeas { vertex= v, helices= hl' }, g') where
      doit :: RandomGen g => HMeas -> (List HMeas, g) -> (List HMeas, g)
      doit h (hs, g) = (hs', g') where
        (h', g') = randomize h g
        hs' = h' : hs
      (hl', g') = foldr doit ([], g) hl

-- calc fitted invariant mass of VHMeas
fitm :: VHMeas -> Number
fitm vm = m where
  MMeas {m= m} = invMass <<< map fromQMeas <<< fitMomenta $ fit vm

testRandom ::Int
             -> VHMeas
             -> IO String
testRandom cnt vm = do
  g <- newStdGen
  let ls :: [Int]
      ls = [ 0 .. (cnt-1) ]
      (ms, _) = foldl doit ([], g) ls where
        doit :: RandomGen g => (List Number, g) -> Int -> (List Number, g)
        doit (ms, g) _ = (ms', g') where
          (vm', g') = randomize vm g
          m = fitm vm'
          ms' = m : ms
  let (m, dm2) = meanVariance $ A.fromList ms

  pure $ "Fit Mass  "
          <> (show <<< invMass <<< map fromQMeas
                              <<< fitMomenta <<< fit $ vm)
          <> "\nMean Mass " <> show (MMeas {m=m, dm=sqrt dm2})
  {-- let hist = histogram binSturges (V.toList hf) --}
  {-- _ <- plot "invMass.png" hist --}

