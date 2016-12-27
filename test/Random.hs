{-# LANGUAGE BangPatterns #-}

module Random ( doRandom ) where

import Control.Parallel.Strategies
import System.Random
import Data.Random.Normal
import Statistics.Sample ( meanVariance )

import qualified Data.Vector.Unboxed as V ( Vector, fromList, toList )

import Data.List ( foldl', unfoldr, mapAccumL, (!!) )

import Types ( XMeas (..), VHMeas (..), HMeas (..), Prong (..), MMeas (..)
             , invMass, q2p, v3,l3,v5,l5 )
import Matrix ( toList, fromList, chol, scalar )
import Fit ( fit, fitw )

import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

doRandom :: Int -> VHMeas -> IO ()
doRandom cnt vm = do
  let Prong _ _ ql _ = fitw vm
  putStrLn $ "Fit Mass  " ++ (show . invMass . map q2p) ql

  g <- newStdGen
  let hf :: V.Vector Double
      hf = V.fromList
          . withStrategy (parBuffer 100 rdeepseq)
          . map fitm
          . take cnt
          . gen vm . normals $ g
      (mean, var) = meanVariance hf
  putStrLn $ "Mean Mass " ++ show (MMeas mean (sqrt var))
  let hist = histogram binSturges (V.toList hf)
  _ <- plot "invMass.png" hist
  return ()

gen :: VHMeas -> [Double] -> [VHMeas]
gen v rs = v' : gen v rs' where
  (v', rs') = randVH v rs

-- calc fitted invariant mass of VHMeas
fitm :: VHMeas -> Double
fitm vm = m where
  Prong _ _ ql _ = fitw vm
  (MMeas m _) = invMass . map q2p $ ql

-- randomize the helices in the supplied VHMeas
-- and return randomized VHMeas and remaining randoms list
randVH :: VHMeas -> [Double] -> (VHMeas, [Double])
randVH (VHMeas v hl) rs = (VHMeas v hl', rs') where
  (rs', hl') = mapAccumL randH rs hl

-- randomize a single helix parameters measurement, based on the cov matrix
-- return randomized helix and "remaining" random numbers
randH :: [Double] -> HMeas -> ([Double], HMeas)
randH (r0:r1:r2:r3:r4:rs) (HMeas h hh w0) = (rs, HMeas h' hh w0) where
  h' = v5 $ zipWith (+) (l5 h) (l5 (chol hh * v5 [r0,r1,r2,r3,r4]))

{-
-- return randomized XMeas vector according to its covariance matrix,
-- using -> 3 normal distributed random numbers
randV :: XMeas -> [Double] -> XMeas
randV (XMeas v vv) rnd = XMeas v' vv where
  v'  = v3 $ zipWith (+) (l3 v) (l3 (chol vv * v3 rnd))

-}

