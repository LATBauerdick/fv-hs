{-# LANGUAGE BangPatterns #-}

module Random ( doRandom ) where

import System.Random
import Data.Random.Normal
import Statistics.Sample ( meanVariance )

import qualified Data.Vector.Unboxed as V ( Vector, fromListN, toList )

import Data.List ( foldl', unfoldr, mapAccumL, (!!) )

import Types ( XMeas (..), VHMeas (..), HMeas (..), Prong (..), MMeas (..)
             , showXMeas, q2p, v3,l3,v5,l5 )
import Coeff ( invMass )
import Matrix ( toList, fromList, chol, scalar )
import Fit ( fit )

import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

-- Function to process our random sequence
process :: [(Double, Double)] -> (Int, Int)
process = foldl' sumInCircle (0, 0)

-- Function to process a running value and a random value, producing a new running value.
sumInCircle :: (Int, Int) -> (Double, Double) -> (Int, Int)
sumInCircle (!ins, !total) (x, y) = (ins + if x*x + y*y < 1.0 then 1 else 0,
                               total + 1)

-- Function to display the running value.
display:: (Int, Int) -> String
display (heads, coins) = "π = "  ++ (show $ 4.0 * fromIntegral heads / fromIntegral coins)

-- function to prepare the random sequence for processing
prep :: [Double] -> [(Double, Double)]
prep (a:b:r) = (a,b):prep r


-- return randomized XMeas vector according to its covariance matrix,
-- using -> 3 normal distributed random numbers
randV :: XMeas -> [Double] -> XMeas
randV (XMeas v vv) rnd = XMeas v' vv where
  v'  = v3 $ zipWith (+) (l3 v) (l3 (chol vv * v3 rnd))

-- randomize a single helix parameters measurment, based on the cov matrix
-- return randomized helix and "remaining" random numbers
randH :: [Double] -> HMeas -> ([Double], HMeas)
randH (r0:r1:r2:r3:r4:rs) (HMeas h hh w0) = (rs, HMeas h' hh w0) where
  h' = v5 $ zipWith (+) (l5 h) (l5 (chol hh * v5 [r0,r1,r2,r3,r4]))

{-
VHMeas v hl <- hSlurp thisFile
doRandom 1000 (VHMeas v (hFilter hl [0,2,3,4,5])) -}
doRandom :: Int -> VHMeas -> IO ()
doRandom cnt vm = do
  let Prong _ _ ql _ = fit vm
  putStrLn $ "Fit Mass  " ++ (show . invMass . map q2p) ql

  g <- newStdGen
  let hf :: V.Vector Double
      hf = V.fromListN cnt $ unfoldr (randomize vm fitMass) . normals $ g
      (mean, var) = meanVariance hf
  putStrLn $ "Mean Mass " ++ show (MMeas mean (sqrt var))
  let hist = histogram binSturges (V.toList hf)
  _ <- plot "invMass.png" hist
  return ()

randomize :: VHMeas -> (VHMeas -> Double) -> [Double] -> Maybe (Double, [Double])
randomize vh f rs = Just (f vh', rs') where
  (vh', rs') = randVH vh rs

-- calc fitted invariant mass of VHMeas
fitMass :: VHMeas -> Double
fitMass vm = m where
  Prong _ _ ql _ = fit vm
  (MMeas m _) = invMass $ map q2p ql

-- randomize the helices in the supplied VHMeas
-- taking normal-distributed random numbers off the (infinite) list of Doubles
-- and return randomized VHMeas and tail of randoms list
randVH :: VHMeas -> [Double] -> (VHMeas, [Double])
randVH (VHMeas v hl) rs = (VHMeas v hl', rs') where
  (rs', hl') = mapAccumL randH rs hl

-- create a list of values for histogramming, randomizing the supplied VHMeas
-- using the list of normaly distributed random numbers
-- mapping each randomized VHMeas -> Double using f
-- also return "remaining" list of random numbers
histVals :: VHMeas -> (VHMeas -> Double) -> [Double] -> [Double]
histVals vh f rs = unfoldr (randomize vh f) rs

{-
  putStrLn . display . process . take cnt . prep $  normals g
  putStrLn $ showXMeas "input  " v
  let v' = randV v $ take 3 ( normals g )
  putStrLn $ showXMeas "smeared" v'

  let (VHMeas _ hl', rnd) = randVH (VHMeas v hl) (normals g)
  let Prong _ _ ql _ = fit v hl'
  showMMeas "Inv Mass " $ invMass $ map q2p ql
  let (VHMeas _ hl', rnd') = randVH (VHMeas v hl) rnd
  let Prong _ _ ql _ = fit v hl'
  showMMeas "Inv Mass " $ invMass $ map q2p ql
  let (VHMeas _ hl', rnd'') = randVH (VHMeas v hl) rnd'
  let Prong _ _ ql _ = fit v hl'
  showMMeas "Inv Mass " $ invMass $ map q2p ql
  return ()

-}

