{-# LANGUAGE BangPatterns #-}

module Random ( rp ) where

import System.Random
import Data.Random.Normal

import Data.List ( foldl', unfoldr, mapAccumL, (!!) )

import Types ( XMeas (..), VHMeas (..), HMeas (..), Prong (..), MMeas (..)
             , showXMeas, showMMeas, q2p, v3,l3,v5,l5 )
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

randVH :: VHMeas HMeas -> [Double] -> (VHMeas HMeas, [Double])
randVH (VHMeas vm hl) rs = (VHMeas vm' hl', rsf) where
  vm' = vm --randV vm rsv
  (rsv, rsh) = splitAt 3 rs
  (rsf, hl') = mapAccumL randH rsh hl

randH :: [Double] -> HMeas -> ([Double], HMeas)
randH (r0:r1:r2:r3:r4:rs) (HMeas h hh w0) = (rs, HMeas h' hh w0) where
  h' = v5 $ zipWith (+) (l5 h) (l5 (chol hh * v5 [r0,r1,r2,r3,r4]))


-- create a list of histogram values (mapping XMeas -> Double using f)
-- as a function of randomized vertices
drs :: XMeas -> (XMeas -> Double) -> [Double] -> [Double]
drs v f rs = unfoldr (dr v f) rs

dr :: XMeas -> (XMeas -> Double) -> [Double] -> Maybe (Double, [Double])
dr v f (r0:r1:r2:rs) = Just (f (randV v [r0,r1,r2]), rs)

getVx :: XMeas -> Double
getVx (XMeas v _) = (!!) (l3 v) 2

drs' :: VHMeas HMeas -> (VHMeas HMeas -> Double) -> [Double] -> [Double]
drs' vh f rs = unfoldr (dr' vh f) rs

dr' :: VHMeas HMeas -> (VHMeas HMeas -> Double) -> [Double] -> Maybe (Double, [Double])
dr' vh f rs = Just (f (vh'), rs') where
  (vh', rs') = randVH vh rs

fm :: VHMeas HMeas -> Double
fm (VHMeas v hl) = m where
  Prong _ _ ql _ = fit v hl
  MMeas m _ = invMass $ map q2p ql

rp :: Int -> VHMeas HMeas -> IO ()
rp cnt (VHMeas v hl) = do
  g <- newStdGen
  putStrLn . display . process . take cnt . prep $  normals g

  let hist = histogram binSturges $ take cnt (drs' (VHMeas v hl) fm (normals g))
  plot "vz.png" hist
  let Prong _ _ ql _ = fit v hl
  showMMeas "Inv Mass " $ invMass $ map q2p ql

  -- putStrLn $ showXMeas "input  " v
  -- let v' = randV v $ take 3 ( normals g )
  -- putStrLn $ showXMeas "smeared" v'
  --
  -- let (VHMeas _ hl', rnd) = randVH (VHMeas v hl) (normals g)
  -- let Prong _ _ ql _ = fit v hl'
  -- showMMeas "Inv Mass " $ invMass $ map q2p ql
  -- let (VHMeas _ hl', rnd') = randVH (VHMeas v hl) rnd
  -- let Prong _ _ ql _ = fit v hl'
  -- showMMeas "Inv Mass " $ invMass $ map q2p ql
  -- let (VHMeas _ hl', rnd'') = randVH (VHMeas v hl) rnd'
  -- let Prong _ _ ql _ = fit v hl'
  -- showMMeas "Inv Mass " $ invMass $ map q2p ql

