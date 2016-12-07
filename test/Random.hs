{-# LANGUAGE BangPatterns #-}

module Random ( rp ) where

import System.Random
import Data.Random.Normal

import Data.List ( foldl', unfoldr, mapAccumL, (!!) )

import Types ( XMeas (..), VHMeas (..), HMeas (..), Prong (..)
             , showXMeas, showMMeas, q2p, v3,l3,v5,l5 )
import Coeff ( invMass )
import Matrix ( toList, fromList, chol, scalar )
import Fit ( fit )

import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

-- return randomized XMeas vector according to its covariance matrix,
-- using -> 3 normal distributed random numbers

randV :: XMeas -> [Double] -> XMeas
randV (XMeas v vv) rnd = XMeas v' vv where
  v'  = v3 $ zipWith (+) (l3 v) (l3 (chol vv * v3 rnd))

randVH :: VHMeas HMeas -> [Double] -> (VHMeas HMeas, [Double])
randVH (VHMeas vm hl) rs = (VHMeas (randV vm rsv) hl', rsf) where
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

rp :: Int -> VHMeas HMeas -> IO ()
rp _cnt (VHMeas v hl) = do
  g <- newStdGen
  -- putStrLn $ showXMeas "input  " v
  -- let v' = randV v $ take 3 ( normals g )
  -- putStrLn $ showXMeas "smeared" v'

  let hist = histogram binSturges imd
  plot "invMass.png" hist

  let (VHMeas _ hl', rnd) = randVH (VHMeas v hl) (normals g)
  let Prong _ _ ql _ = fit v hl'
  showMMeas "Inv Mass " $ invMass $ map q2p ql
  let (VHMeas _ hl', rnd') = randVH (VHMeas v hl) rnd
  let Prong _ _ ql _ = fit v hl'
  showMMeas "Inv Mass " $ invMass $ map q2p ql
  let (VHMeas _ hl', rnd'') = randVH (VHMeas v hl) rnd'
  let Prong _ _ ql _ = fit v hl'
  showMMeas "Inv Mass " $ invMass $ map q2p ql

imd = [1494.8, 1491.9, 1492.1, 1494.2, 1494.5, 1497.0, 1499.0, 1494.9, 1494.8
  , 1497.3, 1494.8, 1499.3, 1492.5, 1492.2, 1502.3, 1489.6, 1494.8, 1494.0
  , 1489.5, 1497.5, 1499.5, 1497.1, 1490.2, 1496.5, 1494.8, 1494.6, 1498.0
  , 1493.8, 1490.5, 1490.4, 1496.4, 1486.7, 1494.8, 1490.0, 1500.1, 1497.1
  , 1495.8, 1494.9, 1490.4, 1491.5, 1494.8, 1501.7, 1495.2, 1492.5, 1494.8
  , 1489.1, 1491.0, 1491.5, 1494.8, 1501.5, 1494.9, 1499.9, 1494.8, 1494.1
  , 1493.3, 1497.5, 1494.8, 1495.7, 1496.1, 1498.7, 1494.8, 1495.2, 1503.6
  , 1489.5, 1494.8, 1494.7, 1491.8, 1498.7, 1494.8, 1490.7, 1501.8, 1502.3
  , 1494.8, 1494.5, 1493.6, 1495.4
  ]
