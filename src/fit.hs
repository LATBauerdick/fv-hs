-- file src/Fit.hs
module Fit ( fit ) where

import Types ( XMeas (..), HMeas (..), QMeas (..)
              , Prong (..), ABh0 (..), Chi2
             ,X3, C33, Q3, H5, C55
             , showXMeas, showHMeas )
import Coeff ( fvABh0, qv2h, hv2q)
import Matrix ( inv, tr, sw, sub, sub2, scalar, toList)
import Text.Printf
import Data.Matrix ( prettyMatrix )
import Debug.Trace ( trace )
debug = flip trace


fit :: XMeas -> [HMeas] -> Prong
fit v0 hl = pr where
  v = kfilter v0 hl
  pr = ksmooth v hl

kfilter :: XMeas -> [HMeas] -> XMeas
kfilter v0 hl = foldl kal v0 hl where
  kal :: XMeas -> HMeas -> XMeas
  kal (XMeas v vv) (HMeas h hh) = kalAdd v (inv vv) h (inv hh) v (hv2q h) 1e6 0
--    `debug` ((showHMeas "kal: add helix " (HMeas h hh)) ++ (showXMeas "\nto vertex " (XMeas v vv)))

kalAdd :: X3 -> C33 -> H5 -> C55 -> X3 -> Q3 -> Double -> Int -> XMeas
kalAdd v0 uu0 h gg ve qe chi20 iter = vm where
      goodEnough :: Double -> Double -> Int -> Bool
      goodEnough chi20 chi2 iter =
        (abs $ chi2 - chi20) < chi2cut || iter > iterMax where
          chi2cut = 0.5
          iterMax = 99 :: Int
      ABh0 aa bb h0 = fvABh0 ve qe
      aaT = tr aa
      bbT = tr bb
      ww = inv $ sw bb gg
      gb = gg - (sw gg (sw bbT ww))
--      gb = gg - (sw (tr (gg * bb)) ww)
      uu = uu0 + (sw aa gb)
      cc = inv uu
      m =  h - h0
      v = cc * (uu0 * v0 + aaT * gb * m)
      dm = m - aa * v
      q = ww * bbT * gg * dm --`debug` ("--> " ++ prettyMatrix uu0 ++ prettyMatrix (sw aa gb) ++ prettyMatrix (cc*uu0*v0) ++ prettyMatrix (cc*aaT*gb*m))
--      ee = - ww * bbT * gg * aa * cc
--      dd = ww + (sw ww (sw bb (sw gg (sw aaT cc))))
      chi2 = scalar $ (sw (dm - bb * q) gg) + (sw (v - v0) uu0)
      ret = (goodEnough chi20 chi2 iter)
--               `debug`  showXMeas (printf ">iter %d, chi2%8.1g, dr=%8.3f, v(x,y,z) ->" (iter::Int) (chi2::Double) (dr::Double)) (XMeas v cc) where dr = sqrt (x^2+y^2+z^2) ; [x,y,z]=toList 3 (v-v0)
      vm  = if ret
               then XMeas v cc
               else kalAdd v0 uu0 h gg v q chi2 (iter +1)

ksmooth :: XMeas -> [HMeas] -> Prong
ksmooth v hl = pr where
  ƒ :: HMeas -> (QMeas, Chi2)
  ƒ h = ksm h v
  qml = map ƒ hl
  (ql, chi2l) = unzip qml
  n = length chi2l
  pr = Prong n v ql chi2l

-- kalman smooth: calculate 3-mom q at kalman filter vertex v
ksm :: HMeas -> XMeas -> (QMeas, Chi2)
ksm (HMeas h hh) (XMeas v vv) = ((QMeas q dd), chi2) -- `debug` ("≫" ++ show chi2)
  where
    ABh0 aa bb h0 = fvABh0 v (hv2q h)
    aaT = tr aa
    bbT = tr bb
    gg = inv hh
    ww = inv $ sw bb gg
    dm = h - h0 - aa * v
    q = ww * bbT * gg * dm
    dd = ww + (sw ww (sw bb (sw gg (sw aaT vv))))
    chi2 = scalar $ (sw (dm - bb * q) gg)

