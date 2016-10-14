-- file src/Fit.hs
module Fit ( fit ) where

import Types ( XMeas (..), HMeas (..), QMeas (..), Prong (..), ABh0 (..) )
import Coeff ( fvABh0 )
import Matrix ( inv, tr, sw, sub, sub2, scalar )
import Debug.Trace ( trace )
debug = flip trace


qMeas :: HMeas -> QMeas
qMeas (HMeas h ch) = q where
  q = QMeas (sub 3 h) (sub2 3 ch)

fit :: XMeas -> [HMeas] -> Prong
fit v0 hl = pr where
  v = kfilter v0 hl
  pr = smooth v hl

kfilter :: XMeas -> [HMeas] -> XMeas
kfilter v0 hl = v where
  v = foldr kal v0 hl

kal :: HMeas -> XMeas -> XMeas
kal (HMeas h hh) (XMeas v0 vv0) = v' `debug` ("." ++ show chi2)
  where
    ABh0 aa bb h0 = fvABh0 v0 h
    aaT = tr aa
    bbT = tr bb
    gg  = inv hh
    uu0 = inv vv0
    ww = inv $ sw bb gg
    gb = gg - (sw gg (sw bbT ww))
    uu = uu0 + (sw aa gb)
    cc = inv uu
    m =  h - h0
    v = cc * (uu0 * v0 + aaT * gb * m)
    dm = m - aa * v
    q = ww * bbT * gg * dm
    ee = - ww * bbT * gg * aa * cc
    dd = ww + (sw ww (sw bb (sw gg (sw aaT cc))))
    chi2 = scalar $ (sw (dm - bb * q) gg) + (sw (v - v0) uu0)
    v' = XMeas v cc


-- kalman smooth: calculate helices hl at kalman filter vertex v
smooth :: XMeas -> [HMeas] -> Prong
smooth v hl =
  Prong 6 v (map qMeas hl) [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
