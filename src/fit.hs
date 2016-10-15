-- file src/Fit.hs
module Fit ( fit ) where

import Types ( XMeas (..), HMeas (..), QMeas (..)
             , Prong (..), ABh0 (..), Chi2 )
import Coeff ( fvABh0, fvh )
import Matrix ( inv, tr, sw, sub, sub2, scalar )
import Debug.Trace ( trace )
debug = flip trace


fit :: XMeas -> [HMeas] -> Prong
fit v0 hl = pr where
  v = kfilter v0 hl
  pr = ksmooth v hl

kfilter :: XMeas -> [HMeas] -> XMeas
kfilter v0 hl = v where
  v = foldr kal v0 hl

ksmooth :: XMeas -> [HMeas] -> Prong
ksmooth v hl = pr where
  ƒ :: HMeas -> (QMeas, Chi2)
  ƒ h = ksm h v
  qml = map ƒ hl
  (ql, chi2l) = unzip qml
  pr = Prong 6 v ql chi2l

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


-- kalman smooth: calculate 3-mom q at kalman filter vertex v
ksm :: HMeas -> XMeas -> (QMeas, Chi2)
ksm (HMeas h hh) (XMeas v vv) = qm `debug` ("^" ++ show chi2)
  where
    ABh0 aa bb h0 = fvABh0 v h
    aaT = tr aa
    m = h - h0
    dm = m - aa * v
    gg = inv hh
    ww = inv $ sw bb gg
    q = ww * (tr bb) * gg * dm
    cc = vv
    dd = (ww + (sw ww (sw bb (sw gg (sw aaT cc)))))
    -- use the simple method to calculate chi2
    chi2 = scalar $ sw (h - (fvh v q)) gg
    qm = (QMeas q dd, chi2)   --(sub 3 h) (sub2 3 hh)

