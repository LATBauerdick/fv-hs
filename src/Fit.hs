-- file src/Fit.hs
module Fit ( fit, fit' ) where

import Types ( XMeas (..), HMeas (..), QMeas (..)
              , Prong (..), ABh0 (..), Chi2
             ,X3, C33, Q3
             )
import Coeff ( fvABh0, hv2q )
import Matrix ( inv, tr, sw, scalar, scale)
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

wght :: Double -> Chi2 -> Chi2 -- weight function with Temperature t
wght t chi2 = w where
  chi2cut = 9.0
  w = 1.0/(1.0 + exp ((chi2-chi2cut)/2.0/t))

fit :: XMeas -> [HMeas] -> Prong
fit v0 hl = pr where
  v = kfilter v0 hl
  pr = ksmooth v hl

fit' :: XMeas -> [HMeas] -> Prong -- fit with annealing function
fit' v0 hl = pr where
  v = kfilter v0 hl
  Prong _ _ _ cl = ksmooth v hl
  wl = fmap (wght 10.0) cl
  v' = kfilterW v0 hl wl --`debug` (printf "%8.1f" $head wl)
  Prong _ _ _ cl' = ksmoothW v' hl wl
  wl' = fmap (wght 1.0) cl'
  v'' = kfilterW v0 hl wl' --`debug` (printf "%8.1f" $head wl)
  pr = ksmoothW v'' hl wl'


kfilterW :: XMeas -> [HMeas] -> [Chi2] -> XMeas
kfilterW v0 hl wl = foldl kal v0 hl' where
  ff (HMeas h hh w2pt) w = HMeas h (scale w (inv hh)) w2pt -- `debug` (printf "%8.1f" w)
  hl' = zipWith ff hl wl
  kal :: XMeas -> HMeas -> XMeas
  kal (XMeas v vv) (HMeas h gg w2pt) = kalAdd v (inv vv) (HMeas h gg w2pt) v (hv2q h v) 1e6 0
--    `debug` ((showHMeas "kal: add helix " (HMeas h hh)) ++ (showXMeas "\nto vertex " (XMeas v vv)))

ksmoothW :: XMeas -> [HMeas] -> [Chi2] -> Prong
ksmoothW v hl wl = pr where
  ff :: HMeas -> Chi2 -> (QMeas, Chi2)
  ff (HMeas h hh w2pt) w = ksm (HMeas h (scale w (inv hh)) w2pt) v
  qml = zipWith ff hl wl
  (ql, chi2l) = unzip qml
  n = length chi2l
  pr = Prong n v ql chi2l

kfilter :: XMeas -> [HMeas] -> XMeas
kfilter = foldl kal where
  kal :: XMeas -> HMeas -> XMeas
  kal (XMeas v vv) (HMeas h hh w2pt) = kalAdd v (inv vv) (HMeas h (inv hh) w2pt) v (hv2q h v) 1e6 0
--    `debug` ((showHMeas "kal: add helix " (HMeas h hh)) ++ (showXMeas "\nto vertex " (XMeas v vv)))

kalAdd :: X3 -> C33 -> HMeas -> X3 -> Q3 -> Double -> Int -> XMeas
kalAdd v0 uu0 (HMeas h gg w2pt) ve qe chi20 iter = vm where
      goodEnough :: Double -> Double -> Int -> Bool
      goodEnough c0 c i = abs (c - c0) < chi2cut || i > iterMax where
          chi2cut = 0.5
          iterMax = 99 :: Int
      ABh0 aa bb h0 = fvABh0 ve qe
      aaT = tr aa
      bbT = tr bb
      ww = inv $ sw bb gg
      gb = gg - sw gg (sw bbT ww)
--      gb = gg - (sw (tr (gg * bb)) ww)
      uu = uu0 + sw aa gb
      cc = inv uu
      m =  h - h0
      v = cc * (uu0 * v0 + aaT * gb * m)
      dm = m - aa * v
      q = ww * bbT * gg * dm --`debug` ("--> " ++ prettyMatrix uu0 ++ prettyMatrix (sw aa gb) ++ prettyMatrix (cc*uu0*v0) ++ prettyMatrix (cc*aaT*gb*m))
--      ee = - ww * bbT * gg * aa * cc
--      dd = ww + (sw ww (sw bb (sw gg (sw aaT cc))))
      chi2 = scalar $ sw (dm - bb * q) gg + sw (v - v0) uu0
      ret = goodEnough chi20 chi2 iter
--               `debug`  showXMeas (printf ">iter %d, chi2%8.1g, dr=%8.3f, v(x,y,z) ->" (iter::Int) (chi2::Double) (dr::Double)) (XMeas v cc) where dr = sqrt (x^2+y^2+z^2) ; [x,y,z]=toList 3 (v-v0)
      vm  = if ret
               then XMeas v cc
               else kalAdd v0 uu0 (HMeas h gg w2pt) v q chi2 (iter +1)

ksmooth :: XMeas -> [HMeas] -> Prong
ksmooth v hl = pr where
  ƒ :: HMeas -> (QMeas, Chi2)
  ƒ (HMeas h hh w2pt) = ksm (HMeas h (inv hh) w2pt) v
  qml = map ƒ hl
  (ql, chi2l) = unzip qml
  n = length chi2l
  pr = Prong n v ql chi2l

-- kalman smooth: calculate 3-mom q at kalman filter vertex v
ksm :: HMeas -> XMeas -> (QMeas, Chi2)
ksm  (HMeas h gg w2pt) (XMeas x cc) = (QMeas q dd w2pt, chi2') -- `debug` ("≫" ++ show chi2)
  where
    ABh0 aa bb h0 = fvABh0 x (hv2q h x)
    aaT = tr aa
    bbT = tr bb
    ww = inv $ sw bb gg

    p = h - h0

    uu = inv cc
    q = ww * bbT * gg * (p - aa * x)
    ee = - cc * aaT * gg * bb * ww
    dd = ww + sw ee uu
    r = p - aa*x - bb*q -- `debug` (prettyMatrix (p - aa*x - bb*q ))
    gb = gg - sw gg ( sw bbT ww )
    uu' =  uu - sw aa gb
    cc' = inv uu'
    x' = cc' * (uu*x - aaT * gb * p)
    dx = x - x'
    chi2' = scalar $ sw dx uu' + sw r gg


