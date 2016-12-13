-- file src/Fit.hs
module Fit ( fit, fitw ) where

import Types (  XMeas (..), HMeas (..), QMeas (..), VHMeas (..)
              , Prong (..), Jaco (..), Chi2
              , X3, C33, Q3
             )
import qualified Coeff ( expand, hv2q )
import Matrix ( inv, tr, sw, scalar, scale)
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

wght :: Double -> Chi2 -> Double -- weight function with Temperature t
wght t chi2 = w where
  chi2cut = 9.0
  w = 1.0/(1.0 + exp ((chi2-chi2cut)/2.0/t))

fitw :: XMeas -> [HMeas] -> Prong -- fit with annealing function
fitw v0 hl = pr where
  v = kfilter' v0 hl
  Prong _ _ _ cl = ksmooth' v hl
  wl = fmap (wght 10.0) cl
  v' = kfilterW v0 hl wl --`debug` (printf "%8.1f" $head wl)
  Prong _ _ _ cl' = ksmoothW v' hl wl
  wl' = fmap (wght 1.0) cl'
  v'' = kfilterW v0 hl wl' --`debug` (printf "%8.1f" $head wl)
  pr  = ksmoothW v'' hl wl'


kfilterW :: XMeas -> [HMeas] -> [Chi2] -> XMeas
kfilterW v0 hl wl = foldl kal v0 hl' where
  ff (HMeas h hh w0) w = HMeas h (scale w (inv hh)) w0 -- `debug` (printf "%8.1f" w)
  hl' = zipWith ff hl wl
  kal :: XMeas -> HMeas -> XMeas
  kal (XMeas v vv) (HMeas h gg w0) = kalAdd' v (inv vv) (HMeas h gg w0) v (Coeff.hv2q h v) 1e6 0
--    `debug` ((showHMeas "kal: add helix " (HMeas h hh)) ++ (showXMeas "\nto vertex " (XMeas v vv)))

ksmoothW :: XMeas -> [HMeas] -> [Chi2] -> Prong
ksmoothW v hl wl = pr where
  ff :: HMeas -> Chi2 -> (QMeas, Chi2)
  ff (HMeas h hh w0) w = ksm' (HMeas h (scale w (inv hh)) w0) v
  qml = zipWith ff hl wl
  (ql, chi2l) = unzip qml
  n   = length chi2l
  pr  = Prong n v ql chi2l

kfilter' :: XMeas -> [HMeas] -> XMeas
kfilter' = foldl kal where
  kal :: XMeas -> HMeas -> XMeas
  kal (XMeas v vv) (HMeas h hh w0) = kalAdd' v (inv vv) (HMeas h (inv hh) w0) v (Coeff.hv2q h v) 1e6 0

kalAdd' :: X3 -> C33 -> HMeas -> X3 -> Q3 -> Double -> Int -> XMeas
kalAdd' v0 uu0 (HMeas h gg w0) ve qe chi2_0 iter = vm where
  Jaco aa bb h0 = Coeff.expand ve qe
  aaT  = tr aa
  bbT  = tr bb
  ww   = inv (sw bb gg)
  gb   = gg - sw gg (sw bbT ww)
  uu   = uu0 + sw aa gb
  cc   = inv uu
  m    =  h - h0
  v    = cc * (uu0 * v0 + aaT * gb * m)
  dm   = m - aa * v
  q    = ww * bbT * gg * dm
  chi2 = scalar $ sw (dm - bb * q) gg + sw (v - v0) uu0
  vm   = if goodEnough chi2_0 chi2 iter
            then XMeas v cc
            else kalAdd' v0 uu0 (HMeas h gg w0) v q chi2 (iter +1)

goodEnough :: Double -> Double -> Int -> Bool
goodEnough c0 c i = abs (c - c0) < chi2cut || i > iterMax where
  chi2cut = 0.5
  iterMax = 99 :: Int

ksmooth' :: XMeas -> [HMeas] -> Prong
ksmooth' v hl = pr where
  ƒ :: HMeas -> (QMeas, Chi2)
  ƒ (HMeas h hh w0) = ksm' (HMeas h (inv hh) w0) v
  qml = map ƒ hl
  (ql, chi2l) = unzip qml
  n   = length chi2l
  pr  = Prong n v ql chi2l

-- kalman smooth: calculate 3-mom q at kalman filter vertex v
ksm' :: HMeas -> XMeas -> (QMeas, Chi2)
ksm'  (HMeas h gg w0) (XMeas x cc) = (QMeas q dd w0, chi2') -- `debug` ("≫" ++ show chi2)
  where
    Jaco aa bb h0 = Coeff.expand x (Coeff.hv2q h x)
    aaT  = tr aa-- (aa ^+)
    bbT  = tr bb
    ww   = inv (sw bb gg)
    p    = h - h0
    uu   = inv cc
    q    = ww * bbT * gg * (p - aa * x)
    ee   = - cc * aaT * gg * bb * ww
    dd   = ww + sw ee uu
    r    = p - aa*x - bb*q
    gb   = gg - sw gg (sw bbT ww)
    uu'  =  uu - sw aa gb
    cc'  = inv uu'
    x'   = cc' * (uu*x - aaT * gb * p)
    dx   = x - x'
    chi2' = scalar (sw dx uu' + sw r gg)

{-
  data Prong = Prong N XMeas [QMeas] [Chi2] ...
  data VHMeas = VHMeas XMeas [HMeas] ...
  instance Monoid VHMeas where ...
-}
fit :: VHMeas -> Prong
fit = ksmooth . kFilter

kFilter :: VHMeas -> VHMeas
ksmooth :: VHMeas -> Prong
kFilter (VHMeas x ps) = VHMeas (foldl kAdd x ps) ps

kAdd :: XMeas -> HMeas -> XMeas
kAdd (XMeas v vv) (HMeas h hh w0) = kAdd' x_km1 p_k x_e q_e 1e6 0 where
  x_km1 = XMeas v (inv vv)
  p_k   = HMeas h (inv hh) w0
  x_e   = v
  q_e   = Coeff.hv2q h v

kAdd' :: XMeas -> HMeas -> X3 -> Q3 -> Double -> Int -> XMeas
kAdd' (XMeas v0 uu0) (HMeas h gg w0) x_e q_e 𝜒2_0 iter = x_k where
  Jaco aa bb h0 = Coeff.expand x_e q_e
  aaT   = tr aa; bbT = tr bb
  ww    = inv (sw bb gg)
  gb    = gg - sw gg (sw bbT ww)
  uu    = uu0 + sw aa gb; cc = inv uu
  m     = h - h0
  v     = cc * (uu0 * v0 + aaT * gb * m)
  dm    = m - aa * v
  q     = ww * bbT * gg * dm
  𝜒2    = scalar $ sw (dm - bb * q) gg + sw (v - v0) uu0
  x_k   = if goodEnough 𝜒2_0 𝜒2 iter
            then XMeas v cc
            else kAdd' (XMeas v0 uu0) (HMeas h gg w0) v q 𝜒2 (iter+1)

ksmooth (VHMeas v hl) = Prong (length ql) v ql chi2l where
  (ql, chi2l) = unzip $ map (ksm v) hl

-- kalman smooth: calculate 3-mom q and chi2 at kalman filter vertex
ksm :: XMeas -> HMeas -> (QMeas, Chi2)
ksm (XMeas x cc) (HMeas h hh w0) = (QMeas q dd w0, chi2) where
    Jaco aa bb h0 = Coeff.expand x (Coeff.hv2q h x)
    aaT  = tr aa
    bbT  = tr bb
    gg   = inv hh
    ww   = inv (sw bb gg)
    p    = h - h0
    uu   = inv cc
    q    = ww * bbT * gg * (p - aa * x)
    ee   = - cc * aaT * gg * bb * ww
    dd   = ww + sw ee uu
    r    = p - aa*x - bb*q
    gb   = gg - sw gg (sw bbT ww)
    uu'  =  uu - sw aa gb
    cc'  = inv uu'
    x'   = cc' * (uu*x - aaT * gb * p)
    dx   = x - x'
    chi2 = scalar (sw dx uu' + sw r gg)

