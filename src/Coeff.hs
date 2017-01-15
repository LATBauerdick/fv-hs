-- file src/coeff.hs
--
module Coeff ( expand, qv2h, hv2q
             ) where

import Types ( M, V3, V5
              , Mom (..)
              , MMeas (..), HMeas (..), QMeas (..), PMeas (..), XMeas (..)
              , Jaco (..) )
import Matrix ( toList, fromList, fromList2 )

import Prelude
import Data.Fixed ( mod' )
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

hv2q :: V5 -> V3 -> V3
hv2q h v = q where
  [xx, yy, _] = toList 3 v
  r = sqrt $ xx*xx + yy*yy
  phi  = atan2 yy xx
  [w0, tl0, psi0, d0, z0] = toList 5 h
  xi = mod' (psi0 - phi + 2.0*pi) (2.0*pi)
  cxi = cos xi
  sxi = sin xi
  q = fromList 3 $
        if w0 /= 0 then
                  let
                    oow0 = 1.0/w0
                    gamma = atan r*cxi/(oow0-r*sxi)
                    in
                    [ w0, tl0, psi0 + gamma ]
                  else [ w0, tl0, psi0 ]



qv2h :: V3 -> V3 -> V5
qv2h q v = h where
  -- calculate and return helix parameters h = [w tl psi0 d0 z0]
  -- from v = [vx vy vz] and q = [w tl psi]
  [xx, yy, z] = toList 3 v
  r = sqrt $ xx*xx + yy*yy
  phi  = atan2 yy xx
  [w, tl, psi] = toList 3 q
  -- some more derived quantities
  xi = mod' (psi - phi + 2.0*pi) (2.0*pi)
  cxi = cos xi
  sxi = sin xi
  h = fromList 5 $
        if (w /= 0) then
                  let
                    oow = 1.0/w
                    gamma = atan r*cxi/(oow-r*sxi)
                    d0 = oow - (oow-r*cxi)/(cos gamma)
                    z0 = z - gamma*oow*tl
                    in
                    [ w, tl, psi - gamma, d0, z0 ]
                  else [w, tl, psi, r*sxi, z]

vmqm2hm :: XMeas -> QMeas -> HMeas
vmqm2hm (XMeas v _) (QMeas q _ w2pt) = HMeas h hh w2pt where
  h = qv2h q v
  hh = undefined

expand :: M -> M -> Jaco
expand v q = Jaco aa bb h0 where
  [xx, yy, z] = toList 3 v
  r = sqrt $ xx*xx + yy*yy
  phi  = atan2 yy xx
  [w, tl, psi] = toList 3 q
  -- some more derived quantities
  xi = mod' (psi - phi + 2.0*pi) (2.0*pi)
  cxi = cos xi
  sxi = sin xi
  oow = 1.0 / w
  rw  = r * w

  gamma = atan $ r*cxi/(oow - r*sxi)
  sg    = sin gamma
  cg    = cos gamma

  -- calculate transformed quantities
  psi0  = psi - gamma
  d0    = oow - (oow - r*sxi)/cg
  z0    = z - tl*gamma/w

  -- calc Jacobian
  [drdx, drdy, rdxidx, rdxidy] =
    if r /= 0 then [xx/r, yy/r, yy/r, -xx/r]
             else [0, 0, 0, 0]
  dgdvar0 =    1.0/(1.0 + rw*rw - 2.0*rw*sxi)
  dgdx    =    dgdvar0*(w*cxi*drdx + w*(rw - sxi)*rdxidx)
  dgdy    =    dgdvar0*(w*cxi*drdy + w*(rw - sxi)*rdxidy)
  dgdw    =    dgdvar0*r*cxi
  dgdpsi  =    dgdvar0*rw*(rw - sxi)

  --  fill matrix:
  -- d w / d r, d phi, d z
  [ a11, a12, a13 ]  =  [ 0.0, 0, 0 ]
  -- d tl / d x, d y, d z
  [ a21, a22, a23 ]  =  [ 0.0, 0, 0 ]
  -- d psi0 / d x, d y, d z
  [ a31, a32, a33 ]  =  [ -dgdx, -dgdy, 0 ]
  -- d d0 / d x, d y, d z
  [ a41, a42, a43 ]  =  [ cxi*rdxidx/cg + sxi*drdx/cg
                            - (oow - r*sxi)*sg*dgdx/cg/cg,
                          cxi*rdxidy/cg + sxi*drdy/cg
                            - (oow - r*sxi)*sg*dgdy/cg/cg,
                           0
                         ]
    -- d z0 / d x, d y, d z
  [ a51, a52, a53 ]  =  [ -tl/w*dgdx, -tl/w*dgdy, 1.0 ]

  -- B
  -- d w / d w, d tl, d psi
  [ b11, b12, b13 ]  =  [ 1.0, 0, 0 ]
  -- d tl / d w, d tl, d psi
  [ b21, b22, b23 ]  =  [ 0, 1.0, 0 ]
  -- d psi0 / d w, d tl, d psi
  [ b31, b32, b33 ]  =  [ -dgdw, 0, 1.0 - dgdpsi ]
  -- d d0 / d w, d tl, d psi
  [ b41, b42, b43 ]  =  [ -oow*oow*(1.0 - 1.0/cg)
                          - (oow - r*sxi)*sg*dgdw/cg/cg,
                          0,
                          r*cxi/cg - (oow - r*sxi)*sg*dgdpsi/cg/cg ]
  -- d z0 / d w, d tl, d psi
  [ b51, b52, b53 ]  =  [ -tl/w*(dgdw - gamma/w),
                          -gamma/w,
                          -tl/w*dgdpsi
                        ]
  [ v01, v02, v03 ] = toList 3 v
  [ q01, q02, q03 ] = [w, tl, psi]
  h0 = fromList 5 [
      0,
      0,
      psi0 - a31*v01 - a32*v02 - b31*q01 - b33*q03,
      d0 - a41*v01 - a42*v02 - b41*q01 - b43*q03,
      z0 - a51*v01 - a52*v02 - a53*v03 - b51*q01 - b52*q02 - b53*q03]
  aa = fromList2 5 3 [a11,a12,a13,a21,a22,a23,a31,a32,a33,a41,a42,a43,a51,a52,a53]
  bb = fromList2 5 3 [b11,b12,b13,b21,b22,b23,b31,b32,b33,b41,b42,b43,b51,b52,b53]
  -- aaT = tr aaTT `debug` ( "v0 --->> " ++ (show v) ++
  --                       "q0 --->> " ++ (show q) ++
  --                       "aa --->> " ++ prettyMatrix aaTT ++
  --                       "bb --->> " ++ prettyMatrix bb ++
  --                       "h0 --->> " ++ (show h0)
  --                         )
  -- aa = tr aaT

