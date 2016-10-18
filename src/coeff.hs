-- file src/coeff.hs
--
module Coeff ( w2pt, fvABh0, fvh, h2p, h2q, q2p, invMass, mass
             ) where

import Types ( M, M33, V3, V5, MMeas (..), HMeas (..), QMeas (..), PMeas (..)
  , ABh0 (..) , w2pt, mπ )
import Matrix ( sub, sub2, toList, fromList, fromList2 )
import Data.Fixed ( mod' )

h2p :: HMeas -> PMeas
h2p hm = (q2p . h2q) hm

h2q :: HMeas -> QMeas -- just drop the d0, z0 part... fix!!!!
h2q (HMeas h ch) = QMeas q cq where
  q = (sub 3 h)
  cq = (sub2 3 ch)

q2p :: QMeas -> PMeas
q2p (QMeas q0 cq0) = (PMeas p0 cp0) where
  m = mπ
  [w,tl,psi0] = toList 3 q0
  sph  = sin psi0
  cph  = cos psi0
  pt   = w2pt / abs w
  px   = pt * cph
  py   = pt * sph
  pz   = pt * tl
  e = sqrt(px^2 + py^2 + pz^2 + m^2)
  ps = w2pt / w
  dpdk = ps*ps/w2pt
  [c11, c12, c13, _, c22, c23, _, _, c33] = toList 9 cq0
  xy = 2.0*ps*dpdk*cph*sph*c13
  sxx = (dpdk*cph)^2 * c11 + (ps*sph)^2 * c33 + xy
  sxy = cph*sph*(dpdk*dpdk*c11 - ps*ps*c33) +
           ps*dpdk*(sph*sph-cph*cph)*c13
  syy = (dpdk*sph)^2 * c11 + (ps*cph)^2 * c33 - xy
  sxz = dpdk*dpdk*cph*tl*c11 -
           ps*dpdk*(cph*c12-sph*tl*c13) -
           ps*ps*sph*c23
  syz = dpdk*dpdk*sph*tl*c11 -
           ps*dpdk*(sph*c12 + cph*tl*c13) +
           ps*ps*cph*c23
  szz = (dpdk*tl)^2 * c11 + ps*ps*c22 -
           2.0*ps*dpdk*tl*c12
  sxe = (px*sxx + py*sxy + pz*sxz)/e
  sye = (px*sxy + py*syy + pz*syz)/e
  sze = (px*sxz + py*syz + pz*szz)/e
  see = (px*px*sxx + py*py*syy + pz*pz*szz +
         2.0*(px*(py*sxy + pz*sxz) + py*pz*syz))/e/e

  p0 = fromList 4 [px,py,pz,e]
  cp0 = fromList2 4 4 [sxx, sxy, sxz, sxe, sxy, syy, syz, sye, sxz, syz, szz, sze, sxe, sye, sze, see]

invMass :: [PMeas] -> MMeas
invMass pl@(h:t) = mass ptot where
  ptot = foldr sumP h t


sumP :: PMeas -> PMeas -> PMeas
sumP (PMeas p1 cp1) (PMeas p2 cp2) = PMeas (p1+p2) (cp1 + cp2)

mass :: PMeas -> MMeas
mass (PMeas p cp) = mm  where
  [px,py,pz,e] = toList 4 p
  [c11, c12, c13, c14, _, c22, c23, c24, _, _, c33, c34, _, _, _, c44]
    = toList 16 cp
  m                = sqrt $ max (e^2-px^2-py^2-pz^2) 0
  sigm0            = px*c11*px + py*c22*py + pz*c33*pz + e*c44*e +
                       2.0*(px*(c12*py + c13*pz - c14*e)
                          + py*(c23*pz - c24*e)
                          - pz*c34*e)
  sigm             = ( sqrt $ max sigm0 0 ) / m
  mm               = MMeas m sigm

fvh :: V3 -> V3 -> V5
fvh v q = h where
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
                    in
                    [ w, tl, psi - gamma, oow - (oow - r*sxi)/( cos gamma), z - gamma/w*tl ]
                  else [w, tl, psi, r*sxi, z]

fvABh0 :: M -> M -> ABh0
fvABh0 v h = ABh0 aa bb h0 where
--  aa = fromList2 5 3 [1.0,0,0, 0,1.0,0, 0,0,1.0, 0,0,0, 0,0,0]
--  bb = fromList2 5 3 [1.0,0,0, 0,1.0,0, 0,0,1.0, 0,0,0, 0,0,0]
--  h0 = h
  [xx, yy, z] = toList 3 v
  r = sqrt $ xx*xx + yy*yy
  phi  = atan2 yy xx
  [w, tl, psi] = toList 3 h
  -- some more derived quantities
  xi = mod' (psi - phi + 2.0*pi) (2.0*pi)
  cxi = cos xi
  sxi = sin xi
  oow = 1.0 / w
  rw  = r * w

  gamma = atan $ r * cxi / (oow - r * sxi)
  sg    = sin gamma
  cg    = cos gamma

  -- calculate transformed quantities
  psi0  = psi - gamma
  d0    = oow - (oow - r * sxi) / cg
  z0    = z - tl * gamma / w

  -- calc Jacobian
  [drdx, drdy, rdxidx, rdxidy] =
    if (r /= 0) then [xx / r, yy / r, yy / r, - xx / r]
               else [0, 0, 0, 0]
  dgdvar0 =    1.0 / (1.0 + (rw * rw) - (2.0 * rw * sxi))
  dgdx    =    dgdvar0 * ((w * cxi * drdx) + (w * (rw - sxi) * rdxidx))
  dgdy    =    dgdvar0 * ((w * cxi * drdy) + (w * (rw - sxi) * rdxidy))
  dgdw    =    dgdvar0 * r * cxi
  dgdpsi  =    dgdvar0 * rw * (rw - sxi)

  --  fill matrix:
  -- d w / d r, d phi, d z
  [ a11, a12, a13 ]  =  [ 0, 0, 0 ]
  -- d tl / d x, d y, d z
  [ a21, a22, a23 ]  =  [ 0, 0, 0 ]
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
  [ b41, b42, b43 ]  =  [ - oow*oow*(1.0 - 1.0/cg)
                            - (oow - r*sxi)*sg*dgdw/cg/cg,
                          0,
                          r*cxi/cg - (oow - r*sxi)*sg*dgdpsi/cg/cg ]
  -- d z0 / d w, d tl, d psi
  [ b51, b52, b53 ]  =  [ tl/w*(gamma/w - dgdw),
                          -gamma/w,
                          -tl/w*dgdpsi
                        ]
  [ v01, v02, v03 ] = toList 3 v
  [ q01, q02, q03 ] = toList 3 h
  h0 = fromList 5 [
      0,
      0,
      psi0 - a31*v01 - a32*v02 - b31*q01 - b33*q03,
      d0 - a41*v01 - a42*v02 - b41*q01 - b43*q03,
      z0 - a51*v01 - a52*v02 - a53*v03 - b51*q01 - b52*q02 - b53*q03]
  aa = fromList2 5 3 [a11,a12,a13,a21,a22,a23,a31,a32,a33,a41,a42,a43,a51,a52,a53]
  bb = fromList2 5 3 [b11,b12,b13,b21,b22,b23,b31,b32,b33,b41,b42,b43,b51,b52,b53]

