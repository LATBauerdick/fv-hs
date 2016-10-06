-- file src/coeff.hs
--
module Coeff ( w2pt, h2p4 ) where

import Data.Matrix ( identity, fromList, fromLists, (!) )
import Fit -- ( M, HVec, PVec )

w2pt :: Double
w2pt = 4.5451703E-03

mπ :: Double
mπ = 0.1395675E0

h2p4 :: HVec -> PVec
h2p4 hh = p where
  m = mπ
  h = hel hh
  ch = cov hh
  w    = h!(1,1)
  tl   = h!(2,1)
  psi0 = h!(3,1)
  sph  = sin psi0
  cph  = cos psi0
  pt   = w2pt / abs w
  px   = pt * cph
  py   = pt * sph
  pz   = pt * tl
  e = sqrt(px^2 + py^2 + pz^2 + m^2)
  ps = w2pt / w
  dpdk = ps*ps/w2pt
  xy = 2.0*ps*dpdk*cph*sph*ch!(1,3)
  sxx = (dpdk*cph)^2 * ch!(1,1) + (ps*sph)^2 * ch!(3,3) + xy
  sxy = cph*sph*(dpdk*dpdk*ch!(1,1) - ps*ps*ch!(3,3)) +
           ps*dpdk*(sph*sph-cph*cph)*ch!(1,3)
  syy = (dpdk*sph)^2 * ch!(1,1) + (ps*cph)^2 * ch!(3,3) - xy
  sxz = dpdk*dpdk*cph*tl*ch!(1,1) -
           ps*dpdk*(cph*ch!(1,2)-sph*tl*ch!(1,3)) -
           ps*ps*sph*ch!(2,3)
  syz = dpdk*dpdk*sph*tl*ch!(1,1) -
           ps*dpdk*(sph*ch!(1,2) + cph*tl*ch!(1,3)) +
           ps*ps*cph*ch!(2,3)
  szz = (dpdk*tl)^2 * ch!(1,1) + ps*ps*ch!(2,2) -
           2.0*ps*dpdk*tl*ch!(1,2)
  sxe = (px*ch!(1,1) + py*ch!(1,2) + pz*ch!(1,3))/e
  sye = (px*ch!(1,2) + py*ch!(2,2) + pz*ch!(2,3))/e
  sze = (px*ch!(1,3) + py*ch!(2,3) + pz*ch!(3,3))/e
  see = (px*px*ch!(1,1) + py*py*ch!(2,2) + pz*pz*ch!(3,3) +
           2.0*(px*(py*ch!(1,2) + pz*ch!(1,3)) + py*pz*ch!(2,3)))/e/e


  p0 = fromList 1 4 [px,py,pz,e ] :: M
  cp0 = fromLists [[sxx, sxy, sxz, sxe], [sxy, syy, syz, sye], [sxz, syz, szz, sze], [sxe, sye, sze, see]] :: M
  p = PVec (p0 , cp0)

