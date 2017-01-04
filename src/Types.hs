-- file src/Types.hs
--
{-# LANGUAGE RankNTypes #-}

module Types (
                M, V, M33, V3, V5, C44
             , XMeas (..), HMeas (..), QMeas (..)
             , PMeas (..), MMeas (..), DMeas (..), Prong (..), VHMeas (..)
  , helicesLens
  , view, over, set
  , vBlowup, hFilter, hRemove
             , Mom (..), Pos (..)
             , X3, C33, Q3, H5, C55
             , Jaco (..), Chi2
             , v3, l3, v5, l5
             , h2p, h2q, q2p
             , mπ, invMass
             ) where

import Text.Printf
import Data.Foldable
import qualified Data.Matrix ( Matrix, getDiag, getCol, toList
                             , fromLists, transpose )
import qualified Data.Vector ( zip, map, fromList, toList, drop )
import qualified Matrix ( scalar, sw, tr, sub, sub2
                        , toList, fromList, fromList2
                        , zero, scaleDiag
                        )

mπ :: Double
mπ = 0.1395675E0


type M     = Data.Matrix.Matrix Double
type V     = Data.Matrix.Matrix Double
type M33   = Data.Matrix.Matrix Double
type V3    = Data.Matrix.Matrix Double
type V5    = Data.Matrix.Matrix Double
type Chi2  = Double
data Prong = Prong { -- a prong results from a vertex fit of N helices
    nProng        :: Int
  , fitVertex     :: XMeas
  , fitMomenta    :: [QMeas] 
  , fitChi2s      :: [Chi2]
  , measurements  :: VHMeas
                   } deriving Show

data VHMeas = VHMeas {
    vertex      :: XMeas
  , helices     :: [HMeas]
                     } deriving Show
instance Monoid VHMeas where
  mappend (VHMeas v hs) (VHMeas _ hs') = VHMeas v ( hs ++ hs' ) -- ???
  mempty = VHMeas (XMeas (Matrix.zero 3 1) (Matrix.zero 3 3)) []

data Jaco = Jaco M M M

type H5 = V
type C55 = M
data HMeas = HMeas H5 C55 Double -- 5-vector and covariance matrix for helix measurement
instance Show HMeas where
  show = showHMeas

type Q3 = V
data QMeas = QMeas Q3 C33 Double -- 3-vector and covariance matrix for momentum measurement
instance Show QMeas where
  show = showQMeas

-- 4-vector and coavariance matrix for momentum px,py,pz and energy
type P4 = V -- four-vector
type C44 = M -- 4x4 covariance matrix
data PMeas = PMeas P4 C44

instance Show PMeas where
  show = showPMeas

instance Monoid PMeas where
  mappend (PMeas p1 cp1) (PMeas p2 cp2) = PMeas (p1+p2) (cp1 + cp2)
  mempty = PMeas (Matrix.zero 4 1) (Matrix.zero 4 4)

-- instance Functor PMeas where
--   fmap f (PMeas p cp) = f p cp

invMass :: [PMeas] -> MMeas
invMass ps = mass . fold $ ps

class Mom a where
  mass :: a -> MMeas
  energy :: a -> Double

instance Mom PMeas where
  mass = pmass
  energy (PMeas p _) = e' where
    [_,_,_,e'] = Matrix.toList 4 p

instance Mom QMeas where
  mass qm = pmass $ q2p qm
  energy qm = e' where
    PMeas p _ = q2p qm
    [_,_,_,e'] = Matrix.toList 4 p

data MMeas = MMeas Double Double -- mass and error
instance Show MMeas where
  show (MMeas m dm) = printf "%8.1f ± %8.1f MeV" (m*1000.0) (dm*1000.0)

-----------------------initial Lens stuff--------------

newtype Identity a = Identity { runIdentity :: a }
instance Functor Identity where
    fmap f (Identity a) = Identity (f a)

newtype Const a b = Const { getConst :: a }
instance Functor (Const a) where
  fmap _ (Const a) = Const a

type Lens s a = forall f.Functor f => (a -> f a) -> s -> f s

over :: Lens s a -> (a -> a) -> s -> s
over ln f s = runIdentity $ ln (Identity . f) s

view :: Lens s a -> s -> a
view ln s = getConst $ ln Const s

set :: Lens s a -> a -> s -> s
set ln x = over ln (const x)

helicesLens :: Lens VHMeas [HMeas]
helicesLens f vm = fmap (\x -> vm { helices = x }) (f ( helices vm ))

vertexLens :: Lens VHMeas XMeas
vertexLens f vm = fmap (\x -> vm { vertex = x }) (f ( vertex vm ))

vBlowup :: Double -> VHMeas -> VHMeas
vBlowup scale vm = over vertexLens (blowup scale) vm where
  blowup :: Double -> XMeas -> XMeas -- blow up diag of cov matrix
  blowup s (XMeas v cv) = XMeas v cv' where
    cv' = Matrix.scaleDiag s $ cv

-- filter list of objects given list of indices in [a]
-- return list with only those b that have  indices that  are in rng [a]
iflt :: ( Eq a, Enum a, Num a ) => [a] -> [b] -> [b]
iflt rng hl =
  [h | (h, i) <- zip hl [0..], i `elem` rng ]

irem :: (Eq a, Enum a, Num a) => a -> [b] -> [b]
irem indx hl = [ h | (h,i) <- zip hl [0..], i /= indx ]

hFilter :: [Int] -> VHMeas -> VHMeas
hFilter is (VHMeas v hl) = VHMeas v (iflt is hl)

hRemove :: Int -> VHMeas -> VHMeas
hRemove indx (VHMeas v hl) = VHMeas v (irem indx hl)
-----------------------positions--------------

type X3 = V
type C33 = M
data XMeas = XMeas X3 C33 -- 3-vector and covariance matrix for position/vertex measurement
instance Show XMeas where
  show = showXMeas
instance Monoid XMeas where
  mappend (XMeas x1 cx1) (XMeas x2 cx2) = XMeas (x1 + x2) (cx1 + cx2)
  mempty = XMeas (Matrix.zero 3 1) (Matrix.zero 3 3)

class Pos a where
  distance :: a -> a -> DMeas -- distance between two poitions

instance Pos XMeas where
  distance x1 x2 = xmDist x1 x2

data DMeas = DMeas Double Double -- distance and error
instance Show DMeas where
  show (DMeas d sd) = printf ("%6.2f ± %6.2f cm")(d::Double) (sd::Double)

v3 :: [Double] -> V3
v3 = Matrix.fromList 3
l3 :: V3 -> [Double]
l3 = Matrix.toList 3
v5 :: [Double] -> V5
v5 = Matrix.fromList 5
l5 :: V5 -> [Double]
l5 = Matrix.toList 5

-- return a string showing vertext position vector with errors
showXMeas :: XMeas -> String
showXMeas (XMeas v cv) = s' where
  s2v        = Data.Vector.map sqrt $ Data.Matrix.getDiag cv
  [x, y, z]  = Matrix.toList 3 v
  [dx,dy,dz] = Data.Vector.toList s2v
  f :: Double -> Double -> String -> String
  f x dx s  = s ++ (printf "%6.2f ± %6.2f" (x::Double) (dx::Double))
  s' = (f z dz) . (f y dy) . (f x dx) $
    "r =" ++ (show $ distance (XMeas v cv) mempty) ++ ", "

-- calculate distance between two vertices
xmDist :: XMeas -> XMeas -> DMeas
xmDist (XMeas v0 vv0) (XMeas v1 vv1) = DMeas d sd where
  [x0, y0, z0] = Matrix.toList 3 v0
  [x1, y1, z1] = Matrix.toList 3 v1

  d    = sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)

  dd   = Data.Matrix.fromLists [[(x0-x1)/d, (y0-y1)/d, (z0-z1)/d]]
  tem0 = Matrix.sw (Matrix.tr dd) vv0
  tem1 = Matrix.sw (Matrix.tr dd) vv1
  sd   = sqrt (Matrix.scalar tem0 + Matrix.scalar tem1)

-- print PMeas as a 4-momentum vector px,py,pz,E with errors
showPMeas :: PMeas -> String
showPMeas (PMeas p cp) = s' where
  sp         = Data.Vector.map sqrt $ Data.Matrix.getDiag cp
  f s (x, dx)  = s ++ printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
  s' = (foldl f "" $ Data.Vector.zip (Data.Matrix.getCol 1 p) sp) ++ " GeV"


-- print HMeas as a 5-parameter helix with errors
showHMeas :: HMeas -> String
showHMeas (HMeas h ch _) = s' where
  sh = Data.Vector.map sqrt $ Data.Matrix.getDiag ch
  s00 = printf "%10.3g ± %10.3g" (x::Double) (dx::Double) where
    x  = head (Data.Matrix.toList h)
    dx = head (Data.Vector.toList sh)
  s' = foldl f s00 (Data.Vector.drop 1 $ Data.Vector.zip (Data.Matrix.getCol 1 h) sh) where
    f s (x, dx)  = s ++ printf "%8.3f ± %8.3f" (x::Double) (dx::Double)

-- print QMeas as a 4-momentum vector with errors, use pt and pz
showQMeas :: QMeas -> String
showQMeas (QMeas q cq w2pt) = s' where
  f :: String -> (Double, Double) -> String
  f s (x, dx) = s ++ printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
  m = mπ
  wp = w2pt
  [w,tl,psi0] = take 3 (Data.Matrix.toList q)
  pt   = wp / abs w
  pz = pt*tl
  psi = psi0*180.0/pi
  e = sqrt(pt*pt  + pz*pz + m*m)
  jj   = Data.Matrix.fromLists [
            [-wp/w/w, -wp/w/w*tl,0, -(pz*pz + pt*pt)/w/e ]
          , [0, wp/w, 0, pt*pt*tl/e]
          , [0, 0, 1.0, 0] ]
  cq'  = (Data.Matrix.transpose jj) * cq* jj
  p'   = Data.Vector.fromList [pt, pz, psi, e]
--    dp'  = Data.Vector.map sqrt $ Data.Matrix.getDiag cq'
  [d1,d2,d3,d4]  = Data.Vector.toList $ Data.Vector.map sqrt $ Data.Matrix.getDiag cq'
  dp' = Data.Vector.fromList [d1,d2,d3*180.0/pi,d4]
  s' = (foldl f "" $ Data.Vector.zip p' dp' ) ++ " GeV"


h2p :: HMeas -> PMeas
h2p hm = (q2p . h2q) hm

h2q :: HMeas -> QMeas -- just drop the d0, z0 part... fix!!!!
h2q (HMeas h ch w2pt) = QMeas q cq w2pt where
  q = Matrix.sub 3 h
  cq = Matrix.sub2 3 ch

pmass :: PMeas -> MMeas
pmass (PMeas p cp) = mm  where
  [px,py,pz,e] = Matrix.toList 4 p
  [c11, c12, c13, c14, _, c22, c23, c24, _, _, c33, c34, _, _, _, c44]
        = Matrix.toList 16 cp
  m     = sqrt $ max (e*e-px*px-py*py-pz*pz) 0
  sigm0 = px*c11*px + py*c22*py + pz*c33*pz + e*c44*e +
            2.0*(px*(c12*py + c13*pz - c14*e)
               + py*(c23*pz - c24*e)
               - pz*c34*e)
  sigm  =  sqrt ( max sigm0 0 ) / m
  mm    = MMeas m sigm

q2p :: QMeas -> PMeas
q2p (QMeas q0 cq0 w2pt) = PMeas p0 cp0 where
  m = mπ
  [w,tl,psi0] = Matrix.toList 3 q0
  sph  = sin psi0
  cph  = cos psi0
  pt   = w2pt / abs w
  px   = pt * cph
  py   = pt * sph
  pz   = pt * tl
  sqr = \x -> x*x
  e = sqrt(px*px + py*py + pz*pz + m*m)
  ps = w2pt / w
  dpdk = ps*ps/w2pt
  [c11, c12, c13, _, c22, c23, _, _, c33] = Matrix.toList 9 cq0
  xy = 2.0*ps*dpdk*cph*sph*c13
  sxx = sqr (dpdk*cph) * c11 + sqr (ps*sph) * c33 + xy
  sxy = cph*sph*(dpdk*dpdk*c11 - ps*ps*c33) +
           ps*dpdk*(sph*sph-cph*cph)*c13
  syy = sqr (dpdk*sph) * c11 + sqr (ps*cph) * c33 - xy
  sxz = dpdk*dpdk*cph*tl*c11 -
           ps*dpdk*(cph*c12-sph*tl*c13) -
           ps*ps*sph*c23
  syz = dpdk*dpdk*sph*tl*c11 -
           ps*dpdk*(sph*c12 + cph*tl*c13) +
           ps*ps*cph*c23
  szz = sqr (dpdk*tl) * c11 + ps*ps*c22 -
           2.0*ps*dpdk*tl*c12
  sxe = (px*sxx + py*sxy + pz*sxz)/e
  sye = (px*sxy + py*syy + pz*syz)/e
  sze = (px*sxz + py*syz + pz*szz)/e
  see = (px*px*sxx + py*py*syy + pz*pz*szz +
         2.0*(px*(py*sxy + pz*sxz) + py*pz*syz))/e/e

  p0 = Matrix.fromList 4 [px,py,pz,e]
  cp0 = Matrix.fromList2 4 4 [sxx, sxy, sxz, sxe, sxy, syy, syz, sye, sxz, syz, szz, sze, sxe, sye, sze, see]

