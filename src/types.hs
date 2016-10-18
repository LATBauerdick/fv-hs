-- file src/Fit.hs
module Types (
  M, V, M33, V3, V5, C44, XMeas (..), HMeas (..), QMeas (..)
             , PMeas (..), MMeas (..), Prong (..), VHMeas (..)
             , ABh0 (..), Chi2
             , showXMeas, showPMeas, showQMeas, showMMeas
             , w2pt, mπ
             ) where

import Text.Printf
import qualified Data.Matrix ( Matrix, getDiag, getCol, toList, zero
                             , fromLists, transpose )
import qualified Data.Vector ( zip, fromList, toList )

w2pt :: Double
w2pt = 4.5451703E-03

mπ :: Double
mπ = 0.1395675E0


type M     = Data.Matrix.Matrix Double
type V     = Data.Matrix.Matrix Double
type M33   = Data.Matrix.Matrix Double
type V3    = Data.Matrix.Matrix Double
type V5    = Data.Matrix.Matrix Double
type N     = Int
type Chi2  = Double
data Prong = Prong N XMeas [QMeas] [Chi2] deriving Show -- a prong results from a vertex fit of N helices

data VHMeas a = VHMeas XMeas [a] deriving Show
instance Monoid (VHMeas a) where
  mappend (VHMeas v1 as1) (VHMeas v2 as2) = VHMeas (v1) ( as1 ++ as2 )
instance Functor VHMeas where -- apply f to each a
  fmap f (VHMeas v (a:as)) = VHMeas v ((f a):(fmap f as))

data ABh0 = ABh0 M M M

type X3 = V
type C33 = M
data XMeas = XMeas X3 C33 deriving Show -- 3-vector and covariance matrix for position/vertex measurement

type H5 = V
type C55 = M
data HMeas = HMeas H5 C55 deriving Show -- 5-vector and covariance matrix for helix measurement

type Q3 = V
data QMeas = QMeas Q3 C33 deriving Show -- 3-vector and covariance matrix for momentum measurement

type P4 = V -- four-vector
type C44 = M -- 4x4 covariance matrix
data PMeas = PMeas P4 C44 deriving Show -- 4-vector and coavariance matrix for momentum px,py,pz and energy

type D = Double
data MMeas = MMeas D D deriving Show -- mass and error

-- show instances -- refactor!!

-- print a vertext position vector with errors
showXMeas :: String -> XMeas -> IO ()
showXMeas s (XMeas v cv) = do
  putStr s
  let
    s2v        = Data.Matrix.getDiag cv
    f (x, s2)  = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
      where dx = sqrt s2
    in mapM_ f $ Data.Vector.zip (Data.Matrix.getCol 1 v) s2v
  putStrLn " cm"

showMMeas :: String -> MMeas -> IO ()
showMMeas s (MMeas m dm) = do
  putStr s
  printf "%8.3f ± %8.3f" (m::Double) (dm::Double)
  putStrLn " GeV"

-- print a 4-momentum vector with errors
showPMeas :: String -> PMeas -> IO ()
showPMeas s (PMeas p cp) = do
  putStr s
  let
    s2p        = Data.Matrix.getDiag cp
    f (x, s2)  = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
      where dx = sqrt s2
    in mapM_ f $ Data.Vector.zip (Data.Matrix.getCol 1 p) s2p
  putStrLn " GeV"


-- print a 4-momentum vector with errors
showQMeas :: String -> QMeas -> IO ()
showQMeas s (QMeas q cq) = do
  putStr s
  let
    m = mπ
    wp = w2pt
    [w,tl,psi0] = take 3 (Data.Matrix.toList q)
    pt   = wp / abs w
    pz = pt*tl
    psi = psi0*180.0/pi
    e = sqrt(pt^2  + pz^2 + m^2)

    dpdk = pt*pt/wp
    [c11, c12, c13, _, c22, c23, _, _, c33] = take 9 $ Data.Matrix.toList cq
    dpt = sqrt $ (dpdk*dpdk) * c11 + pt*pt * c33
    dpz = sqrt $ (dpdk*tl)^2 * c11 + pt*pt*c22 - 2.0*wp/w*dpdk*tl*c12
    dpsi = (sqrt c22)*(180.0/pi)
    sptz = 0 -- fix!!!!!!!!!!!!!
    de = sqrt(pt*pt*dpt*dpt + pz*pz*dpz*dpz + 2.0*pt*pz*sptz)/e

    p'           = Data.Vector.fromList [pt, pz, psi, e]
    dp'         = Data.Vector.fromList [dpt,dpz,dpsi,de]
    f (x, dx)    = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
    in mapM_ f $ Data.Vector.zip p' dp'
  putStrLn " GeV"
    -- cc = cos phi
    -- ss = sin phi
    -- xy = 1.0/(1+py*py/px/px)
    -- j_xyz = Data.Matrix.fromLists [
    --                     [ cc, 0, -pt*ss, 0 ],
    --                     [ ss, 0, pt*cc, 0 ],
    --                     [ 0, 1.0, 0, 0 ],
    --                     [ 0, 0, 0, 1.0 ]
    --                   ]
    -- j_ptphiz = Data.Matrix.fromLists [
    --                       [ px/pt, py/pt, 0, 0 ],
    --                       [ 0, 0, 1.0, 0 ],
    --                       [ -py/px/px*xy, 1.0/px*xy, 0, 0 ],
    --                       [ 0, 0, 0, 1.0 ]
    --                      ]
    -- cp' = (Data.Matrix.transpose j_ptphiz) * cp * j_ptphiz
    -- s2p          = Data.Matrix.getDiag cp'
    -- [s0,s1,s2,s3]= Data.Vector.toList s2p

