-- file src/Fit.hs
module Fit ( M, XVec (..), HVec (..), QVec (..), PVec (..), cov, hel ) where

import Data.Matrix ( Matrix )

type M = Matrix Double

data XVec = XVec (M, M) deriving Show -- 3-vector and covariance matrix for position/vertex measurement

data HVec = HVec (M, M) deriving Show -- 5-vector and covariance matrix for helix measurement
hel :: HVec -> M
hel hv = h where HVec (h, _) = hv
cov :: HVec -> M
cov hv = c where HVec (_, c) = hv

data QVec = QVec (M, M) deriving Show -- 3-vector and covariance matrix for momentum measurement

data PVec = PVec (M, M) deriving Show -- 4-vector and coavariance matrix for momentum px,py,pz and energy
