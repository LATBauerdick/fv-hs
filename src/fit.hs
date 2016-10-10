-- file src/Fit.hs
module Fit ( fit ) where

import Data.Matrix ( Matrix, submatrix )
import Data.Vector ( take )
import Types ( XVec (..), HVec (..), QVec (..), Prong (..) )
import Debug.Trace ( trace )
debug = flip trace


qvec :: HVec -> QVec
qvec (HVec h ch) = q where
  q = QVec (Data.Vector.take 3 h) (submatrix 1 3 1 3 ch)

fit :: XVec -> [HVec] -> Prong
fit v0 hl = pr where
  v = filt v0 hl
  pr = smooth v hl

filt :: XVec -> [HVec] -> XVec
filt v0 hl = v where
  ƒ :: HVec -> XVec -> XVec
  ƒ h (XVec v cv) = (XVec v cv) `debug` (".")
  v = foldr ƒ v0 hl



-- kalman smooth: calculate helices hl at kalman filter vertex v
smooth :: XVec -> [HVec] -> Prong
smooth v hl =
  Prong 6 v (map qvec hl) [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
