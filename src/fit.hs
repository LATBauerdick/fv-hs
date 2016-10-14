-- file src/Fit.hs
module Fit ( fit ) where

import Data.Matrix ( Matrix, submatrix )
import Data.Vector ( take )
import Types ( XMeas (..), HMeas (..), QMeas (..), Prong (..), ABh0 (..) )
import Coeff ( fvABh0 )
import Matrix ( inv, tr, sw, cv, tov, mATBA, mABAT, mATB, mAMB, mAPB, vAMB, vAPB, mAv )
import Debug.Trace ( trace )
debug = flip trace


qMeas :: HMeas -> QMeas
qMeas (HMeas h ch) = q where
  q = QMeas (Data.Vector.take 3 h) (submatrix 1 3 1 3 ch)

fit :: XMeas -> [HMeas] -> Prong
fit v0 hl = pr where
  v = kfilter v0 hl
  pr = smooth v hl

kfilter :: XMeas -> [HMeas] -> XMeas
kfilter v0 hl = v where
  v = foldr kal v0 hl

kal :: HMeas -> XMeas -> XMeas
kal (HMeas h_ hh) (XMeas v0_ vv0) = vm' `debug` (".") where
  h = cv h_
  v0 = cv v0_
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
  vv'' = inv . inv $ vv0
  v_ = tov v
  vm' = XMeas v_ vv''

    -- (loop [v0 v0, U0 U0, A A, B B, h0 h0, 𝜒20 1e10, iter 0]
    --   (let [
    --         ;; -- W = (B^T.G.B)^(-1)
    --         ;; -- GB = G - G.B.W.B^T.G^T
    --         ;; -- C = U^-1 = (U0 + A^T.GB.A)^-1
    --         U            (fvAPB U0 (fvsATBA A GB))
    --         C            (fvInverse U)   ;; check for editing of singular values?
    --         ;; -- m = h - h0
    --         m      (fvAMB h h0)
    --         ;; -- v = C. (U0.v0 + A^T.GB.(h-h0) )
    --         v      (fvAB C (fvAPB (fvAB U0 v0) (fvATB A (fvAB GB m))))
    --         ;; -- dm = h - h0 - A.v
    --         dm     (fvAMB m (fvAB A v))
    --         ;; -- q = W.BT.(G.dm)
    --         q      (fvAB  W  (fvATB  B  (fvAB  G  dm)))
    --         ;; -- D
    --         D      (fvAPB W (fvsATBA W (fvsATBA B (fvsATBA G (fvsABAT A C)))))
    --         ;; -- E
    --         E     (fvNegA  (fvATBC  W  (fvATBC  B  G  A)  C))
    --         ;; -- chi2 - (dm - B.q)T.G. (dm - B.q) + (v - v0)T.U0. (v-v0)
    --         𝜒2    (+ (scalar  (fvsATBA  (fvAMB  dm  (fvAB  B  q))  G))
    --                  (scalar  (fvsATBA  (fvAMB  v  v0)  U0)))
    --         ]
    --     (when-not fvLog (print ".")) ;; progress bar
    --     (when fvLog
    --       (print "Filter iteration " iter "yields chi2 " (format "%9.3g" 𝜒2))
    --       (fvPMerr " at x,y,z," v C))
    --     (if-not (goodEnough? 𝜒20 𝜒2)
    --       (do
    --         (let [ [A1 B1 h01]     (fvABh0 v (fvq h v))] ;; recalc derivs at v
    --           (recur v0 U0 A1 B1 h01 𝜒2 (inc iter))))
    --       [v U q (fvQ H) 𝜒2])))))



-- kalman smooth: calculate helices hl at kalman filter vertex v
smooth :: XMeas -> [HMeas] -> Prong
smooth v hl =
  Prong 6 v (map qMeas hl) [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
