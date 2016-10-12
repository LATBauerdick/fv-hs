module Matrix ( inv, vATBA, mATBA, mABAT, mATB, mAMB, mAPB
              , vAMB, vAPB, mAv ) where

import Debug.Trace ( trace )
import Data.Matrix ( inverse, identity, nrows, transpose, elementwise
                   , rowVector, colVector, getCol, multStd2
                   , getMatrixAsVector )
import Types ( M, V )

debug = flip trace

vAPB :: V -> V -> V
vAPB a b = b -- zipWith (\ a b -> a+b) a b

vAMB :: V -> V -> V
vAMB a b = b -- zipWith (-) a b

mAv :: M -> V -> V
mAv a v = getCol 1 $ multStd2 a (colVector v)

mAPB :: M -> M -> M
mAPB a b =  elementwise (+) a b

mAMB :: M -> M -> M
mAMB a b =  elementwise (-) a b

vATBA :: V -> M -> M
vATBA a b =
  multStd2 (rowVector a) (multStd2 b (colVector a))

mATB :: M -> M -> M
mATB a b = multStd2 (transpose a) b

mABAT :: M -> M -> M
mABAT a b = multStd2 a (multStd2 b (transpose a))

mATBA :: M -> M -> M
mATBA a b = multStd2 (transpose a) (multStd2 b a)

inv' :: M -> M
inv' m = either invErr id (Data.Matrix.inverse m)
  where invErr s = (identity $ nrows m) `debug` ("🚩" ++ s)

inv :: M -> M
inv m = f e where
  e = Data.Matrix.inverse m
  f :: Either String M -> M
  f (Right m') = m' `debug` ("^" ++ show mx ) where
    mx = (* 100.0) . maximum  . getMatrixAsVector $ mxx
    mxx = elementwise (/)
                      (elementwise (-) (multStd2 m m') (identity (nrows m)))
                      m
  f (Left s) = (identity $ nrows m) `debug` ("🚩" ++ s)

