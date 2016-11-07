module Matrix ( inv, tr, sw, sub, sub2, scalar, scaleDiag, scale
              , toList, fromList, fromList2 ) where

import Debug.Trace ( trace )
import Text.Printf
import qualified Data.Matrix ( Matrix, inverse
                             , identity, nrows, transpose, elementwise
                             , rowVector, colVector, getCol, multStd2
                             , zero, scaleMatrix
                             , diagonal, getDiag
                   , getMatrixAsVector, submatrix, toList, fromList, (!) )
import qualified Data.Vector (Vector)

type M     = Data.Matrix.Matrix Double
type V     = Data.Matrix.Matrix Double

debug = flip trace

-- vectors are column-wise, represented as matrix of dimension nx1
sub :: Int -> M -> M
sub n v = Data.Matrix.submatrix 1 n 1 1 v
sub2 :: Int -> M -> M
sub2 n m = Data.Matrix.submatrix 1 n 1 n m
scalar :: M -> Double
scalar m = m Data.Matrix.! (1,1)
toList :: Int -> M -> [Double]
toList n m = take n $ Data.Matrix.toList m

fromList :: Int -> [Double] -> M
fromList rows ds = Data.Matrix.fromList rows 1 ds -- column vector to list
fromList2 :: Int -> Int -> [Double] -> M
fromList2 rows cols ds = Data.Matrix.fromList rows cols ds

scale :: Double -> M -> M
scale s = Data.Matrix.scaleMatrix s

scaleDiag :: Double -> M -> M
scaleDiag s = (Data.Matrix.diagonal 0.0 . Data.Matrix.getDiag . Data.Matrix.scaleMatrix  s)

tr :: M -> M
tr m = Data.Matrix.transpose m

sw :: M -> M -> M
sw a b = (Data.Matrix.transpose a) * b* a

inv' :: M -> M
inv' m = either invErr id (Data.Matrix.inverse m)
  where invErr s = (Data.Matrix.identity $ Data.Matrix.nrows m) `debug` ("🚩" ++ s)

inv :: M -> M
inv m = f e where
  e = Data.Matrix.inverse m
  f :: Either String M -> M
  f (Left s) = (Data.Matrix.identity $ Data.Matrix.nrows m) `debug` ("🚩" ++ s) -- can't invert
  f (Right m') = fdeb mx where
    mxx = Data.Matrix.elementwise (/)
                      (Data.Matrix.elementwise (-) (Data.Matrix.multStd2 m m') (Data.Matrix.identity (Data.Matrix.nrows m)))
                      m
    mx = (* 1000.0) . maximum  . Data.Matrix.getMatrixAsVector $ mxx
    fdeb mx
      | mx < 1.0 = m'
      | otherwise = let
                        sx :: String; sx = printf "%8.3f" (mx :: Double)
                    in m' --`debug` ("^" ++ "inv: max " ++ sx ++ " permille" )

