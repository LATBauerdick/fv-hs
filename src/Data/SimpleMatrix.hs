{-# LANGUAGE PostfixOperators #-}

module Data.SimpleMatrix ( inv, invMaybe, det, tr, sw, chol
              , sub, sub2, scalar, scaleDiag, diagonal, scale
              , toList, fromList, fromList2, zero ) where

import Prelude
import Debug.Trace ( trace )
import Text.Printf
import qualified Data.Matrix as M (
                                Matrix, inverse, cholDecomp, detLU
                              , identity, nrows, transpose, elementwise
                              , rowVector, colVector, getCol, multStd2
                              , zero, scaleMatrix
                              , diagonal, getDiag, diagonalList
                              , getMatrixAsVector, submatrix
  , toList, fromList, (!) , (<|>), (<->), combineRows, getElem, ncols
  , splitBlocks, switchRows, scaleRow
                                  )
import qualified Data.Vector (Vector)
import Data.Maybe ( listToMaybe )
import Control.Monad.Error

type M     = M.Matrix Double
type V     = M.Matrix Double

debug = flip trace

-- vectors are column-wise, represented as matrix of dimension nx1
sub :: Int -> M -> M
sub n v = M.submatrix 1 n 1 1 v
sub2 :: Int -> M -> M
sub2 n m = M.submatrix 1 n 1 n m
scalar :: M -> Double
scalar m = m M.! (1,1)

toList :: Int -> M -> [Double]
toList n m = take n $ M.toList m

fromList :: Int -> [Double] -> M
fromList rows ds = M.fromList rows 1 ds -- column vector to list
fromList2 :: Int -> Int -> [Double] -> M
fromList2 rows cols ds = M.fromList rows cols ds

zero :: Int -> Int -> M
zero rows cols = M.zero rows cols
scale :: Double -> M -> M
scale s = M.scaleMatrix s

diagonal :: [Double] -> M
diagonal d = M.diagonalList (length d) 0.0 d

scaleDiag :: Double -> M -> M
scaleDiag s = (M.diagonal 0.0 . M.getDiag . M.scaleMatrix  s)

tr :: M -> M
tr = M.transpose

(^+) :: M -> M
(^+) = M.transpose

sw :: M -> M -> M
sw a b = (tr a) * b * a

chol :: M -> M
chol a = M.cholDecomp a

det :: M -> Double
det m = M.detLU m

-- This is the type of our Inv error representation.
data InvError = Err { quality::Double, reason::String }

-- We make it an instance of the Error class
instance Error InvError where
  noMsg    = Err 0 "Inversion Error"
  strMsg s = Err 0 s

invMaybe :: M -> Maybe M
invMaybe m = case invsm m of
               Right im -> Just im
               Left s -> Nothing `debug` ("Error in Matrix.invsm: " ++ s)

inv :: M -> M
inv m =  let (Right m') = do { invm m } `catchError` printError
          in m'
         where
           one = (M.identity $ M.nrows m)
           printError :: InvError -> InvMonad M
           printError e = return one `debug` ("Error in Matrix.inv: " ++ (show (quality e)) ++ ": " ++ (reason e))

-- For our monad type constructor, we use Either InvError
-- which represents failure using Left InvError or a
-- successful result of type a using Right a.
type InvMonad = Either InvError

invm :: M -> InvMonad M
invm m = case invsm m of
            Right m'  -> return m'
            Left s    -> throwError (Err 0 ("In Matrix.invm: " ++ s)) -- `debug` "yyyyyyyy"

-- inverse of a square matrix, from Data.Matrix with fix
--   Uses naive Gaussian elimination formula.
invsm ::  M -> Either String M
invsm m = rref'd >>= return . M.submatrix 1 n (n + 1) (n * 2) where
            n = M.nrows m
            adjoinedWId = m M.<|> M.identity n
            rref'd = rref adjoinedWId

rref :: M -> Either String M
rref m = rm where
    rm = case ref m of
           Right r -> rrefRefd r
           Left s -> Left s
    rrefRefd mtx
      | M.nrows mtx == 1    = Right mtx
      | otherwise =
            let
                resolvedRight = foldr (.) id (map resolveRow [1..col-1]) mtx
                    where
                    col = M.nrows mtx
                    resolveRow n = M.combineRows n (-M.getElem n col mtx) col
                top = M.submatrix 1 (M.nrows resolvedRight - 1) 1 (M.ncols resolvedRight) resolvedRight
                top' = rrefRefd top
                bot = M.submatrix (M.nrows resolvedRight) (M.nrows resolvedRight) 1 (M.ncols resolvedRight) resolvedRight
            in top' >>= return . (M.<-> bot)

ref :: M -> Either String M
ref mtx
        | M.nrows mtx == 1
            = Right clearedLeft
        | goodRow == 0
            = Left ("In Matrix.ref: Attempt to invert a non-invertible matrix") -- `debug` "xxxxxxxx"
        | otherwise =
            let
                (tl, tr, bl, br) = M.splitBlocks 1 1 clearedLeft
                br' = ref br
            in case br' of 
                  Right br'' -> Right ((tl M.<|> tr) M.<-> (bl M.<|> br''))
                  Left s -> Left s
    where
      goodRow = case listToMaybe (filter (\i -> M.getElem i 1 mtx /= 0) [1..M.nrows mtx]) of -- ERROR in orig: ncols
                  Nothing   -> 0
                  Just x -> x
      sigAtTop = M.switchRows 1 goodRow mtx
      normalizedFirstRow = M.scaleRow (1 / M.getElem 1 1 mtx) 1 sigAtTop
      clearedLeft = foldr (.) id (map combinator [2..M.nrows mtx]) normalizedFirstRow where
        combinator n = M.combineRows n (-M.getElem n 1 normalizedFirstRow) 1


inv' :: M -> M
inv' m = either invErr id (M.inverse m)
  where invErr s = (M.identity $ M.nrows m) `debug` ("🚩" ++ s)

inv''' :: M -> M
inv''' m = f e where
  e = M.inverse m
  f :: Either String M -> M
  f (Left s) = (M.identity $ M.nrows m) `debug` ("🚩" ++ s) -- can't invert
  f (Right m') = fdeb mx where
    mxx = M.elementwise (/)
                      (M.elementwise (-) (M.multStd2 m m') (M.identity (M.nrows m)))
                      m
    mx = (* 1000.0) . maximum  . M.getMatrixAsVector $ mxx
    fdeb mx
      | mx < 1.0 = m'
      | otherwise = let
                        sx :: String; sx = printf "%8.3f" (mx :: Double)
                    in m' `debug` ("^" ++ "inv: max " ++ sx ++ " permille" )
