{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE ExplicitForAll #-}
--{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
--{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RankNTypes #-}
--{-# LANGUAGE RebindableSyntax #-}
--{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE OverloadedLists #-}
--{-# LANGUAGE NamedFieldPuns #-}

module Data.Cov
    where

import Prelude
import qualified Data.Vector.Unboxed as A
  ( Vector, length, fromList, toList, unsafeIndex, create, replicate
  , singleton, map, foldl, zipWith )
import qualified Data.Vector.Unboxed.Mutable as MA
  ( new, unsafeWrite, unsafeRead, unsafeTake )
import Control.Loop ( numLoop )
import Data.Foldable ( sum )
import Data.Maybe ( Maybe (..), fromJust )
import Control.Monad (guard)

-- import Data.Int ( toNumber, ceil )
-- import Math ( abs, sqrt )
-- import Unsafe.Coerce  as Unsafe.Coerce ( unsafeCoerce )

import qualified Data.SimpleMatrix as M
  ( Matrix
  , transpose
  , fromArray, fromArray2, toArray
  )

import Stuff

newtype Dim3 = DDim3 Int
newtype Dim4 = DDim4 Int
newtype Dim5 = DDim5 Int
class DDim a where
  ddim :: a -> Int
instance DDim Dim3 where
  ddim _ = 3
instance DDim Dim4 where
  ddim _ = 4
instance DDim Dim5 where
  ddim _ = 5
instance DDim a where
  ddim _ = undefined

-- this requires OverlappingTypeVariables
-- class Dim a where
--   dim :: a -> Int
-- instance Dim (Cov a) where
--   dim ccc = n where
--     xx = undefined::a
--     n = ddim xx
-- instance Dim (Cov Dim3) where
--   dim ccc = 3
-- instance Dim (Cov Dim4) where
--   dim ccc = 4
-- instance Dim (Cov Dim5) where
--   dim ccc = 5

newtype Cov a   = Cov { vc :: Array Number }
newtype Jac a b = Jac { vj :: Array Number }
newtype Vec a   = Vec { vv :: Array Number }
type Cov3 = Cov Dim3
type Cov4 = Cov Dim4
type Cov5 = Cov Dim5
{-- type Jac43 = Jac Dim4 Dim3 --}
type Jac53 = Jac Dim5 Dim3
type Jac33 = Jac Dim3 Dim3
type Jac34 = Jac Dim3 Dim4
type Jac35 = Jac Dim3 Dim5
type Jac44 = Jac Dim4 Dim4
type Jac55 = Jac Dim5 Dim5
type Vec3 = Vec Dim3
type Vec4 = Vec Dim4
type Vec5 = Vec Dim5
data Jacs = Jacs
            { aa :: Jac53
            , bb :: Jac53
            , h0 :: Vec5}

-- access to arrays of symmetrical matrices
uGet :: Array Number -> Int -> Int -> Int -> Number
uGet a w i j | i <= j     = uidx a ((i-1)*w - (i-1)*(i-2) `div` 2 + j-i)
             | otherwise = uidx a ((j-1)*w - (j-1)*(j-2) `div` 2 + i-j)
indV :: Int -> Int -> Int -> Int
indV w i0 j0 = (i0*w+j0)
indVs :: Int -> Int -> Int -> Int
indVs w i0 j0 | i0 <= j0  = (i0*w - i0*(i0-1) `div` 2 + j0-i0)
              | otherwise = (j0*w - j0*(j0-1) `div` 2 + i0-j0)

-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-- Mat to give behavior to Cov and Vec and Jac
-- ability to convert to and from Matrix and Array
-- while keeping info about dimensionality
-- also define Semiring and Ring functions
--

class Mat a where
  val :: a -> Array Number
  fromArray :: Array Number -> a
  toArray :: a -> Array Number
  elementwise :: (Number -> Number -> Number) -> a -> a -> a

instance Mat (Cov a) where
  val (Cov {vc=v}) = v
  fromArray a = c' where
    l = A.length a
    c' = case l of
      6   -> Cov {vc= a}
      10  -> Cov {vc= a}
      15  -> Cov {vc= a}
      _   -> Cov {vc= let
          n = floor . sqrt . fromIntegral $ l
          iv = indV n
        in A.fromList $ do -- only upper triangle
          i0 <- [0 .. (n-1)]
          j0 <- [i0 .. (n-1)]
          pure $ uidx a (iv i0 j0) }

  toArray c@(Cov {vc=v}) = v' where
    l = A.length v
    n = case l of
      6  -> 3
      10 -> 4
      15 -> 5
      _  -> error $ "matacov toArray not supported " <> show l
    iv = indVs n
    v' = A.fromList $ do
      i0 <- [0..(n-1)]
      j0 <- [0..(n-1)]
      pure $ uidx v (iv i0 j0)
  elementwise f (Cov {vc=va}) (Cov {vc=vb}) = (Cov {vc=vc}) where
    vc = A.zipWith f va vb
instance Mat (Vec a) where
  val (Vec {vv=v}) = v
  fromArray a = Vec {vv= a}
  toArray (Vec {vv=v}) = v
  elementwise f (Vec {vv=va}) (Vec {vv=vb}) = (Vec {vv=vc}) where
    vc = A.zipWith f va vb
instance Mat (Jac a b) where
  val (Jac {vj=v}) = v
  fromArray a = Jac {vj= a}
  toArray (Jac {vj=v}) = v
  elementwise f (Jac {vj=va}) (Jac {vj=vb}) = (Jac {vj=vc}) where
    vc = A.zipWith f va vb

class Mat1 a where
  toMatrix :: a -> M.Matrix
instance Mat1 (Cov a) where
  toMatrix (Cov {vc=v}) = case A.length v of
                            6  -> M.fromArray2 3 3 v
                            10 -> M.fromArray2 4 4 v
                            15 -> M.fromArray2 5 5 v
                            _ -> error $ "mat1Cova toMatrix "
                                          <> show (A.length v)
instance Mat1 (Vec a) where
  toMatrix (Vec {vv=v}) = M.fromArray (A.length v) v
instance Mat1 (Jac Dim5 Dim3) where
  toMatrix (Jac {vj=v}) = M.fromArray2 5 3 v `debug` "WTF??? 5 3"
instance Mat1 (Jac Dim3 Dim5) where
  toMatrix (Jac {vj=v}) = M.fromArray2 3 5 v `debug` "WTF??? 3 5"
instance Mat1 (Jac a b) where
  toMatrix j@(Jac {vj=v}) = case A.length v of
                              9  -> M.fromArray2 3 3 v
                              16 -> M.fromArray2 4 4 v
                              25 -> M.fromArray2 5 5 v
                              12 -> M.fromArray2 3 4 v `debug` "this should not have happened ??????????????????? 4 3"
                              15 -> M.fromArray2 5 3 v `debug` "this should not have happened ??????????????????? 5 3"
                              _  -> error $ "mat1Jacaa toMatrix "
                                          <> show (A.length v)

--{{{
--}}}
-----------------------------------------------------------------
-- | funcitons for symetric matrices: Cov
-- | type class SymMat
class SymMat a where
  inv :: a -> a                -- | inverse matrix
  invMaybe :: a -> Maybe a     -- | Maybe inverse matrix
  det :: a -> Number           -- | determinant
  diag :: a -> Array Number    -- | Array of diagonal elements
-- chol :: Cov a -> Jac a a                -- | Cholsky decomposition

instance SymMat (Cov Dim3) where
  inv m = uJust (invMaybe m)
  invMaybe (Cov {vc=v}) = _inv v where
    _inv :: Array Number -> Maybe (Cov Dim3)
    _inv [a11,a12,a13,a22,a23,a33] = do
      let det = (a33*a12*a12 - 2.0*a13*a23*a12 + a13*a13*a22
                +a11*(a23*a23 - a22*a33))
      guard $ (abs det) > 1.0e-50
      let
          b11 = (a23*a23 - a22*a33)/det
          b12 = (a12*a33 - a13*a23)/det
          b13 = (a13*a22 - a12*a23)/det
          b22 = (a13*a13 - a11*a33)/det
          b23 = (a11*a23 - a12*a13)/det
          b33 = (a12*a12 - a11*a22)/det
          v' = [b11,b12,b13,b22,b23,b33]
      pure $ Cov {vc=v'}
  det (Cov {vc=v}) = _det $ A.toList v where
    _det :: [Number] -> Number
    _det [a,b,c,d,e,f] = a*d*f - a*e*e - b*b*f + 2.0*b*c*e - c*c*d
  diag (Cov {vc=v}) = _diag v where
    _diag :: Array Number -> Array Number
    _diag [a11,_,_,a22,_,a33] = A.fromList [a11,a22,a33]
--  chol a = choldc a
instance SymMat (Cov Dim4) where
  inv m = uJust (invMaybe m)
  invMaybe (Cov {vc=v}) = _inv v where
    _inv :: Array Number -> Maybe (Cov Dim4)
    _inv [a,b,c,d,e,f,g,h,i,j] = do
      let det = (a*e*h*j - a*e*i*i - a*f*f*j + 2.0*a*f*g*i - a*g*g*h
            - b*b*h*j + b*b*i*i - 2.0*d*(b*f*i - b*g*h - c*e*i + c*f*g)
            + b*c*(2.0*f*j - 2.0*g*i) + c*c*(g*g - e*j) + d*d*(f*f - e*h))
      guard $ (abs det) > 1.0e-50
      let a' = (-j*f*f + 2.0*g*i*f - e*i*i - g*g*h + e*h*j)/det
          b' = (b*i*i - d*f*i - c*g*i + d*g*h + c*f*j - b*h*j)/det
          c' = (c*g*g - d*f*g - b*i*g + d*e*i - c*e*j + b*f*j)/det
          d' = (d*f*f - c*g*f - b*i*f - d*e*h + b*g*h + c*e*i)/det
          e' = (-j*c*c + 2.0*d*i*c - a*i*i - d*d*h + a*h*j)/det
          f' = (f*d*d - c*g*d - b*i*d + a*g*i + b*c*j - a*f*j)/det
          g' = (g*c*c - d*f*c - b*i*c + b*d*h - a*g*h + a*f*i)/det
          h' = (-j*b*b + 2.0*d*g*b - a*g*g - d*d*e + a*e*j)/det
          i' = (i*b*b - d*f*b - c*g*b + c*d*e + a*f*g - a*e*i)/det
          j' = (-h*b*b + 2.0*c*f*b - a*f*f - c*c*e + a*e*h)/det
      pure $ fromArray [a',b',c',d',e',f',g',h',i',j']
  det (Cov {vc=v}) = _det v where
    _det :: Array Number -> Number
    _det [a,b,c,d,e,f,g,h,i,j] =
        (a*e*h*j - a*e*i*i - a*f*f*j + 2.0*a*f*g*i - a*g*g*h
          - b*b*h*j + b*b*i*i - 2.0*d*(b*f*i - b*g*h - c*e*i + c*f*g)
          + b*c*(2.0*f*j - 2.0*g*i) + c*c*(g*g - e*j) + d*d*(f*f - e*h))
    _det _ = undefined
  diag (Cov {vc=v}) = _diag v where
    _diag :: Array Number -> Array Number
    _diag [a11,_,_,_,a22,_,_,a33,_,a44] = A.fromList [a11,a22,a33,a44]
--  chol a = choldc a
instance SymMat (Cov Dim5) where
  inv m = cholInv m
  invMaybe m = Just (cholInv m)
  det (Cov {vc=v}) = _det v where
    _det :: Array Number -> Number
    _det [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o] =
      a*f*j*m*o - a*f*j*n*n - a*f*k*k*o + 2.0*a*f*k*l*n - a*f*l*l*m
      - a*g*g*m*o + a*g*g*n*n + 2.0*a*g*h*k*o - 2.0*a*g*h*l*n - 2.0*a*g*i*k*n
      + 2.0*a*g*i*l*m - a*h*h*j*o + a*h*h*l*l + 2.0*a*h*i*j*n - 2.0*a*h*i*k*l
      - a*i*i*j*m + a*i*i*k*k - b*b*j*m*o + b*b*j*n*n + b*b*k*k*o
      - 2.0*b*b*k*l*n + b*b*l*l*m + 2.0*b*c*g*m*o - 2.0*b*c*g*n*n - 2.0*b*c*h*k*o
      + 2.0*b*c*h*l*n + 2.0*b*c*i*k*n - 2.0*b*c*i*l*m - 2.0*b*d*g*k*o
      + 2.0*b*d*g*l*n + 2.0*b*d*h*j*o - 2.0*b*d*h*l*l - 2.0*b*d*i*j*n
      + 2.0*b*d*i*k*l + 2.0*b*e*g*k*n - 2.0*b*e*g*l*m - 2.0*b*e*h*j*n
      + 2.0*b*e*h*k*l + 2.0*b*e*i*j*m - 2.0*b*e*i*k*k - c*c*f*m*o + c*c*f*n*n
      + c*c*h*h*o - 2.0*c*c*h*i*n + c*c*i*i*m + 2.0*c*d*f*k*o - 2.0*c*d*f*l*n
      - 2.0*c*d*g*h*o + 2.0*c*d*g*i*n + 2.0*c*d*h*i*l - 2.0*c*d*i*i*k
      - 2.0*c*e*f*k*n + 2.0*c*e*f*l*m + 2.0*c*e*g*h*n - 2.0*c*e*g*i*m
      - 2.0*c*e*h*h*l + 2.0*c*e*h*i*k - d*d*f*j*o + d*d*f*l*l + d*d*g*g*o
      - 2.0*d*d*g*i*l + d*d*i*i*j + 2.0*d*e*f*j*n - 2.0*d*e*f*k*l - 2.0*d*e*g*g*n
      + 2.0*d*e*g*h*l + 2.0*d*e*g*i*k - 2.0*d*e*h*i*j - e*e*f*j*m + e*e*f*k*k
      + e*e*g*g*m - 2.0*e*e*g*h*k + e*e*h*h*j
    _det _ = undefined
  diag (Cov {vc=v}) = _diag v where
    _diag :: Array Number -> Array Number
    _diag [a,_,_,_,_,b,_,_,_,c,_,_,d,_,e] = A.fromList [a,b,c,d,e]
--  chol a = choldc a

class MulMat a b c | a b -> c where
  (*.) :: a -> b -> c
infixr 7 *.
instance MulMat (Cov a) (Cov a) (Jac a a) where
  (*.) (Cov {vc= va}) (Cov {vc= vb}) = Jac {vj= vc} where
    na = case A.length va of
              6  -> 3
              10 -> 4
              15 -> 5
              _  -> error $ "mulMatCC wrong length of Cov v "
                            <> show (A.length va)
    vc = A.create $ do
      v <- MA.new $ na * na
      let ixa = indVs na
          ixb = indVs na
          ixc = indV na
      numLoop 0 (na-1) $ \i0 ->
        numLoop 0 (na-1) $ \j0 ->
          MA.unsafeWrite v (ixc i0 j0) $
          sum [ (uidx va (ixa i0 k0)) * (uidx vb (ixb k0 j0))
                 | k0 <- [0 .. na-1] ]
      pure v
instance MulMat (Jac a b) (Cov b) (Jac a b) where
  (*.) j@(Jac {vj= va}) c@(Cov {vc= vb}) = Jac {vj= vc} where
    nb = case A.length vb of
              6  -> 3
              10 -> 4
              15 -> 5
              _  -> error $ "mulMatJC wrong length of Cov v "
                            <> show (A.length vb)
    na = (A.length va) `div` nb
    vc :: Array Number
    vc = A.create $ do
      v <- MA.new $ na * nb
      let ixa = indV nb
          ixb = indVs na
          ixc = indV na
      numLoop 0 (na-1) $ \i0 ->
        numLoop 0 (nb-1) $ \j0 ->
          MA.unsafeWrite v (ixc i0 j0) $
          sum [ (uidx va (ixa i0 k0)) * (uidx vb (ixb k0 j0))
                 | k0 <- [0 .. nb-1] ]
      pure v
instance MulMat (Cov a) (Jac a b) (Jac a b) where
  (*.) c@(Cov {vc= va}) j@(Jac {vj= vb}) = Jac {vj= vc} where
    na = case A.length va of
              6  -> 3
              10 -> 4
              15 -> 5
              _  -> error $ "mulMatCJ wrong length of Cov v "
                            <> show (A.length va)
    nb = (A.length vb) `div` na
    vc = A.create $ do
      v <- MA.new $ na * nb
      let ixa = indVs na
          ixb = indV na
          ixc = indV na
      numLoop 0 (na-1) $ \i0 -> 
        numLoop 0 (nb-1) $ \j0 -> 
          MA.unsafeWrite v (ixc i0 j0) $
          sum [ (uidx va (ixa i0 k0)) * (uidx vb (ixb k0 j0)) 
                 | k0 <- [0 .. na-1] ]
      pure v
instance MulMat (Jac a b) (Vec b) (Vec a) where
  (*.) j@(Jac {vj= va}) v@(Vec {vv=vb}) = Vec {vv=vc} where
    nb = A.length vb
    na = (A.length va) `div` nb
    vc = A.create $ do
      v <- MA.new $ na
      let ixa = indVs nb
          ixb = indV 1
          ixc = indV 1

      numLoop 0 (na-1) $ \i0 -> 
        MA.unsafeWrite v i0 $
        sum [ (uidx va (ixa i0 k0)) * (uidx vb k0 )
                 | k0 <- [0 .. nb-1] ]
      pure v
instance MulMat (Jac a b) (Jac b a) (Jac a a) where -- Dim3 x Dim5
  (*.) (Jac {vj= va}) (Jac {vj= vb}) = Jac {vj= vc} where
    nb = case A.length va of
              12 -> 4
              15 -> 5
              9  -> 3
              16 -> 4
              25 -> 5
              _  -> error $ "mulMatJJ can only do 3x5 * 5x3, 3x4 * 4*3, or squares"
                            <> show (A.length vb)
    na = (A.length va) `div` nb
    vc :: Array Number
    vc = A.create $ do
      v <- MA.new $ na * nb
      let ixa = indV nb
          ixb = indV na
          ixc = indV na
      numLoop 0 (na-1) $ \i0 ->
        numLoop 0 (na-1) $ \j0 ->
          MA.unsafeWrite v (ixc i0 j0) $
          sum [ (uidx va (ixa i0 k0)) * (uidx vb (ixb k0 j0))
                 | k0 <- [0 .. nb-1] ]
      pure v
instance MulMat (Cov a) (Vec a) (Vec a) where
  (*.) (Cov {vc= va}) (Vec {vv=vb}) = Vec {vv=vc} where
    nb = A.length vb
    na = nb
    vc = A.create $ do
      v <- MA.new $ na
      let ixa = indVs na
          ixb = indV 1
          ixc = indV 1
      numLoop 0 (na-1) $ \i0 -> 
        MA.unsafeWrite v i0 $
        sum [ (uidx va (ixa i0 k0)) * (uidx vb k0 )
                 | k0 <- [0 .. na-1] ]
      pure v
instance MulMat (Vec a) (Vec a) Number where
  (*.) (Vec {vv=va}) (Vec {vv=vb}) = A.foldl (+) 0.0 $ A.zipWith (*) va vb
class TrMat a b | a -> b where
  tr :: a -> b
instance TrMat (Cov a) (Cov a) where
  tr c = c
instance TrMat (Jac a b) (Jac b a) where
  tr j@(Jac {vj=va}) = Jac {vj=vc} where
    l = A.length va
    na = case l of
              9 -> 3
              15 -> 5
              16 -> 4
              25 -> 5
              _  -> error $ "trMatJ: sorry, can't do anything but 5x3 and square "
                            <> show (A.length va)
    nb = l `div` na
    vc = A.create $ do
      v <- MA.new $ na*nb
      let ixa = indV nb
          ixc = indV na
      numLoop 0 (nb-1) $ \i0 ->
        numLoop 0 (na-1) $ \j0 ->
          MA.unsafeWrite v (ixc j0 i0) $ 12.0 -- uidx va (ixa i0 j0)
      pure v
class SW a b c | a b -> c where
  (.*.) :: a -> b -> c
--(.*.) = sw
infixl 7 .*.
instance SW (Vec a) (Cov a) Number where
  (.*.) v c = undefined
--  sw v c = n where
--    mv = toMatrix v
--    mc = toMatrix c
--    mc' = M.transpose mv * mc * mv
--    n = uidx (M.toArray mc') 0
instance SW (Cov a) (Cov a) (Cov a) where
  (.*.) ca cb = undefined
--  sw c1 c2 = c' where
--    j' = c1 *. c2 *. c1
--    c' = fromArray $ toArray j'
instance SW (Jac a b) (Cov a) (Cov b) where
  (.*.) (Jac {vj= va}) (Cov {vc= vb}) = Cov {vc= v'} where
    l = A.length vb
    n = case l of
              6  -> 3
              10 -> 4
              15 -> 5
              _  -> error $ "swJac: don'w know how to " <> show l
    m = (A.length va) `div` n -- > mxn * nxn * nxm -> mxm
    vint :: Array Number
    vint = A.create $ do
      v <- MA.new $ n*m
      let ixa = indVs n
      let ixb = indV m
      let ixc = indV m
      numLoop 0 (n-1) $ \i0 ->
        numLoop 0 (m-1) $ \j0 ->
          MA.unsafeWrite v (ixc i0 j0) $
            sum [ (uidx vb (ixa i0 k0)) * (uidx va (ixb k0 j0))
              | k0 <- [0 .. ( n-1)] ]
      pure v
    v' = A.create $ do
      v <- MA.new $ n*m
      let ixa = indV m
          ixb = indV m
          ixc = indVs m
      numLoop 0 (m-1) $ \i0 ->
        numLoop i0 (m-1) $ \j0 ->
          MA.unsafeWrite v (ixc i0 j0) $
            sum [ (uidx va (ixa k0 i0 )) * (uidx vint (ixb k0 j0))
              | k0 <- [0 .. (n-1)] ]
      pure v
-------------------------------------------------------
-------------------------------------------------------
---- NUMERICAL INSTANCE

instance Num (Cov a) where
  fromInteger i = Cov {vc= A.singleton <<< fromInteger $ i}
  negate (Cov {vc=v}) = Cov {vc=A.map negate v}
  abs (Cov {vc=v}) = Cov {vc=A.map abs v}
  signum (Cov {vc=v}) = Cov {vc=A.map signum v}
  (+) = elementwise (+)
  (*) = error "cannot multiply Cov*Cov to return a Cov, use *. instead"
instance Num (Vec a) where
  fromInteger i = Vec {vv= A.singleton <<< fromInteger $ i}
  negate (Vec {vv=v}) = Vec {vv=A.map negate v}
  abs (Vec {vv=v}) = Vec {vv=A.map abs v}
  signum (Vec {vv=v}) = Vec {vv=A.map signum v}
  (+) = elementwise (+)
  (*) = error "cannot multiply Vec*Vec to return a Vec, use *. instead"
instance Num (Jac a b) where
  fromInteger i = Jac {vj= A.singleton <<< fromInteger $ i}
  negate (Jac {vj=v}) = Jac {vj=A.map negate v}
  abs (Jac {vj=v}) = Jac {vj=A.map abs v}
  signum (Jac {vj=v}) = Jac {vj=A.map signum v}
  (+) = elementwise (+)
  (*) = error "cannot multiply Jac*Jac to return a Jac, use *. instead"

instance Show (Cov a) where
  show c = "Show (Cov a) \n" <> (show $ toMatrix c)
instance Show (Vec a) where
  show c = "Show (Vec a) \n" <> (show $ toMatrix c)
instance Show (Jac a b) where
  show c = "Show (Jac a b) \n" <> (show $ toMatrix c)

instance Semiring (Cov Dim3) where
  add (Cov {vc= va}) (Cov {vc= vb}) = Cov {vc= A.zipWith (+) va vb}
  zero = Cov {vc= A.replicate 6 0.0 }
  mul (Cov {vc= va}) (Cov {vc= vb}) = error "------------> mul cov3 * cov3 not allowed"
  one = Cov { vc= [1.0, 0.0, 0.0, 1.0, 0.0, 1.0] }
instance Ring (Cov Dim3) where
  sub (Cov {vc= va}) (Cov {vc= vb}) = Cov {vc= A.zipWith (-) va vb}

instance Semiring (Cov Dim4) where
  add (Cov {vc= va}) (Cov {vc= vb}) = Cov {vc= A.zipWith (+) va vb}
  zero = Cov {vc= A.replicate 10 0.0 }
  mul (Cov {vc= va}) (Cov {vc= vb}) = error "------------> mul cov4 * cov4 not allowed"
  one = Cov { vc= [1.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0] }
instance Ring (Cov Dim4) where
  sub (Cov {vc= va}) (Cov {vc= vb}) = Cov {vc= A.zipWith (-) va vb}

instance Semiring (Cov Dim5) where
  add (Cov {vc= va}) (Cov {vc= vb}) = Cov {vc= A.zipWith (+) va vb}
  zero = Cov {vc= A.replicate 15 0.0 }
  mul (Cov {vc= va}) (Cov {vc= vb}) = error "------------> mul cov5 * cov5 not allowed"
  one = Cov {vc= [1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0] }
instance Ring (Cov Dim5) where
  sub (Cov {vc= va}) (Cov {vc= vb}) = Cov {vc= A.zipWith (-) va vb}

instance Semiring (Jac a b) where
  add (Jac {vj= va}) (Jac {vj= vb}) = Jac {vj= A.zipWith (+) va vb}
  zero = undefined
  mul = undefined
  one = undefined
instance Ring (Jac a b) where
  sub (Jac {vj= va}) (Jac {vj= vb}) = Jac {vj= A.zipWith (-) va vb}

-- -- Instances for Vec -- these are always column vectors
-- instance Semiring (Vec Dim3) where
--   add (Vec {v= v1}) (Vec {v= v2}) = Vec {v= A.zipWith (+) v1 v2}
--   zero = Vec {v= A.replicate 3 0.0 }
--   mul (Vec {v= v1}) (Vec {v= v2}) = undefined
--   one = Vec { v= A.replicate 3 1.0 }

-- instance Semiring (Vec Dim4) where
--   add (Vec {v= v1}) (Vec {v= v2}) = Vec {v= A.zipWith (+) v1 v2}
--   zero = Vec {v= A.replicate 4 0.0 }
--   mul (Vec {v= v1}) (Vec {v= v2}) = undefined
--   one = Vec { v= A.replicate 4 1.0 }

-- instance Semiring (Vec Dim5) where
--   add (Vec {v= v1}) (Vec {v= v2}) = Vec {v= A.zipWith (+) v1 v2}
--   zero = Vec {v= A.replicate 5 0.0 }
--   mul (Vec {v= v1}) (Vec {v= v2}) = undefined
--   one = Vec { v= A.replicate 5 1.0 }

instance Semiring (Vec a) where
  add (Vec {vv= va}) (Vec {vv= vb}) = Vec {vv= A.zipWith (+) va vb}
  {-- zero = error "error calling zero for Vec a" -- Vec {v= A.replicate 5 0.0 } --}
  zero = Vec {vv= A.replicate 5 0.0 } `debug` "xxxxxxxxxxx>>> called Vec zero"
  mul = undefined
  {-- one = error "error calling one for Vec a" -- Vec { v= A.replicate 5 1.0 } --}
  one = Vec { vv= A.replicate 5 1.0 } `debug` "xxxxxxxxxxx>>> called Vec one"
instance Ring (Vec a) where
  sub (Vec {vv= va}) (Vec {vv= vb}) = Vec {vv= A.zipWith (-) va vb}

scaleDiag :: Number -> Cov3 -> Cov3
scaleDiag s (Cov {vc=v}) = (Cov {vc= _sc $ A.toList v}) where
  _sc :: [Number] -> Array Number
  _sc [a,_,_,b,_,c] = A.fromList [s*a,0.0,0.0,s*b,0.0,s*c]
  _sc _ = undefined

subm :: Int -> Vec5 -> Vec3
subm n (Vec {vv=v}) = Vec {vv= _subm v} where
  _subm :: Array Number -> Array Number
  _subm [a,b,c,_,_] = A.fromList [a,b,c]
  _subm _ = undefined

subm2 :: Int -> Cov5 -> Cov3
subm2 n (Cov {vc=v}) = Cov {vc= _subm2 v} where
  _subm2 :: Array Number -> Array Number
  _subm2 [a,b,c,_,d,e,_,_,f] = A.fromList [a,b,c,d,e,f]
  _subm2 _ = undefined


-- CHOLESKY DECOMPOSITION

-- | Simple Cholesky decomposition of a symmetric, positive definite matrix.
--   The result for a matrix /M/ is a lower triangular matrix /L/ such that:
--
--   * /M = LL^T/.
--
--   Example:
--
-- >            (  2 -1  0 )   (  1.41  0     0    )
-- >            ( -1  2 -1 )   ( -0.70  1.22  0    )
-- > choldx     (  0 -1  2 ) = (  0.00 -0.81  1.15 )
--
-- Given a positive-deﬁnite symmetric matrix a[1..n][1..n],
-- this routine constructs its Cholesky decomposition,
-- A = L · L^T
-- The Cholesky factor L is returned in the lower triangle of a,
-- except for its diagonal elements which are returned in p[1..n].


choldc :: forall a. Cov a -> Jac a a
choldc (Cov {vc= a}) = Jac {vj= a'} where
  ll = A.length a
  n = case ll of
        6  -> 3
        10 -> 4
        15 -> 5
  a' = A.create $ do -- make a Array of n x n +n space for diagonal +1 for summing
    arr <- MA.new $ (ll+n+1)
    -- loop over input array using Numerical Recipies algorithm (chapter 2.9)
    let ixa = indVs n
        ixarr = indV n
    numLoop 0 (n-1) $ \i0 -> do
      numLoop i0 (n-1) $ \j0 -> do
        let aij = uidx a (ixa i0 j0)
        _ <- if i0==j0 then MA.unsafeWrite arr (ll + i0) aij
                         else MA.unsafeWrite arr (ixarr j0 i0) aij
        numLoop 0 i0 $ \k0 -> do
          maik <- MA.unsafeRead arr (ixarr i0 k0)
          majk <- MA.unsafeRead arr (ixarr j0 k0)
          maij <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                            else MA.unsafeRead arr (ixarr j0 i0)
          let sum = (maij) - (maik) * (majk)
          _ <- if i0==j0 then MA.unsafeWrite arr (ll+i0) sum
                           else MA.unsafeWrite arr (ixarr j0 i0) sum
          msum <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                            else MA.unsafeRead arr (ixarr j0 i0)
          let 
              sum = if i0==j0 && msum < 0.0
                        then error ("choldInv: not a positive definite matrix "
                                     <> show a)
                        else msum
          p_i' <- MA.unsafeRead arr (ll+i0)
          let 
              p_i = if i0 == j0 then sqrt sum else p_i'
          _ <- if i0==j0 then MA.unsafeWrite arr (ll+i0) p_i
                         else MA.unsafeWrite arr (ixarr j0 i0) (sum/p_i)
          pure ()

     -- copy diagonal back into array
    numLoop 0 (n-1) $ \i0 -> do
      aii <- MA.unsafeRead arr (ll+i0)
      _ <- MA.unsafeWrite arr (ixarr i0 i0) aii
      pure ()
    pure $ MA.unsafeTake ll arr

-- | Matrix inversion using Cholesky decomposition
-- | based on Numerical Recipies formula in 2.9
--
cholInv :: forall a. Cov a -> Cov a
cholInv (Cov {vc= a}) = Cov {vc= a'} where
  ll = A.length a
  n = case ll of
        6  -> 3
        10 -> 4
        15 -> 5
  l = A.create $ do -- make a Array of n x n +n space for diagonal +1 for summing
    arr <- MA.new $ (ll+n+1)
    -- loop over input array using Numerical Recipies algorithm (chapter 2.9)
    let ixa = indVs n
        ixarr = indV n
    numLoop 0 (n-1) $ \i0 -> do
      numLoop i0 (n-1) $ \j0 -> do
        let aij = uidx a (ixa i0 j0)
        _ <- if i0==j0 then MA.unsafeWrite arr (ll + i0) aij
                         else MA.unsafeWrite arr (ixarr j0 i0) aij
        numLoop 0 i0 $ \k0 -> do
          maik <- MA.unsafeRead arr (ixarr i0 k0)
          majk <- MA.unsafeRead arr (ixarr j0 k0)
          maij <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                            else MA.unsafeRead arr (ixarr j0 i0)
          let sum = (maij) - (maik) * (majk)
          _ <- if i0==j0 then MA.unsafeWrite arr (ll+i0) sum
                           else MA.unsafeWrite arr (ixarr j0 i0) sum
          msum <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                            else MA.unsafeRead arr (ixarr j0 i0)
          let 
              sum = if i0==j0 && msum < 0.0
                        then error ("choldInv: not a positive definite matrix "
                                     <> show a)
                        else msum
          p_i' <- MA.unsafeRead arr (ll+i0)
          let 
              p_i = if i0 == j0 then sqrt sum else p_i'
          _ <- if i0==j0 then MA.unsafeWrite arr (ll+i0) p_i
                         else MA.unsafeWrite arr (ixarr j0 i0) (sum/p_i)
          pure ()

     -- invert L -> L^(-1)
    numLoop 0 (n-1) $ \i0 -> do
      mp_i <- MA.unsafeRead arr (ll+i0)
      _ <- MA.unsafeWrite arr (ixarr i0 i0) (1.0/mp_i)
      numLoop (i0+1) (n-1) $ \j0 -> do
        _ <- MA.unsafeWrite arr (ll+n) 0.0
        numLoop i0 j0 $ \k0 -> do
          majk <- MA.unsafeRead arr (ixarr j0 k0)
          maki <- MA.unsafeRead arr (ixarr k0 i0)
          sum <- MA.unsafeRead arr (ll+n)
          _ <- MA.unsafeWrite arr (ll+n) ((sum) - (majk) * (maki))
          pure ()
        msum <- MA.unsafeRead arr (ll+n)
        mp_j <- MA.unsafeRead arr (ll+j0)
        _ <- MA.unsafeWrite arr (ixarr j0 i0) ((msum)/(mp_j))
        pure ()
    pure arr
      
  a' = A.create $ do
    v <- MA.new $ n * (n+1) `div` 2
    let ixa = indVs n
        ixb = indVs n
        ixc = indVs n
    numLoop 0 (n-1) $ \i0 -> 
      numLoop i0 (n-1) $ \j0 -> 
        MA.unsafeWrite v (ixc i0 j0) $
          sum [ (uidx l (ixa k0 i0)) * (uidx l (ixb k0 j0)) 
               | k0 <- [0 .. n-1] ]
    pure v


--C version Numerical Recipies 2.9
--for (i=1;i<=n;i++) {
--  for (j=i;j<=n;j++) {
--    for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
--    if (i == j) {
--      if (sum <= 0.0) nrerror("choldc failed");
--      p[i]=sqrt(sum);
--    } else a[j][i]=sum/p[i];
--  }
--}
-- In this, and many other applications, one often needs L^(−1) . The lower
-- triangle of this matrix can be efﬁciently found from the output of choldc:
--for (i=1;i<=n;i++) {
--  a[i][i]=1.0/p[i];
--  for (j=i+1;j<=n;j++) {
--    sum=0.0;
--    for (k=i;k<j;k++) sum -= a[j][k]*a[k][i];
--    a[j][i]=sum/p[j];
--  }
--}

testCov2 :: String
testCov2 = s where
  xc3 :: Cov Dim3
  xc3 = Cov {vc= [1.0 .. 6.0]}
  xj3 :: Jac Dim3 Dim3
  xj3 = Jac {vj= [1.0 .. 9.0]}
  xj31 :: Jac Dim3 Dim3
  xj31 = Jac {vj= [1.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0]}
  xj32 :: Jac Dim3 Dim3
  xj32 = Jac {vj= [0.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,0.0]}
  xj33 :: Jac Dim3 Dim3
  xj33 = Jac {vj= [1.0 .. 9.0]}
  xj53 :: Jac Dim5 Dim3
  xj53 = Jac {vj= [1.0 .. 15.0]}
  xvc3 = toArray xc3
  xv3 = fromArray [1.0,1.0,1.0] :: Vec3
  xv5 = fromArray [1.0,1.0,1.0,1.0,1.0] :: Vec5
  s =  "Test Cov 2----------------------------------------------\n"
    <> "Vec *. Vec = " <> show (v3 *. v3) <> "\n"
--    <> "Cov *. Cov = " <> show ((one::Cov3) *. inv (one::Cov3)) <> "\n"
--    <> "Vec + Vec = " <> show (v5 + v5) <> "\n"
--    <> "chol Cov = " <> show (chol (one::Cov5)) <> "\n"
--    <> "Vec .*. Cov = " <> show (v5 .*. inv (one::Cov5)) <> "\n"
    <> "xc3 :: Cov Dim3 " <> show xc3
    <> show (toArray $ xc3) <> "\n"
    <> "xj53 ---> " <> show xj53
    <> "xj53 *. xc3 ---> " <> show (xj53 *. xc3)
    <> show xvc3 <> "\n"
        {-- <> show md <> "\n" --}
        {-- <> show mm3 --}
        {-- <> show mm5 --}
        {-- <> "exp v3 " <> show ( (v3 + v3) |.| v3 ) <> "\n" --}
        {-- <> show (j53) --}
        {-- <> show (tr j53) --}
        {-- <> "tj3 " <> show tj3 --}
        {-- <> "vv3 " <> show vv3 --}
        {-- <> show (v3 |*| c3) --}
        {-- <> "\n(tr j53 .*. c3)" <> show (tr j53 .*. c3) --}
        {-- <> "(tr j53 ||| v5)" <> show (tr j53 ||| v5) --}
        {-- <> show (c3 ** (inv c3)) --}
        {-- <> show (c4 ** (inv c4)) --}
    <> "chol: -----------------\n"
    <> "A = L * L^T         " <> show ch3
    <> "L                   " <> show (choldc ch3)
        {-- <> "L * L^T             " <> show ((choldc ch3) *. tr (choldc ch3)) --}
    <> "A^(-1) = L' * L'^T  " <> show (inv ch3)
    <> "A^(-1) from cholInv " <> show (cholInv ch3)
    <> "A = L * L^T         " <> show ch5
    <> "L                   " <> show (choldc ch5)
        {-- <> "L * L^T             " <> show ((choldc ch5) *. tr (choldc ch5)) --}
    <> "A^(-1) = L' * L'^T  " <> show (inv ch5)
    <> "A^(-1) from cholInv " <> show (cholInv ch5)
    <> "det this            " <> show (det ch5)
        {-- <> "chol" <> show ch5 --}
        {-- <> show (choldc ch5) <> show ( choldc ch5 *. tr (choldc ch5)) --}
    <> "\n" -- <> testCov2 
  c3 :: Cov3
  c3 = fromArray [1.0,2.0,3.0,4.0,5.0,6.0]
  c4 :: Cov4
  c4 = fromArray [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
  c5 :: Cov5
  c5 = fromArray [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]
  c50 :: Cov5
  c50 = fromArray [15.0,14.0,13.0,12.0,11.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0]
  c50m :: Cov5
  c50m = fromArray $ M.toArray $ toMatrix c50
  c51 :: Cov5
  c51 = one
  v3 :: Vec3
  v3 = fromArray [10.0,11.0,12.0]
  v5 :: Vec5
  v5 = fromArray [10.0,11.0,12.0,13.0,14.0]
  j53 :: Jac53
  j53 = fromArray [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]
--  tj3 :: Cov3
--  tj3 = j53 .*. c5
--  vv5 :: Vec5
--  vv5 = j53 *. v3
--  vv3 :: Vec3
--  vv3 = tr j53 *. j53 *. c3 *. v3

  m3 :: M.Matrix
  m3 = M.fromArray2 3 3 [1.0,2.0,3.0,2.0,4.0,5.0,3.0,5.0,6.0]
--  mm3 = (m3+m3)*m3
  m5 :: M.Matrix
  m5 = M.fromArray2 5 5 [1.0,2.0,3.0,4.0,5.0, 2.0,6.0,7.0,8.0,9.0
                        ,3.0,7.0,10.0,11.0,12.0, 4.0,8.0,11.0,13.0,14.0
                        ,5.0,9.0,12.0,14.0,15.0]
--  mm5 = (m5+m5)*m5
  ch3 :: Cov3
  ch3 = fromArray [2.0, -1.0, 0.0, 2.0, -1.0, 2.0]
  cch3 = choldc ch3
  ich3 = cholInv ch3

  ch5 :: Cov5
  ch5 = fromArray [2.0, -1.0, 0.0, 0.0, 0.0, 2.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0]
