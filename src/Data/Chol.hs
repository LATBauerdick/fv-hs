--{-# LANGUAGE EmptyDataDecls #-}
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

{-# LANGUAGE DisambiguateRecordFields #-}
{-# LANGUAGE NamedFieldPuns #-}

module Data.Chol ( doCholdc, doCholInv )
  where

import Prelude.Extended

import qualified Data.Vector.Unboxed as A ( create )
import qualified Data.Vector.Unboxed.Mutable as MA (
  new, unsafeWrite, unsafeRead, unsafeTake )
import Control.Monad ( void )
import Control.Loop ( numLoop )

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

doCholdc :: Array Number -> Int -> Array Number
doCholdc a n = a' where
  ll = n*n
  a' = A.create $ do -- make a Array of n x n +n space for diagonal +1 for summing
    arr <- MA.new $ (ll+n+1)
    -- loop over input array using Numerical Recipies algorithm (chapter 2.9)
    let ixa = indVs n
        ixarr = indV n
    numLoop 0 (n-1) $ \i0 -> do
      numLoop i0 (n-1) $ \j0 -> do
        let aij = uidx a (ixa i0 j0)
        void $ if i0==j0 then MA.unsafeWrite arr (ll + i0) aij
                        else MA.unsafeWrite arr (ixarr j0 i0) aij
        numLoop 0 i0 $ \k0 -> do
          aik <- MA.unsafeRead arr (ixarr i0 k0)
          ajk <- MA.unsafeRead arr (ixarr j0 k0)
          maij <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                            else MA.unsafeRead arr (ixarr j0 i0)
          let s = maij - aik * ajk
          void $ if i0==j0 then MA.unsafeWrite arr (ll+i0) s
                          else MA.unsafeWrite arr (ixarr j0 i0) s
        msum <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                         else MA.unsafeRead arr (ixarr j0 i0)
        let s = if i0==j0 && msum < 0.0
                        then error ("choldc: not a positive definite matrix " <> show a)
                        else msum
        p_i' <- MA.unsafeRead arr (ll+i0)
        let p = if (i0 == j0) then (sqrt s) else s/p_i'
        void $ if i0==j0 then MA.unsafeWrite arr (ll+i0) p
                        else MA.unsafeWrite arr (ixarr j0 i0) p
        pure ()

    -- copy diagonal back into array
    numLoop 0 (n-1) $ \i0 -> do
      aii <- MA.unsafeRead arr (ll+i0)
      void $ MA.unsafeWrite arr (ixarr i0 i0) aii
      pure ()

    pure $ MA.unsafeTake ll arr

-- | Matrix inversion using Cholesky decomposition
-- | based on Numerical Recipies formula in 2.9
--
doCholInv :: Array Number -> Int -> Array Number
doCholInv a n = a' where
  ll = n*n
  l = A.create $ do -- make a Array of n x n +n space for diagonal +1 for summing
    arr <- MA.new $ (ll+n+1)
    -- loop over input array using Numerical Recipies algorithm (chapter 2.9)
    let ixa = indVs n
        ixarr = indV n
    numLoop 0 (n-1) $ \i0 -> do
      numLoop i0 (n-1) $ \j0 -> do
        let aij = uidx a (ixa i0 j0)
        void $ if i0==j0 then MA.unsafeWrite arr (ll + i0) aij
                      else MA.unsafeWrite arr (ixarr j0 i0) aij
        numLoop 0 i0 $ \k0 -> do
          aik <- MA.unsafeRead arr (ixarr i0 k0)
          ajk <- MA.unsafeRead arr (ixarr j0 k0)
          maij <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                            else MA.unsafeRead arr (ixarr j0 i0)
          let s = maij - aik * ajk
          void $ if i0==j0 then MA.unsafeWrite arr (ll+i0) s
                          else MA.unsafeWrite arr (ixarr j0 i0) s
        msum <- if i0==j0 then MA.unsafeRead arr (ll+i0)
                         else MA.unsafeRead arr (ixarr j0 i0)
        let s = if i0==j0 && msum < 0.0
                        then error ("cholInv: not a positive definite matrix "
                                     <> show a)
                        else msum
        p_i' <- MA.unsafeRead arr (ll+i0)
        let p = if i0 == j0 then sqrt s else s/p_i'
        void $ if i0 ==j0 then MA.unsafeWrite arr (ll+i0) p
                        else MA.unsafeWrite arr (ixarr j0 i0) p
        pure ()

    -- invert L -> L^(-1)
    numLoop 0 (n-1) $ \i0 -> do
      p_i <- MA.unsafeRead arr (ll+i0)
      void $ MA.unsafeWrite arr (ixarr i0 i0) (1.0/p_i)
      numLoop (i0+1) (n-1) $ \j0 -> do
        void $ MA.unsafeWrite arr (ll+n) 0.0
        numLoop i0 j0 $ \k0 -> do
          ajk <- MA.unsafeRead arr (ixarr j0 k0)
          aki <- MA.unsafeRead arr (ixarr k0 i0)
          s <- MA.unsafeRead arr (ll+n)
          void $ MA.unsafeWrite arr (ll+n) (s - ajk * aki)
        msum <- MA.unsafeRead arr (ll+n)
        p_j <- MA.unsafeRead arr (ll+j0)
        void $ MA.unsafeWrite arr (ixarr j0 i0) (msum/p_j)
    pure arr

  a' = fromList $ do
    let idx = indV n
    i0 <- range 0 (n-1)
    j0 <- range i0 (n-1)
    let aij = sum $ do
                  k0 <- range 0 (n-1)
                  pure $ (uidx l (idx k0 i0)) * (uidx l (idx k0 j0))
    pure $ aij

--C version Numerical Recipies 2.9
--for (i=1;i<=n;i++) {
--  for (j=i;j<=n;j++) {
--    for (s=a[i][j],k=i-1;k>=1;k--) s -= a[i][k]*a[j][k];
--    if (i == j) {
--      if (s <= 0.0) nrerror("choldc failed");
--      p[i]=sqrt(s);
--    } else a[j][i]=s/p[i];
--  }
--}
-- In this, and many other applications, one often needs L^(−1) . The lower
-- triangle of this matrix can be efﬁciently found from the output of choldc:
--for (i=1;i<=n;i++) {
--  a[i][i]=1.0/p[i];
--  for (j=i+1;j<=n;j++) {
--    s=0.0;
--    for (k=i;k<j;k++) s -= a[j][k]*a[k][i];
--    a[j][i]=s/p[j];
--  }
--}

