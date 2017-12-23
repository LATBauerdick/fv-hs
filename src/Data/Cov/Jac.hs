
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Jac where
  import Prelude.Extended
  data Jac a b = Jac { v :: Array Number, nr :: Int }
