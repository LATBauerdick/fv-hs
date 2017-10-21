
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Jac where
  import Prelude
  import Stuff
  data Jac a b = Jac { v :: Array Number, nr :: Int }
