
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Cov where
import Prelude.Extended
import qualified Data.Vector.Unboxed as A (
    length, singleton, map, zipWith )

newtype Cov a = Cov { v :: Array Number }

instance Num (Cov a) where
  fromInteger i = Cov {v= A.singleton <<< fromInteger $ i}
  negate (Cov {v=v_}) = Cov {v=A.map negate v_}
  abs (Cov {v=v_}) = Cov {v=A.map abs v_}
  signum (Cov {v=v_}) = Cov {v=A.map signum v_}
  (+) (Cov {v= va}) (Cov {v= vb}) = (Cov {v= vc}) where vc = A.zipWith (+) va vb
  (*) = error "cannot multiply Cov*Cov to return a Cov, use *. instead"

instance Show (Cov a) where
  show c = toString $ ( "Show (Cov a) \n") <> showMatrix c where
    showMatrix (Cov {v=v_}) = let
      makeSymMat :: Int -> Array Number -> Array Number
      makeSymMat n vs = fromList $ do
        let iv = indVs n
        i <- range 0 (n-1)
        j <- range 0 (n-1)
        pure $ uidx vs (iv i j)
      in
        case A.length v_ of
                            6  -> prettyMatrix 3 3 $ makeSymMat 3 v_
                            10 -> prettyMatrix 4 4 $ makeSymMat 4 v_
                            15 -> prettyMatrix 5 5 $ makeSymMat 5 v_
                            _ -> error $ "showCova showMatrix "
