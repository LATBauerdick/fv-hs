module LibSpec where

import Test.Hspec
import Test.QuickCheck
import Control.Exception (evaluate)

import FVT.Input (hSlurp)
import FV.Types ( VHMeas (..), HMeas (..) )

main :: IO ()
main = hspec spec

spec :: Spec
spec =
  describe "FVT.Test" $ do
    it "works" $ do
      VHMeas _ hl <- hSlurp "dat/tav-0.dat"
      let HMeas _ _ w = head hl
      w `shouldBe` 1.0
--    prop "ourAdd is commutative" $ \x y ->
--      ourAdd x y `shouldBe` ourAdd y x
