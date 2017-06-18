module LibSpec where

import Test.Hspec ( Spec, hspec, describe, it, shouldBe )
import Test.QuickCheck ( property )
--import Control.Exception (evaluate)
import Control.Monad ( (<=<) )
import Data.List ( sort )
import Data.Maybe ( mapMaybe )

import Test.Input ( hSlurp )
-- import Test.Cluster ( doCluster, fsmw )
import Test.FVT ( testFVT )
import Data.Cov ( testCov2 )
import FV.Types ( VHMeas (..), HMeas (..), MCtruth (..) )


import Stuff

main :: IO ()
main = hspec spec

spec :: Spec
spec =
  describe "Lib Tests" $ do
    describe "Cov" $ do
      it "testCov works" $ do
        let  s = testCov2
        putStrLn s
        -- putStrLn $ testCov 0
        (head s) `shouldBe` 'T'

    -- describe "FSMW Tests" $ do
    --   it "FSMW5 works" $ do
    --     let xs = [5.0, 6.2, 7.0, 8.0, 9.0]
    --         n = length xs
    --     fsmw n xs `shouldBe` 6.6

    --   it "FSMW4 works" $ do
    --     let xs = [5.0, 6.2, 7.0, 8.0]
    --         n = length xs
    --     fsmw n xs `shouldBe` 6.6

    --   it "FSMW3 works" $ do
    --     let xs = [5.0, 6.2, 7.0]
    --         n = length xs
    --     fsmw n xs `shouldBe` 6.6

    --   it "FSMW property: returns value in bound" $ do
    --     let fff :: [Double] -> Bool
    --         fff [] = True
    --         fff xs = f >= head sxs && f <= last sxs where
    --           sxs = sort xs
    --           n = length xs
    --           f = fsmw n sxs -- `debug` show n
    --     property $ fff


    --   it "FSMW works" $ do
    --     let xs = [5.0, 6.0, 7.0, 7.1, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 14.1] ++ [15.0 .. 26.0] ++ [26.05, 26.1, 26.15, 26.17, 26.2, 26.3, 26.4] ++ [27.0..127.0]
    --         n = length xs
    --     fsmw n xs `shouldBe` 26.16

    describe "FV Test" $ do
      -- it "doCluster works" $ do
      --   _ <- doCluster . fst <=< hSlurp $ "dat/tav-0.dat"
      --   (1 :: Int) `shouldBe` 1

      it "hSlurp works" $ do
        ds <- readFile "dat/tav-4.dat"
        let VHMeas _ hl = uJust $ hSlurp ds
            HMeas _ _ w = head hl
        w `shouldBe` 0.0114

      it "hSlurp works" $ do
        ds <- readFile "dat/tr05129e001412.dat"
        let VHMeas _ hl = uJust <<< hSlurp $ ds
            HMeas _ _ w = head hl
        w `shouldBe` 4.5451703e-3

      it "test 1 works" $ do
        ds <- readFile "dat/tr05129e001412.dat"
        _ <- testFVT [0,2,3,4,5] <<< uJust <<< hSlurp $ ds
        (1 :: Int) `shouldBe` 1

      -- it "test p works" $ do
      --   _ <- test ["p"]
      --   (1 :: Int) `shouldBe` 1

--    prop "ourAdd is commutative" $ \x y ->
--      ourAdd x y `shouldBe` ourAdd y x
