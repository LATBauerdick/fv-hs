module LibSpec where

import Prelude.Extended
import Unsafe ( unsafeHead )
import qualified Data.Text as T ( head ) 

import Test.Hspec ( Spec, hspec, describe, it, shouldBe )
import Test.QuickCheck ( property )
--import Control.Exception (evaluate)
-- import Control.Monad ( (<=<) )
import Data.List ( sort )
-- import Data.Maybe ( mapMaybe )
-- import Data.Foldable ( traverse_ )

import Test.Input ( hSlurp )
import Test.Cluster ( doCluster, fsmw )
import Test.FVT ( testFVT )
import Data.Cov ( testCov2 )
import FV.Types ( VHMeas (..), HMeas (..), hFilter, vBlowup, fromHMeas )
import Test.Random ( testRandom )

main :: IO ()
main = hspec spec

spec :: Spec
spec =
  describe "Lib Tests" $ do
    describe "Cov" $ do
      it "testCov------------------------" $ do
        let  s = testCov2
        Prelude.Extended.putStrLn s
        -- Prelude.Extended.putStrLn $ testCov 0
        T.head s `shouldBe` 'T'

    describe "FSMW Tests" $ do
      it "FSMW5 works" $ do
        let xs = [5.0, 6.2, 7.0, 8.0, 9.0]
            n = length xs
        fsmw n xs `shouldBe` 6.6

      it "FSMW4 works" $ do
        let xs = [5.0, 6.2, 7.0, 8.0]
            n = length xs
        fsmw n xs `shouldBe` 6.6

      it "FSMW3 works" $ do
        let xs = [5.0, 6.2, 7.0]
            n = length xs
        fsmw n xs `shouldBe` 6.6

      it "FSMW property: returns value in bound" $ do
        let fff :: [Double] -> Bool
            fff [] = True
            fff xs = f >= unsafeHead sxs && f <= last sxs where
              sxs = sort xs
              n = length xs
              f = fsmw n sxs -- `debug` show n
        property $ fff


      it "FSMW works" $ do
        let xs = [5.0, 6.0, 7.0, 7.1, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 14.1] ++ [15.0 .. 26.0] ++ [26.05, 26.1, 26.15, 26.17, 26.2, 26.3, 26.4] ++ [27.0..127.0]
            n = length xs
        fsmw n xs `shouldBe` 26.16

    describe "FV Test" $ do
      it "hSlurp CMS event" $ do
        ds <- Prelude.Extended.readFile "dat/tav-4.dat"
        let VHMeas _ hl = uJust $ hSlurp ds
            HMeas _ _ w = unsafeHead hl
        w `shouldBe` 0.0114

      it "hSlurp Aleph event" $ do
        ds <- Prelude.Extended.readFile "dat/tr05129e001412.dat"
        let VHMeas _ hl = uJust <<< hSlurp $ ds
            HMeas _ _ w = unsafeHead hl
        w `shouldBe` 4.5451703e-3

      it "--Test FVT 1" $ do
        _ <- testFVT [0,2,3,4,5] <<< uJust <<< hSlurp =<< Prelude.Extended.readFile "dat/tr05129e001412.dat"
        1 `shouldBe` 1

      -- it "--Test FVT 2" $ do
      --   ds <- readFile "dat/tr05158e004656.dat"
      --   _ <- testFVT [0,1,2,4,5] <<< uJust <<< hSlurp $ ds
      --   1 `shouldBe` 1

      -- it "--Test FVT 3" $ do -- this file has an additional track that makes the fit bad
      --   ds <- readFile "dat/tr07849e007984.dat"
      --   _ <- testFVT [0,2,3,4,5] <<< uJust <<< hSlurp =<< readFile "dat/tr07849e007984.dat"
      --   ds <- readFile "dat/tr07849e007984.dat"
      --   let vm = uJust <<< hSlurp $ ds
      --   -- pr <- showProng . fitw . hFilter l5 . vBlowup 10000.0 $ vm
      --   -- let vf = fitVertex $ pr
      --   --     h = head . helices . hFilter [6] $ vm
      --   --     Just (_, chi2, _) =  ksm vf h
      --   -- Prelude.Extended.putStrLn $ printf "chi2 of track 6 w/r to fit vertex is %8.1f" (chi2::Double)
      --   1 `shouldBe` 1

      it "--Test Cluster" $ do
        let
            showMomentum :: HMeas -> String
            showMomentum h = "pt,pz,fi,E ->" <> (show <<< fromHMeas) h
            showHelix :: HMeas -> String
            showHelix h = "Helix ->" <> show h
        ds <- Prelude.Extended.readFile "dat/tav-4.dat"
        let vm = uJust $ hSlurp ds
        -- traverse_ (Prelude.Extended.putStrLn <<< showHelix) $ helices vm
        -- traverse_ (Prelude.Extended.putStrLn <<< showMomentum) $ helices vm
        s <- doCluster vm
        Prelude.Extended.putStrLn s
        -- _ <- doCluster . fst <=< hSlurp $ "dat/tav-0.dat"
        (1 :: Int) `shouldBe` 1

      it "--Test Random 5,000" $ do
        out <- testRandom 5 <<< hFilter [0,2,3,4,5] <<< vBlowup 10000.0
                            <<< uJust <<< hSlurp
                            =<< Prelude.Extended.readFile "dat/tr05129e001412.dat"
        Prelude.Extended.putStrLn out
        (1 :: Int) `shouldBe` 1

      -- it "test p works" $ do
      --   _ <- test ["p"]
      --   (1 :: Int) `shouldBe` 1

--    prop "ourAdd is commutative" $ \x y ->
--      ourAdd x y `shouldBe` ourAdd y x
