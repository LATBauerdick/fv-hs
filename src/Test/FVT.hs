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

module Test.FVT (testFVT) where

import Prelude.Extended
-- import Control.Monad.Eff (Eff)
-- import Control.Monad.Eff.Console (CONSOLE, log, logShow)
-- import Control.Monad.Eff.Random ( RANDOM )
import Data.Monoid ( mempty )
-- import Data.Tuple ( Tuple(..) )
-- import Data.Array ( length, zip, foldl )
import Data.Foldable ( traverse_ )

import FV.Types
  ( VHMeas, HMeas, QMeas
  , XMeas, Prong (..), Chi2 (Chi2)
  , vertex, helices,  hFilter, fromHMeas, fromQMeas, vBlowup, distance, invMass
  )

import FV.Fit ( fit )

type String = Text

showMomentum :: HMeas -> String
showMomentum h = "pt,pz,fi,E ->" <> (pack <<< show <<< fromHMeas) h
showHelix :: HMeas -> String
showHelix h = "Helix ->" <> (pack <<< show) h
showProng :: Prong -> String
showProng (Prong {nProng= n, fitVertex= v, fitMomenta= ql, fitChi2s= cl}) =
  let
      showCl :: String -> List Chi2 -> String
      showCl = foldl (\s (Chi2 x) -> s <> to1fix x)
      Chi2 chi2tot = sum cl
      sc = "chi2tot ->" <> to1fix chi2tot <> ", ndof " <> (pack <<< show $ n*2)
      sd = ", r ->" <> (pack <<< show $ distance v mempty)
      scl = showCl ", chi2s ->" cl
      sm = ", Mass ->" <> (pack <<< show $ invMass (map fromQMeas ql))
  in sc <> sd <> scl <> sm

testFVT :: List Int -> VHMeas -> IO ()
testFVT l5 vm = do
  let hel = helices vm
  traverse_ (putStrLn <<< showHelix) hel
  traverse_ (putStrLn <<< showMomentum) hel
  doFitTest vm l5
  putStrLn $ showProng <<< fit <<< hFilter l5 <<< vBlowup 10000.0 $ vm
  pure ()

doFitTest :: VHMeas
            -> List Int
            -> IO ()
doFitTest vm' l5 = do
  let vm = vBlowup 10000.0 vm'
  let showLen xs = show $ length xs
      showQChi2 :: (QMeas, Chi2) -> String
      showQChi2 (qm, (Chi2 chi2)) = "q"
                                <> " chi2 ->" <> to1fix chi2
                                <> " pt,pz,fi,E ->"
                                <> (pack <<< show) qm

  putStrLn $           "initial vertex position -> " <> show ((vertex vm)::XMeas)

  let pl         = map (fromQMeas <<< fromHMeas) $ helices vm
  putStrLn $ "Inv Mass " <> showLen pl <> " helix" <> show (invMass pl)
  let pl5        = map (fromQMeas <<< fromHMeas) (helices <<< hFilter l5 $ vm)
  putStrLn $ "Inv Mass " <> showLen pl5 <> " helix" <> show (invMass pl5)

  putStrLn $           pack "Fitting Vertex --------------------"
  let -- pr = fit vm
      Prong {fitVertex= vf, fitMomenta= ql, fitChi2s= cl} = fit vm
  putStrLn $           "Fitted vertex -> " <> show vf
  traverse_ (putStrLn <<< showQChi2) $ zip ql cl
  putStrLn $ "Inv Mass " <> show (length ql) <> " fit"
                    <> show (invMass (map fromQMeas ql))

  let m5 = invMass <<< map fromQMeas <<< iflt l5 $ ql
  putStrLn $ "Inv Mass " <> show (length l5) <> " fit" <> show m5

  putStrLn $           pack "Refitting Vertex-----------------"
  let Prong {fitVertex=fv, fitMomenta=fqs, fitChi2s=fcs, nProng=np} = fit <<< hFilter l5 $ vm
  putStrLn $           "Refitted vertex -> " <> show fv
  traverse_ (putStrLn <<< showQChi2) $ zip fqs fcs
  putStrLn $           "Inv Mass " <> (pack <<< show) np <> " refit" 
                       <> (pack <<< show <<< invMass <<< map fromQMeas $ fqs)
  putStrLn $           "Final vertex -> " <> (pack <<< show) fv
  putStrLn $           pack "end of doFitTest------------------------------------------"

