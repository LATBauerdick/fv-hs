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

import FV.Types
  ( VHMeas, HMeas, QMeas
  , XMeas, Prong (..), Chi2 (Chi2)
  , vertex, helices,  hFilter, fromHMeas, fromQMeas, vBlowup, distance, invMass
  )

import FV.Fit ( fit )

showMomentum :: HMeas -> Text
showMomentum h = "pt,pz,fi,E ->" <> (tshow <<< fromHMeas) h
showHelix :: HMeas -> Text
showHelix h = "Helix ->" <> tshow h
showProng :: Prong -> Text
showProng (Prong {nProng= n, fitVertex= v, fitMomenta= ql, fitChi2s= cl}) =
  let
      showCl :: Text -> List Chi2 -> Text
      showCl = foldl (\s (Chi2 x) -> s <> to1fix x)
      Chi2 chi2tot = sum cl
      sc = "chi2tot ->" <> to1fix chi2tot <> ", ndof " <> (tshow $ n*2)
      sd = ", r ->" <> (tshow $ distance v mempty)
      scl = showCl ", chi2s ->" cl
      sm = ", Mass ->" <> (tshow $ invMass (map fromQMeas ql))
  in sc <> sd <> scl <> sm

testFVT :: List Int -> VHMeas -> IO ()
testFVT l5 vm = do
  let hel = helices vm
  traverse_ (putTextLn <<< showHelix) hel
  traverse_ (putTextLn <<< showMomentum) hel
  doFitTest vm l5
  putTextLn $ showProng <<< fit <<< hFilter l5 <<< vBlowup 10000.0 $ vm
  pure ()

doFitTest :: VHMeas
            -> List Int
            -> IO ()
doFitTest vm' l5 = do
  let vm = vBlowup 10000.0 vm'
  let showLen xs = tshow $ length xs
      showQChi2 :: (QMeas, Chi2) -> Text
      showQChi2 (qm, (Chi2 chi2)) = "q"
                                <> " chi2 ->" <> to1fix chi2
                                <> " pt,pz,fi,E ->"
                                <> tshow qm

  putTextLn $           "initial vertex position -> " <> tshow ((vertex vm)::XMeas)

  let pl         = map (fromQMeas <<< fromHMeas) $ helices vm
  putTextLn $ "Inv Mass " <> showLen pl <> " helix" <> tshow (invMass pl)
  let pl5        = map (fromQMeas <<< fromHMeas) (helices <<< hFilter l5 $ vm)
  putTextLn $ "Inv Mass " <> showLen pl5 <> " helix" <> tshow (invMass pl5)

  putTextLn $           "Fitting Vertex --------------------"
  let -- pr = fit vm
      Prong {fitVertex= vf, fitMomenta= ql, fitChi2s= cl} = fit vm
  putTextLn $           "Fitted vertex -> " <> tshow vf
  traverse_ (putTextLn <<< showQChi2) $ zip ql cl
  putTextLn $ "Inv Mass " <> tshow (length ql) <> " fit"
                    <> tshow (invMass (map fromQMeas ql))

  let m5 = invMass <<< map fromQMeas <<< iflt l5 $ ql
  putTextLn $ "Inv Mass " <> tshow (length l5) <> " fit" <> tshow m5

  putTextLn $           "Refitting Vertex-----------------"
  let Prong {fitVertex=fv, fitMomenta=fqs, fitChi2s=fcs, nProng=np} = fit <<< hFilter l5 $ vm
  putTextLn $           "Refitted vertex -> " <> tshow fv
  traverse_ (putTextLn <<< showQChi2) $ zip fqs fcs
  putTextLn $           "Inv Mass " <> (tshow) np <> " refit" 
                       <> (tshow <<< invMass <<< map fromQMeas $ fqs)
  putTextLn $           "Final vertex -> " <> tshow fv
  putTextLn $           "end of doFitTest------------------------------------------"

