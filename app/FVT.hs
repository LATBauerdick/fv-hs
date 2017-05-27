-- file: test.hs
--
module Main ( main ) where

import Prelude
import System.Environment
import System.Exit

import Test.Test ( test )

main :: IO ()
main = getArgs >>= parse

parse :: [String] -> IO ()
parse ["-h"] = usage   >> exit
parse ["-v"] = version >> exit
parse []     = test ["1"]
parse args   = test args

usage :: IO ()
usage   = putStrLn "Usage: fvt [-vh] [test# ..]"
version :: IO ()
version = putStrLn "Haskell fvt 0.1"
exit :: IO ()
exit    = exitSuccess

