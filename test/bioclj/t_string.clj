(ns bioclj.t-string
  (:require [clojure.java.io :as io])
  (:use [midje.sweet]
        [bioclj.string]
        [gloss core io]))

(defn read-bytes [file]
  (with-open [input (io/input-stream file)
              output (new java.io.ByteArrayOutputStream)]
    (io/copy input output)
    (.toByteArray output)))

(facts "actg-codec"
       (let [acgt-16 (gloss.io/decode acgt-codec (read-bytes "test/data/ACGTx16.bin"))
             ;acgt-1 (gloss.io/decode acgt-codec (read-bytes "test/data/ACGTx1.bin"))
             ;acgt-9 (gloss.io/decode acgt-codec (read-bytes "test/data/ACGTx9.bin"))
             ]

         (fact "can read bytes when num bytes % 8 = 0")
         (fact "can read bytes when num bytes % 8 != 0")
       (fact "can convert strings to 64b integer array")))

