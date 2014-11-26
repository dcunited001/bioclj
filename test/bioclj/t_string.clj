(ns bioclj.t-string
  (:require [clojure.java.io :as io])
  (:use [midje.sweet]
        [bioclj.string]))

(facts "reduce-acgt-64b && acgt-64b-to-str"
  (let [seq1 "ACGTACGTACGTACGTACGTACGTACGTACGT"
        long1 1953184666628070171]
    (fact "should be reversible"
          (reduce-acgt-64b seq1) => long1
          (acgt-64b-to-str 32 long1) => (vec seq1))
    (fact "can encode nucleotide sequences of multiple lengths"
          (acgt-64b-to-str 4 long1) => (vec "ACGT")
          (acgt-64b-to-str 8 long1) => (vec "ACGTACGT"))))

(facts "acgt-to-64b"
  (let [seq1 "ACGTACGTACGTACGTACGTACGTACGTACGTTCGATCGATCGATCGATCGATCGATCGATCGAACGTACGT"
        longs1 (acgt-to-64b seq1)]
    (fact "should be reversible"
          (acgt-from-64b (first longs1) (rest longs1)) => (vec seq1))))

(facts "hamming-b64"
       (let [s1 (reduce-acgt-64b "ACGT")
             s2 (reduce-acgt-64b "AAAA")
             s3 (reduce-acgt-64b "AAAATTTT")
             s4 (reduce-acgt-64b "AATTTTGG")]
         (fact "can compare sequences of equal length accurately"
               (hamming-64b s1 s2) => 3
               (hamming-64b s3 s4) => 4)
         (fact "can compare sequences of varied length, given k"
               (hamming-64b 4 s2 s3) => 0
               (hamming-64b 4 s1 s4) => 2)))

(facts "neighborhood-64b-acgt")
