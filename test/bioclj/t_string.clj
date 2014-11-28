(ns bioclj.t-string
  (:require [clojure.java.io :as io])
  (:use [midje.sweet]
        [bioclj.string]))

(facts "acgt-str-to-64b && acgt-64b-to-str"
       (let [seq1 "ACGTACGTACGTACGTACGTACGTACGTACGT"
             long1 1953184666628070171]
         (fact "should be reversible"
               (acgt-str-to-64b seq1) => long1
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
       (let [s1 (acgt-str-to-64b "ACGT")
             s2 (acgt-str-to-64b "AAAA")
             s3 (acgt-str-to-64b "AAAATTTT")
             s4 (acgt-str-to-64b "AATTTTGG")]
         (fact "can compare sequences of equal length accurately"
               (hamming-64b s1 s2) => 3
               (hamming-64b s3 s4) => 4)
         (fact "can compare sequences of varied length, given k"
               (hamming-64b 4 s2 s3) => 0
               (hamming-64b 4 s1 s4) => 2)))

(facts "neighborhood-64b-acgt"
       (let [s1 (acgt-str-to-64b "AA")
             n1 (sort (map acgt-str-to-64b '("AC" "AT" "AG" "CA" "TA" "GA" "AA")))
             s2 (acgt-str-to-64b "AAA")
             n2 (sort (map acgt-str-to-64b '("CAA" "TAA" "GAA" "AAA" "AGA" "ACA" "ATA" "AAG" "AAT" "AAC")))]
         (fact "generates a matching neighborhood of an ACGT string"
               (sort (:b64 (neighborhood-acgt-64b 2 1 s1))) => n1
               (sort (:b64 (neighborhood-acgt-64b 3 1 s2))) => n2)
         (fact "generates a neighborhood with the appropriate counts"
               (count (:b64 (neighborhood-acgt-64b 2 1 s1))) => 7
               (count (:b64 (neighborhood-acgt-64b 2 2 s1))) => 16
               (count (:b64 (neighborhood-acgt-64b 3 1 s2))) => 10
               (count (:b64 (neighborhood-acgt-64b 3 2 s2))) => 37
               (count (:b64 (neighborhood-acgt-64b 3 3 s2))) => 64)
         (fact "neighborhoods include the input string itself"
               (.contains (:b64 (neighborhood-acgt-64b 2 1 s1)) s1) => true
               (.contains (:b64 (neighborhood-acgt-64b 3 2 s2)) s2) => true)))

;(count  (neighborhood-acgt-64b 32 3 (reduce-acgt-64b "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")))
;=> 138481 @ 00:01.000
;(count  (neighborhood-acgt-64b 32 4 (reduce-acgt-64b "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")))
;=> 3051241 @ 00:12.000
;(count  (neighborhood-acgt-64b 32 5 (reduce-acgt-64b "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")))
;=> 51985609 @ 05:00.000
;; not bad
;; if i just had a GPU method for filtering nils from lists, i could get 32,5 down to seconds
