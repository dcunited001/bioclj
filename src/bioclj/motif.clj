(ns bioclj.motif
  (:require [clojure.core.reducers :as r]
            [clj-biosequence.core :as bio.core]
            [clj-biosequence.alphabet :as bio.alpha]
            [clojure.math.numeric-tower :as maths]))

(defn kmer-occurances
  "outputs the average number of occurrences of a random kmer within n random {ACTG} strings of L length"
  [k n L]
  (let [num-positions (inc (- L k))]
    (/ (* n num-positions 1.0) (maths/expt 4 k))))


(defn motif-enumeration
  []
  )

(defn score-motif
  [motifs]
  )

(defn count-motif
  [motifs]
  )

(defn profile-motif
  [motifs]
  )

