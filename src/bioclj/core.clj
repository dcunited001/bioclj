(ns bioclj.core
  (:require [clojure.math.combinatorics :as com]))

(defn combinate-domain
  "Permutes a string or array of chars and returns an array of strings representing the domain"
  [chars]
  ;; the algorithm this way is very inefficient.  hash is always most efficient.
  ;; plus, clojure already has the com/permutate function
  (com/selections chars (count chars)))

(defn kmer-frequency
  "returns a hash where keys are kmers and values are frequency.  accepts k and text."
  [k text]
  (persistent!
    (let [idxstop (- (count text) k)]
      (reduce
        #(let [kmer (subs text %2 (+ %2 k))
               cnt (%1 kmer 0)]
          (assoc! %1 kmer (inc cnt)))
        (transient {})
        (range 0 (inc idxstop))
        ))))

(defn sorted-freq-kmer
  [freqmap]
  (into (sorted-map-by #(compare [(get freqmap %2) %2]
                                 [(get freqmap %1) %1])) freqmap))