(ns bioclj.core
  (:require [clojure.math.combinatorics :as com]))

(defn combinate-domain
  "Permutes a string or array of chars and returns an array of strings representing the domain"
  [chars]
  ;; the algorithm this way is very inefficient.  hash is always most efficient.
  ;; plus, clojure already has the com/permutate function
  (com/selections chars (count chars)))

(defn kmer-op
  "inspects text, running a fn with all substrings of length k
  f takes 3 params [hash index substring]"
  [f k text]
  (persistent!
    (let [idxstop (- (count text) k)]
      (reduce
        #(let [kmer (subs text %2 (+ %2 k))]
          (assoc! %1 kmer (f %1 %2 kmer)))
        (transient {})
        (range 0 (inc idxstop))))))

(defn kmer-frequency
  "returns a hash where keys are kmers and values are frequency.  accepts k and text."
  [k text]
  (kmer-op #(inc (%1 %3 0)) k text))

(defn kmer-indices
  [k text]
  (kmer-op #(conj (%1 %3 []) %2) k text))

(defn sorted-freq-kmer
  [freqmap]
  (into (sorted-map-by #(compare [(get freqmap %2) %2]
                                 [(get freqmap %1) %1])) freqmap))

(def pairs {\A \T \T \A \C \G \G \C \a \T \t \A \c \G \g \C})

(defn reverse-complement
  "Returns the reverse complement.  Reverses the string and maps pairs on it."
  [s]
  (apply str (map pairs (reverse (vec s)))))
