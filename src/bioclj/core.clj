(ns bioclj.core
  (:require [clojure.math.combinatorics :as com]))

(defn combinate-domain
  "Permutes a string or array of chars and returns an array of strings representing the domain"
  [chars]
  ;; the algorithm this way is very inefficient.  hash is always most efficient.
  ;; plus, clojure already has the com/permutate function
  (com/selections chars (count chars)))

;; convert permuted domain to strings, sort strings by dictionary order,
;; use strings as keys to hash that can be looked up by integers
;; this will allow us to get stats on hamming-distance near-miss strings that weren't in the text

;; TODO: hmmm.. how to do this for a power set, including strings of multiple lengths?
;; i guess that's not too hard. indexes are in the set of 0 to ∑(4^k) - 1, from k = 1 to K
;; ... gee thanks Coursera Analytics Combonatorics class that I barely understand lulz

(defn domain-index
  "returns an index in a combinated domain, given a string.
  combinated domain is an ordered map that can look up by keys or values."
  [dom s])

(defn domain-val
  "returns an value, corresponding to the integer index of a combinated domain"
  [dom i])

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

(defn kmer-clumps
  [k L n text]
  ;; iterate over vec values in kmer-indices and increment if distance between n kmer's is >= L
  )

(defn sorted-freq-kmer
  [freqmap]
  (into (sorted-map-by #(compare [(get freqmap %2) %2]
                                 [(get freqmap %1) %1])) freqmap))

(def pairs {\A \T \T \A \C \G \G \C \a \T \t \A \c \G \g \C})

(defn reverse-complement
  "Returns the reverse complement.  Reverses the string and maps pairs on it."
  [s]
  (apply str (map pairs (reverse (vec s)))))

(defn hamming-distance
  "return the number of mismatches in two equally sized strings
  if not equally sized, return -1"
  [a b]
  (if (= (count a) (count b))
    (reduce
      #(if (not= (get a %2) (get b %2)) (inc %1) %1)
      0
      (range 0 (count a)))
    -1
    ))

(defn kmer-hammers
  "get the kmer-indices of length k, then permutate the keys to get the combos of strings of length k,
  for each string pair calc hamming distance and add it's match indexes to the each strings' total,
  if that hamming-distance is less than d, returning a list of kmers of size k and each kmer's hamming match occurances

  .. hmmm is this most efficient? also, the substring we're comparing against may not actually appear as a substring in the text...
  "
  [k d text]
  true)

;; find the most frequent k-mers with up to d mismatches
;; find the most frequent k-mers with up to d mismatches and includeing reverse complements...