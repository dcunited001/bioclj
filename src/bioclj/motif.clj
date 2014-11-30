(ns bioclj.motif
  (:use [bioclj.string])
  (:require [clojure.core.reducers :as r]
            [clj-biosequence.core :as bio.core]
            [clj-biosequence.alphabet :as bio.alpha]
            [clojure.math.numeric-tower :as maths]))

(defn kmer-occurances
  "outputs the average number of occurrences of a random kmer within n random {ACTG} strings of L length"
  [k n L]
  (let [num-positions (inc (- L k))]
    (/ (* n num-positions 1.0) (maths/expt 4 k))))

(defn seq-neighbors
  "unions all the neighbors from each kmer in seq"
  [base-hood seq]
  (reduce
    (fn [neighbors kmer]
      (apply (partial conj neighbors) (:b64 (neighborhood-transform base-hood kmer))))
    #{}
    (:b64 seq)))

(defn motif-enumeration
  [k d seqs]

  ;1) init the first seq with a sorted set of all the neighbors of each kmer
  ;2) run through the rest of the seqs
  ;3) add each kmer's neighbors to a sorted set, if that kmer hasn't been run
  ;4) either get the union of the set, or iterate over the last set and delete values that don't exist in the new set
  ;5) return the composite set

  (let [base-hood (base-neighborhood-acgt k d)

        seqs-kmers (if (string? (first seqs))
                     (map (partial acgt-get-64b-kmers k) seqs)
                     seqs)]
    (reduce
      (fn [motifs seq-kmers]
        (let [seq-nbors (seq-neighbors base-hood seq-kmers)]
          (clojure.set/intersection motifs seq-nbors)))
      (seq-neighbors base-hood (first seqs-kmers))
      seqs-kmers)))

;; additional interesting notes for optimization:
;; if hdist(kmer1, kmer2) > 2d
;;   then neighborhood(kmer1) ^ neighborhood(kmer2) = null set
;; if hdist(kmer1, kmer2) = 2d
;;   then neighborhood(kmer1) ^ neighborhood(kmer2) = {predictable set?} with size d^2?
;;   the set may be somewhat predictable using neighborhood transforms? and ... ?
;; if there are multiple calls to retreive neighborhoods for the same kmer,
;;   this could either be memoized
;;   or could somehow utilize idempotency to minimize function calls
;; for example, if you generate motif on 10 sequences
;;   and those sequences each contain the same exact kmer at least once,
;;   then the function call to this kmer's neighborhood can safely be removed from the calculation for each seq
;;   then or'd into the result later on
;; also, for each kmer in each sequence, their neighborhoods are interrelated since k1 and k2 overlap
;;   but this probably isn't useful in simplifying the problem
;; and neighborhood for AAAA,CCCC is similar because hdist(AAAA,CCCC) = k and because AAAA - CCCC = 1111
;;   neighborhoods for CCGG,GGTT are similar in the same way because CCGG - GGTT = 1111
;; the intersections of neighborhoods AACC,AAGG is similar to CCGG,CCTT because (AACC-AAGG) = (CCGG-CCTT) = 0011
;;   but unfortunately, this problem involves the conjunction of neighborhoods
;;   there might be more to this that may be useful with optimizations
;; and there's probably some interesting stuff to pull in from functional analysis here
;;   since hdist is a distance function, though this isn't a metric space .. AFAIK
;; and there may be a way to utilize the change in neighborhood set from the base neighborhood, a ∂-N(kmer,d) function,
;;   to simplify the disjunction of neighborhoods

(defn score-motif
  [motifs]

  )

(defn count-motif
  [motifs]

  )

(defn profile-motif
  [motifs]

  )
