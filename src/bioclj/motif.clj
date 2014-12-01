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


(defn median-string-max-distance
  [k dna]
  (* (count dna) k))

(defn median-string-distance
  [dna-seqs pattern]
  (let [k (:k (first dna-seqs))]
    (r/fold
      4
      (fn foldsum ([] 0) ([x y] (+ x y)))
      (fn sum-minimums [total dna-seq]
        (+ total
           (reduce
             #(let [dist (hamming-64b k pattern %2)]
               (if (< dist %1)
                 dist
                 %1))
             Double/POSITIVE_INFINITY
             (:b64 dna-seq))))
      dna-seqs)))


;; instead of using the set of all kmers as a starting domain
;; may be possible to use the disjunction of seq-neighbors for all seqs in dna
;; as a starting domain instead
(defn median-string
  [k dna]
  (let [dna-seqs (if (string? (first dna))
                   (map (partial acgt-get-64b-kmers k) dna)
                   dna)]
    (reduce
      #(let [shifted-pattern (bit-shift-left %2 (- 64 (* 2 k)))
             dist (median-string-distance dna-seqs shifted-pattern)]
        (if (< dist (:min %1))
          {:median shifted-pattern :min dist}
          %1))
      {:median -1 :min Double/POSITIVE_INFINITY}
      (range 0 (maths/expt 4 k)))))

(defn motif-count
  [motifs]

  )

(defn motif-profile
  [motifs]

  )

(defn motif-profile-most-probable-kmer [dna k profile])
(defn motif-most-probable-to-kmers [dna mpk])

;;TODO: profile record and methods for operating on them
(defn motif-profile-consensus-kmer-generator
  [dna k profile]
  (reduce
    (fn [nucs pos]
      (conj nucs
            (reduce
              (fn [nuc i]
                (let [new-max (nth profile (+ (* k i) pos))]
                  (if (> new-max (:max nuc))
                    {:max new-max :nuc [i]}
                    (if (= new-max (:max nuc))
                      (assoc nuc :nuc (conj (:nuc nuc) i))
                      nuc))))
              {:max -1 :nuc []}
              (range 0 4))))
    []
    (range k)))

(defn motif-profile-top-consensus-kmer
  "converts the array of hashed returned above to a nucleotide"
  [dna mpk]
  ;;TODO: alter this to return the first nucleotides in dna
  (reduce
    #(let [m (nth mpk %2)
           nucleo (nth (acgt-64b-shifted-nucleotides %2) (first (:nuc m)))]
      (bit-or %1 nucleo))
    0
    (range (count mpk))))
