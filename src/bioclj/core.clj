(ns bioclj.core
  (:use [criterium.core])
  (:require [clojure.math.combinatorics :as com]
            [clojure.core.reducers :as r]))

(defn now [] (java.util.Date.))

(defn combinate-domain
  "Permutes a string or array of chars and returns an array of strings representing the domain"
  [chars n]
  ;; the algorithm this way is very inefficient.  hash is always most efficient.
  ;; plus, clojure already has the com/permutate function
  (com/selections chars n))

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

(defn kmer-hamming-count
  [kmer d text]
  (kmer-op #(if (<= (hamming-distance kmer %3 :dmax (inc d)) d)
             (inc (%1 %3 0)) (%1 %3 0))
           (count kmer)
           text))

(defn kmer-clumps
  "filters the kmer-indices down by number of occurrances, then passes to kmer-clumps [kindex L n],
  which further filters the kmer-indices passed in and removes indices that not occuring n times within a window of L distance"
  ([L n k text] (kmer-clumps L n (filter #(>= (count (second %)) n) (kmer-indices k text))))
  ([L n kindex]
   (persistent!
     (reduce
       #(let [key (first %2)
              ;; ind is the set of matching indices to inspect
              ind (second %2)
              clumped (vec (reduce-kv
                             (fn [arr i val]
                               ;;TODO: fix 9999999! lol
                               (let [clump-end (get ind (+ i (dec n)) 9999999)]
                                 (if (<= (- clump-end val) L)
                                   ;; ugh there's still got to be a better way to do this!
                                   (apply (partial conj arr)
                                          (subvec ind (.indexOf ind val)
                                                  (inc (.indexOf ind clump-end))))
                                   arr)
                                 )
                               )
                             (sorted-set)
                             ind))]
         (if (not-empty clumped) (assoc! %1 key clumped) %1))
       (transient {})
       kindex)
     )
   )
  )

(defn sorted-freq-kmer
  [freqmap]
  (into (sorted-map-by #(compare [(get freqmap %2) %2]
                                 [(get freqmap %1) %1])) freqmap))

(def base-pairs {\A \T \T \A \C \G \G \C \a \T \t \A \c \G \g \C})

(defn reverse-complement
  "Returns the reverse complement.  Reverses the string and maps pairs on it."
  [s]
  (apply str (map base-pairs (reverse (vec s)))))

(defn hamming-distance
  "return the number of mismatches in two equally sized strings, returning distance d
  use dmax to ensure that hamming distance terminates early, dmax defaults to (count a)
  if dmax terminates early (inc dmax) is returned, this is assuming that result will be filtered out anyways"
  ;; user must ensure strings are equal sized before calling
  ([a b & {:keys [d dmax] :or {d 0 dmax (count a)}}]
   (if (not (string? a))
     (if (< d dmax)
       (if (not-empty a)
         (recur (rest a) (rest b)
                {:d    (+ d (or (and (= (take 1 a) (take 1 b)) 0) 1))
                 :dmax dmax})
         d)
       (inc d))                                             ;; returns dmax + 1, when dmax is exceeded
     (hamming-distance (seq a) (seq b) :d d :dmax dmax)
     )))

(defn hammify-domain
  "permute a domain and return matches with hamming distance less than d for strings of length k"
  ;([d k] (hammify-domain (map (partial apply str) (com/selections "ACGT" k)) d k))
  ([domain d]
     (r/fold
       (fn combinef
         ([] {})
         ([x y] (merge x y)))
       (fn [ham-dom s1]
         (assoc ham-dom s1 (into [] (r/filter #(<= (hamming-distance s1 %1 :dmax d) d) domain)))
         (prn s1)
         ham-dom)
       domain)
   ))

(defn kmer-hammers
  "get the kmer-indices of length k, then com/selections the keys to get the combos of strings of length k,
  for each string pair calc hamming distance and add it's match indexes to the each strings' total,
  if that hamming-distance is less than d, returning a list of kmers of size k and each kmer's hamming match occurances

  .. hmmm is this most efficient? also, the substring we're comparing against may not actually appear as a substring in the text..."
  [k d text]
  (let [ki (kmer-indices k text)
        kmers (keys ki)
        ham-dom (hammify-domain kmers d)]
    (persistent!
      (reduce
        (fn [hammers k]
          (assoc! hammers k (mapcat ki (ham-dom k))))
        (transient {})
        kmers)
      )))

(defn kmer-hammer
  "kmer-hammers, but for one string, because i'm short on time, 15 mins to go"
  ([kmer d text]
   (let [kmers (kmer-indices (count kmer) text)]
     (filter #(<= (hamming-distance kmer (first %1)) d) kmers))))


(defn kmer-frequency-with-mismatches
  "returns a hash where keys are kmers and values are the frequency of that word and it's close mismatches with hamming distance < d"
  [k d text]
  (sorted-freq-kmer
    (reduce #(assoc %1 (first %2) (count (second %2)))
            {}
            (kmer-hammers k d text))
    ))

;; find the most frequent k-mers with up to d mismatches
;; find the most frequent k-mers with up to d mismatches and includeing reverse complements...
