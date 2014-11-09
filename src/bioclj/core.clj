(ns bioclj.core
  (:use [criterium.core])
  (:require [clojure.math.combinatorics :as com]
            [clojure.core.reducers :as r]
            [clj-biosequence.core :as bio.core]
            [clj-biosequence.alphabet :as bio.alpha]
            ))

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

(def base-pairs {\A \T \T \A \C \G \G \C \a \T \t \A \c \G \g \C})

(def nucleotides '(\A \G \T \C))
(def rna-nucleotides '(\A \G \T \C))
(def integer-mass-table {\G 57 \A 71 \S 87 \P 97 \V 99 \T 101 \C 103 \I 113 \L 113 \N 114
                         \D 115 \K 128 \Q 128 \E 129 \M 131 \H 137 \F 147 \R 156 \Y 163 \W 186})

(defn reverse-complement
  "Returns the reverse complement.  Reverses the string and maps pairs on it."
  [s]
  (apply str (map base-pairs (reverse (vec s)))))

;; TODO: try converting the strings to integers, then compare bitwise
(defn hamming-distance
  "return the number of mismatches in two equally sized strings, returning distance d
  use dmax to ensure that hamming distance terminates early, dmax defaults to (count a)
  if dmax terminates early (inc dmax) is returned, this is assuming that result will be filtered out anyways"
  ;; user must ensure strings are equal sized before calling
  ([a b & {:keys [d dmax] :or {d 0 dmax (count a)}}]
   (if (or (not (string? a)) (not (string? b)))
     ;; TODO: fix dmax so you don't have to call with :dmax (inc d)
     (if (< d dmax)
       (if (not-empty a)
         (recur (rest a) (rest b)
                {:d    (+ d (or (and (= (take 1 a) (take 1 b)) 0) 1))
                 :dmax dmax})
         d)
       (inc d))                                             ;; returns dmax + 1, when dmax is exceeded
     (hamming-distance (seq a) (seq b) :d d :dmax dmax)
     )))

(defn format-char-list [char-list]
  (map (partial apply str) char-list))

;;TODO: figure out why neighborhood-recur isn't correct
(defn neighborhood-recur
  [s d]
  ;; seems to be no way to use recur here, not that it would make much of a difference in this alg
  (if (> d 0)
    (if (> (count s) 1)
      (let [sub (rest s)
            subhood (neighborhood-recur sub d)]
        (reduce
          (fn [thishood t]
            (if (<= (hamming-distance sub t :dmax (inc d)) d)
              (apply conj thishood (map #(conj t %1) nucleotides))
              (conj thishood (conj t (first sub)))))
          '()
          subhood))
      (map list nucleotides))
    s))

(defn neighborhood-iter
  [s d]
  ;; i actually think this will be faster, when using the appropriate data structures
  )

(defn sorted-freq-kmer
  [freqmap]
  (into (sorted-map-by #(compare [(get freqmap %2) %2]
                                 [(get freqmap %1) %1])) freqmap))

;; TODO: refactor to macros?
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

(defn kmer-op-with-near-misses
  [f foldmergef k d text]
  (let [idxstop (- (count text) k)]
    (r/fold
      (fn combine-indexes
        ([] {})
        ([x y] (merge-with foldmergef x y)))
      (fn fold-substrings [h i]
        (let [kmer (subs text i (+ i k))]
          (f h i kmer)))
      (range 0 (inc idxstop))
      )))

(defn kmer-frequency-with-near-misses
  "returns a hash where keys are kmers and values are the frequency of that word and it's close mismatches with hamming distance < d"
  [k d text]
  (sorted-freq-kmer
    (kmer-op-with-near-misses
      (fn [h i kmer]
        (merge-with + h
                    (r/fold
                      (fn combinef ([] {}) ([x y] (merge x y)))
                      (fn fold-neighbor-indices [nh n] (assoc nh n 1))
                      (format-char-list (neighborhood-recur kmer d)))))
      + k d text)))

(defn kmer-indices-with-near-misses
  [k d text]
  (kmer-op-with-near-misses
    (fn [h i kmer]
      (merge-with concat h
                  (r/fold
                    (fn combinef ([] {}) ([x y] (merge x y)))
                    (fn fold-neighbor-indices [nh n] (assoc nh n (list i)))
                    (format-char-list (neighborhood-recur kmer d)))))
    concat k d text))

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
                                   arr)))
                             (sorted-set)
                             ind))]
         (if (not-empty clumped) (assoc! %1 key clumped) %1))
       (transient {})
       kindex)
     )))

(defn hammify-domain
  "permute a domain and return matches with hamming distance less than d for strings of length k"
  ;([d k] (hammify-domain (map (partial apply str) (com/selections "ACGT" k)) d k))
  ([domain d]
   (prn (now))
   (persistent!
     (reduce
       (fn [ham-dom s1]
         (assoc! ham-dom s1 (filter #(<= (hamming-distance s1 %1 :dmax (inc d)) d) domain)))
       (transient {})
       domain)
     )))

(defn hammify-domain-parallel
  [domain d]
  (r/fold
    (fn combinef
      ([] {})
      ([x y] (merge x y)))
    (fn [ham-dom s1]
      (assoc ham-dom s1 (into [] (r/filter #(<= (hamming-distance s1 %1 :dmax (inc d)) d) domain)))
      ham-dom)
    domain))

(defn hammify-domain-nonoccuring
  "takes the strings that compose the domain of a kmer-indices result (the keys)
  and contrasts them with the combinated domain of all possible kmers of length k, (subtracts them)
  then filters the remaining kmers in k by their hamming distances.

  this is such an inefficient operation that i can't believe it's required for our assignment
  and yet, i can't think of how you can identify kmers that are missing in your set,
  without starting with a complete set of all of them

  nevermind. by taking the domain of the result of kmer-hammers (DKH), i can compare that against the combinated-domain (CD),
  but that's still massive: |DKH| x |CD|. might be better to take DKH and generate all the possible strings from it, since |CD|=65K for k=8.

  by using memoized hamming-distances, i can further reduce time. of course, for k=8, that's 4GB of RAM, just to store the keys"
  []
  ;(into [] (r/filter #(<= (hamming-distance ))) (combinate-domain "ACGT" k)
  )

(defn kmer-hammers
  "get the kmer-indices of length k, then com/selections the keys to get the combos of strings of length k,
  for each string pair calc hamming distance and add it's match indexes to the each strings' total,
  if that hamming-distance is less than d, returning a list of kmers of size k and each kmer's hamming match occurances

  .. hmmm is this most efficient? also, the substring we're comparing against may not actually appear as a substring in the text..."
  ;; ugh .. rewrite to inspect each substring, then expand that substring into its neighborhood, then add to the list of indexes
  [k d text]
  (let [ki (kmer-indices k text)
        kmers (keys ki)
        ;kmers (format-char-list (r/mapcat #(neighborhood-recur %1 d) (keys ki))) ;; prohibitively expensive
        ham-dom (hammify-domain-parallel kmers d)]
    (persistent!
      (reduce
        (fn [hammers k]
          (assoc! hammers k (r/mapcat ki (ham-dom k))))
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

;; TODO: find the most frequent k-mers with up to d mismatches and includeing reverse complements...
;; TODO: rewrite everything to use r/fold that doesn't mutate state

(def codon-table (bio.alpha/codon-tables 1))

(defn transcribe-rna
  [rna &
   {:keys [aminos]
    :or   {aminos ""}}]
  (if (empty? rna)
    aminos
    (recur (nthrest rna 3) {:aminos (str aminos (bio.alpha/codon->aa (take 3 rna) codon-table))})))

(defn to-char-array [abc] (if (string? abc)
                            (map char (.getBytes abc))
                            abc))

;; convert to char arrays only
(defmulti rna-encodes-peptide? (fn [rna peptide] (map string? [rna peptide])))
(defmethod rna-encodes-peptide? :default [rna peptide]
  (rna-encodes-peptide? (to-char-array rna) (to-char-array peptide)))
(defmethod rna-encodes-peptide? [false false] [rna peptide]
  [rna peptide]
  (or (empty? rna)
      (and (= (bio.alpha/codon->aa (take 3 rna) codon-table)
              (first peptide))
           (recur (nthrest rna 3) (nthrest peptide 1)))))

;; TODO: add behavior for stop codon
(defn find-encoded-peptides
  [rna peptide &
   {:keys [find-rev reverse]
    :or   {find-rev true reverse false}}]
  (if find-rev
    (into {} (filter
               (comp not empty? val)
               (merge-with
                 concat
                 (find-encoded-peptides rna peptide :find-rev false)
                 (find-encoded-peptides (reverse-complement rna) peptide :find-rev false :reverse true))))
    (kmer-op #(if (rna-encodes-peptide? %3 peptide)
               ;; report back negative indices for the reverse complement
               (conj (%1 %3 []) (if reverse (* -1 %2) %2))
               (%1 %3 []))
             (* 3 (count peptide))
             rna)))

(defn format-encoded-peptides
  [peptides-map]
  (mapcat
    (fn [[strand idx]]
      (map
        #(if (> %1 0) strand
                      (apply str (reverse-complement strand)))
        idx))
    peptides-map))


