(ns bioclj.string
  (:require [clojure.core.reducers :as r]
            [clj-biosequence.core :as bio.core]
            [clj-biosequence.alphabet :as bio.alpha]
            [clojure.pprint :refer (cl-format) :as pp]
            [clojure.math.numeric-tower :as maths]))

;; TODO: rename to bioclj.string?
;; TODO: look into process builder for implementing library with GPU calls.
;; TODO: check out kibit

(def acgt-2 {\A 0 \a 0 97 0 65 0
             \C 1 \c 1 99 1 67 1
             \G 2 \g 2 103 2 71 2
             \T 3 \t 3 116 3 84 3
             \U 3 \u 3 117 3 85 3})
(def b2-to-acgt {0 \A
                 1 \C
                 2 \G
                 3 \T})

(def two-multiples (map #(* 2 %1) (range 32)))
(def acgt-bitshifts (map #(long (maths/expt 2 %1)) two-multiples))
(def acgt-index-bitmasks (map #(bit-flip (bit-flip 0 %1) (inc %1)) (reverse two-multiples)))

;; 000000000
;; 110000000
;; 001100000
;; 000011000
;; ....

(def acgt-lbitmasks (reduce #(conj %1 (bit-or (last %1)
                                              (bit-flip (bit-flip 0 %2) (inc %2))))
                            [(bit-flip (bit-flip 0 62) 63)]
                            (reverse (take 31 two-multiples))))

;; 00000000
;; 00000011
;; 00001100
;; 00110000
;; ....

(def acgt-rbitmasks (reduce #(conj %1 (bit-or (last %1)
                                              (bit-flip (bit-flip 0 %2) (inc %2))))
                            [(bit-flip (bit-flip 0 0) 1)]
                            (nthrest two-multiples 1)))

(defn split-kmers [kmer-str] (clojure.string/split kmer-str #" "))

(defn complement-64b [k kmer]
  (bit-and (bit-not kmer) (acgt-lbitmasks (dec k))))

;; k=4 0000111100000000 => 0000111100000000
;; k=4 0101101100000000 => 0010010100000000
;; k=5 0100110011000000 => 0011001101000000

;; TODO: implement reverse-complement-64b
;(defn reverse-complement-64b [k kmer])

;; contains the values acgt, shifted for each index
(def acgt-64b-shifted-nucleotides
  (vec (map
         (fn [pos]
           (vec (map #(bit-shift-left %1 (- 62 (* 2 pos))) [0 1 2 3])))
         (range 32))))

;given a length and a vector of 64b longs, creates a record that represents a ACGT sequence
(defprotocol VariableBitOps
  (sub64b [k i] "return the i'th element of k length in a 64b encoded sequence as a 64b int"))
;;TODO: include another method that outputs every sub for length k as 64b ints
;;TODO: implement another hamming-b64 for comparing 2 Acgt64's and automatically identifying differing regions

(defrecord Acgt64Contiguous [L orig-seq b64]
  ;;TODO: add the original string to the fields to avoid unnecessary conversions
  VariableBitOps
  (sub64b [k i] []))

(defn acgt-str-to-64b
  "converts up to 32 bases into a 64-bit long"
  [bases]
  (reduce
    #(bit-or %1 (bit-shift-left (acgt-2 (nth bases %2)) (- 62 (* 2 %2))))
    0
    (range (count bases))))

(defn acgt-to-64b
  "returns a vector of longs that represent the encoded ACTG string, with the original number of bases prepended"
  ;; appends the original count of bases to the front of the vector
  ([s] (apply (partial conj [(count s)]) (acgt-to-64b [] s)))
  ([longs s]
    (let [b64 (acgt-str-to-64b (take 32 s))
          rest-bases (nthrest s 32)]
      (if (not-empty rest-bases)
        (recur (conj longs b64) rest-bases)
        (conj longs b64)))))

(defn acgt-64b-to-str
  [k long]
  (reduce #(conj %1 (b2-to-acgt (bit-and 3 (bit-shift-right long (- 62 %2)))))
          []
          (take k two-multiples)))

(defn acgt-from-64b
  "takes a vector of longs and returns the original ACGT string"
  ([k longs] (acgt-from-64b k longs []))
  ([k longs chr-list]

    (if (> k 0)
      (let [this-k (or (and (> k 32) 32) k)
            chr32 (acgt-64b-to-str this-k (first longs))]
        (recur (- k 32) (rest longs) (apply (partial conj chr-list) chr32)))
      chr-list)))

(defn acgt-str-from-64b
  [k longs] (apply str (acgt-from-64b k longs '())))

(defn acgt-truncate-to-k
  "truncates a 64b long to a using the acgt-lbitmasks"
  [k a]
  (bit-and a (nth acgt-lbitmasks (dec k))))

(defrecord Acgt64Kmers [k orig-seq b64]
  VariableBitOps
  (sub64b [k i] []))

(defn acgt-get-64b-kmers
  [k dna]
  (Acgt64Kmers. k
                dna
                (let [idxstop (- (count dna) k)]
                  (reduce
                    #(let [next-nucleotide (acgt-2 (get dna (+ %2 (dec k))))
                           shifted-nucleotide (nth (acgt-64b-shifted-nucleotides (dec k)) next-nucleotide)
                           kmer (bit-or (bit-shift-left (last %1) 2) shifted-nucleotide)]
                      (conj %1 kmer))
                    [(acgt-str-to-64b (subs dna 0 k))]
                    (range 1 (inc idxstop))))))

;; for 2bit inputs a & b
;; 1) xor a,b => x
;; 2) (dup x)*2 => y
;; 3) init array of 10101010 bytes => z
;; 4) (y & z) | x => hamming distance, for 2 bits
;; 5) run hammingweight to get distance (00111100 => 4 / 2 => 2)
(def hamming-64b-b2-magic-xor (reduce #(bit-or %1 (bit-flip 0 (inc %2))) 0 two-multiples))
(defn hamming-64b
  "returns the hamming-distance between two b64 encoded nucleotide sequences"
  ([a b]
    (let [x (bit-xor a b)]
      (. Long bitCount
         (bit-and (bit-or x
                          (bit-and (bit-shift-left x 1)
                                   hamming-64b-b2-magic-xor))
                  hamming-64b-b2-magic-xor))))
  ([k a b]
    (hamming-64b (acgt-truncate-to-k k a) (acgt-truncate-to-k k b))))

;;TODO: implement hamming-weight with bit count
;; either use: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
;; or use an algorithm with an 8b or 16b lookup table
;; or just use (. Long bitCount (long 1)) to call to Java's default implementation (may require -XX:+UsePopCountInstruction jvmOpt)

;; TODO: b64 Levenschtien Distance? ........
;; .. i think levenschtein is more useful with amino seqs anyways,
;; nucleotide deletions & additions cause terrible things to happen

(defn format64
  [long]
  (if (= long -9223372036854775808)
    ;; on 64b minimum, cl-format breaks
    (str "-" (pp/cl-format nil "~63,'0',B" 0))
    (pp/cl-format nil "~64,'0',B" long)))

;; TODO: hamming distance for 6bit amino's
;; use similar process to above, but recursive & repeat 4 times
;; end up with annoying 6b array 000000111111111111000000 => 0110
;; .. is this even worthwhile to implement?? .. 10 amino's per 64b int plus 4-bits .. hmmm
;; would the increase in performance from 8b => 6b offset the overhead?
;; how to include different alphabets?

(defprotocol NeighborOps
  (neighbors-to-str [_] "converts the list of longs into strings")
  (neighbors-one-str [_] "formats output to be submitted"))

(defrecord Acgt64bNeighborhood [k d orig-seq b64]
  NeighborOps
  (neighbors-to-str [this] (->> b64
                                (map (partial acgt-64b-to-str k))
                                (map (partial apply str))))
  (neighbors-one-str [this] (println (clojure.string/join "\n" (neighbors-to-str this)))))

;; 11......
;; 0011....
;; 1111....
;; 000011..
;; 001111..
;; 111111..
;; generates a map of bitmasks to use for the neighborhood algorithm
(def neighborhood-64b-masks
  (reduce
    (fn [h i]
      (let [end (- 62 (* i 2))]
        (assoc h i
               (reduce
                 (fn [a k]
                   (let [start (+ end k 2)]
                     (conj a (bit-flip (bit-flip (last a) start) (inc start)))))
                 [(bit-flip (bit-flip 0 end) (inc end))]
                 (take i two-multiples)))))
    {}
    (range 32)))

;;might be possible to process this 8+ bits at a time, truncate the results & remove dupes.. hmmmmm
;; if so, then that algorithm might be worth parallelizing with GPU
(defn neighborhood-acgt-64b
  "recursively generates the actg neighborhood of a 64b encoded actg string of up to 32 chars"
  ([k d lng]
    (if (<= d 0)
      [lng]
      (Acgt64bNeighborhood. k d lng (neighborhood-acgt-64b k d lng k []))))
  ([k d lng pos longs]
    (let [nucleotides (nth acgt-64b-shifted-nucleotides (dec pos))]
      (if (= pos 1)
        nucleotides
        (let [sub (bit-and lng (nth (neighborhood-64b-masks (dec pos)) (dec pos)))
              subhood (neighborhood-acgt-64b k d sub (dec pos) longs)]
          (vec (r/fold 4
                       (fn combinef ([] []) ([x y] (concat x y)))
                       (fn [a n]
                         (apply (partial conj a)
                                (reduce
                                  (fn [aa s]
                                    (let [sn (bit-xor s n)
                                          hd (hamming-64b k sub sn)]
                                      (if (<= hd d)
                                        (conj aa sn)
                                        aa)))
                                  []
                                  subhood)))
                       nucleotides)))))))

(defn base-neighborhood-acgt
  [k d]
  ;; idn't that zero just purdy
  (neighborhood-acgt-64b k d 0))

(def neighbor-transform-magic-numbers
  {:left  hamming-64b-b2-magic-xor
   :right (bit-not hamming-64b-b2-magic-xor)})

(defn neighbor-transform
  "given a kmer and the base-neighborhood neighbor (nbor), performs a kind of 2-bit addition to transform that neighbor to the neighborhood for kmer"
  [kmer nbor]
  (let [el-magic (:left neighbor-transform-magic-numbers)
        r-magic (:right neighbor-transform-magic-numbers)
        x (unchecked-add (bit-and kmer r-magic)
             (bit-and nbor r-magic))
        y (unchecked-add (bit-and kmer el-magic)
             (bit-and nbor el-magic))
        z (bit-xor x y)]
    (bit-or (bit-and x r-magic)
            (bit-and z el-magic))))

(defn neighborhood-transform
  "given a base neighborhood for k & d, transforms to the neighborhood for kmer.  GPU Acceleratable :)"
  [base-neighborhood kmer-lng]
  (Acgt64bNeighborhood.
    (:k base-neighborhood)
    (:d base-neighborhood)
    kmer-lng
    (map (partial neighbor-transform kmer-lng) (:b64 base-neighborhood))))