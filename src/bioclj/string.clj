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

(def two-multiples (map #(* 2 %1) (range 32)))
(def acgt-bitshifts (map #(long (maths/expt 2 %1)) two-multiples))
(def acgt-lbitmasks (reduce #(conj %1 (bit-or (last %1)
                                              (bit-flip (bit-flip 0 %2) (inc %2))))
                            [(bit-flip (bit-flip 0 62) 63)]
                            (reverse (take 31 two-multiples))))
(def acgt-rbitmasks (reduce #(conj %1 (bit-or (last %1)
                                              (bit-flip (bit-flip 0 %2) (inc %2))))
                            [(bit-flip (bit-flip 0 0) 1)]
                            (nthrest two-multiples 1)))

(defn reduce-acgt-64b
  "converts up to 32 bases into a 64-bit long"
  [bases]
  ;(prn (map #(nth bases %1) (range (count bases))))
  (reduce
    #(bit-or %1 (bit-shift-left (acgt-2 (nth bases %2)) (- 62 (* 2 %2))))
    0
    (range (count bases))))

(defn acgt-to-64b
  "returns a vector of longs that represent the encoded ACTG string, with the original number of bases prepended"
  ;; appends the original count of bases to the front of the vector
  ([s] (apply (partial conj [(count s)]) (acgt-to-64b [] s)))
  ([longs s]
   (let [b64 (reduce-acgt-64b (take 32 s))
         rest-bases (nthrest s 32)]
     (if (not-empty rest-bases)
       (recur (conj longs b64) rest-bases)
       (conj longs b64)))))

;;TODO: implement hamming-weight with bit count
;; either use: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
;; or use an algorithm with an 8b or 16b lookup table
;; or just use (. Long bitCount (long 1)) to call to Java's default implementation (may require -XX:+UsePopCountInstruction jvmOpt)

(defn acgt-from-64b
  "takes a vector of longs and returns the original ACGT string"
  [s]
  )

(defn format64
  [long]
  (pp/cl-format nil "2r~64,'0',B" long))

(defn byter-2
  "converts a string of ACGT to a 2 bit byte array"
  [s]
  (bytes (byte-array (map (comp byte acgt-2 int) s))))

;; for 2bit inputs a & b
;; 1) xor a,b => x
;; 2) (dup x)*2 => y
;; 3) init array of 10101010 bytes => z
;; 4) (y & z) | x => hamming distance, for 2 bits
;; 5) run hammingweight to get distance (00111100 => 4 / 2 => 2)

(defn byter-6
  "converts a string of aminos to a 6 bit byte array"
  [s])

;; use similar process to above, but recursive & repeat 4 times
;; end up with annoying 6b array 000000111111111111000000 => 0110
