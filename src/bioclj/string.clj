(ns bioclj.string
  (:use [gloss core io])
  (:require [clojure.core.reducers :as r]
            [clj-biosequence.core :as bio.core]
            [clj-biosequence.alphabet :as bio.alpha]))

;; TODO: rename to bioclj.string?
;; TODO: look into process builder for implementing library with GPU calls.
;; TODO: check out kibit

(def acgt-2 {97 0
             65 0
             99 1
             67 1
             103 2
             71 2
             116 3
             84 3})

(defcodec acgt-codec (repeated :int64-be :prefix :none))

(defn acgt-to-64b
  [s])

(defn acgt-from-64b
  [s])

(defn byter-2
  "converts a string of ACGT to a 2 bit byte array"
  [s]
  (bytes (byte-array (map (comp byte acgt-2 int) s))))

;; for 2bit inputs a & b
;; 1) xor a,b => x
;; 2) (dup x)*2 => y
;; 3) init array of 10101010 bytes => z
;; 4) (y & z) | x => hamming distance, for 2 bits
;; 5) factor result to step down from 2 bit result to 1bit array (00111100 => 0110)

(defn byter-6
  "converts a string of aminos to a 6 bit byte array"
  [s])

;; use similar process to above, but recursive & repeat 4 times
;; end up with annoying 6b array 000000111111111111000000 => 0110
