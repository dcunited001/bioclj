(ns bioclj.core
   (:require [clojure.math.combinatorics :as com]))

(defn combinate-domain
  "Permutes a string or array of chars and returns an array of strings representing the domain"
  [chars]
  ;; the algorithm this way is very inefficient.  hash is always most efficient.
  ;; plus, clojure already has the com/permutate function
  (com/combinations chars (count chars)))

(defn count-str [str in])
(defn iterate-str [in fn])
