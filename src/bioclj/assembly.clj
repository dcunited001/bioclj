(ns bioclj.assembly
  (:use [bioclj.string])
  (:require [clojure.core.reducers :as r]
            [clj-biosequence.core :as bio.core]
            [clj-biosequence.alphabet :as bio.alpha]
            [clojure.math.numeric-tower :as maths]))

(defn composition-64b [k dna]
  (sort (acgt-get-64b-kmers k dna)))

(defn composition [k dna]
  (sort
    (let [idxstop (- (count dna) k)]
      (map
        #(subs dna %1 (+ %1 k))
        (range (inc idxstop))))))

(defn string-spelled-by-a-genome-path [k kmers]
  (reduce
    #(conj %1 (last %2))
    (mapv char (first kmers))
    (rest kmers)))

;; when constructing the graph in with binary kmers,
;; when you can assume that the kmers are only shifted by one,
;; - then you can leftshift kmer-a, then && it with kmer-b
;; - and you should get 11110000, all ones to the left
;; - and therefore a manageably finite number of possibilities
;; when you have to compare kmers can be shifted by varible lengths
;; - you can preprocess kmer-a to produce kmer-a1, kmer-a2, kmer-a3
;; - then simultaneously compare && all kmer-a's with kmer-b
;; - and the same one's to the left will emerge
;; - obviously, this should be restricted to some minimum overlap
;; when the kmers themselves are variable length
;; - this isn't really that big of a problem when constructing the initial graph
;; these additional variabilities complicate the rest of the problem
;; - however, in some ways, they simplify some of the graph manipulations later on

(defn overlap-graph [kmers]
  (let [k (count (first kmers))
        k-1 (dec k)]
    (reduce
      (fn [graph left-mer]
        (let [adj-nodes (reduce
                          (fn [matches right-mer]
                            (if (= (subs left-mer 1 k) (subs right-mer 0 k-1))
                              (conj matches right-mer)
                              matches))
                          []
                          kmers)]
          (if ((comp not empty?) adj-nodes)
            (conj graph [left-mer adj-nodes])
            graph)))
      []
      kmers)))

(defn adjacency-graph-to-str [og]
  (map
    (fn [node] (str (first node) " -> " (first (second node))))
    og))

(defn debruijin-graph-to-str [dbg]
  (map
    (fn [node] (str (first node) " -> " (clojure.string/join "," (second node))))
    dbg))

(defn construct-debruijin-graph [k dna]
  (let [kmers (composition k dna)
        k-1 (dec k)]
    (reduce
      (fn [graph kmer]
        (let [lmer (subs kmer 0 k-1)
              rmer (subs kmer 1 k)
              this-node (get graph lmer {})
              this-rmer-count (inc (get this-node rmer 0))]
          (assoc graph lmer (assoc this-node rmer this-rmer-count))))
      {}
      kmers)))

