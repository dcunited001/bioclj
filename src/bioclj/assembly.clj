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

(defn construct-debruijin-graph [k kmers]
  (let [k-1 (dec k)]
    (reduce
      (fn [graph kmer]
        (let [lmer (subs kmer 0 k-1)
              rmer (subs kmer 1 k)
              this-node (get graph lmer {})
              this-rmer-count (inc (get this-node rmer 0))]
          (assoc graph lmer (assoc this-node rmer this-rmer-count))))
      {}
      kmers)))

(defn debruijin-graph-to-str [dbg]
  (map
    (fn [lmer]
      (let [lmer-map (get dbg lmer)
            rmers (sort (keys lmer-map))
            concat-rmers (flatten (map #(repeat (get lmer-map %1) %1) rmers))]
        (str lmer " -> " (clojure.string/join "," concat-rmers))))
    (sort (keys dbg))))

(defn parse-graph [gstr]
  (let [node-strings (clojure.string/split gstr #"\s*\n\s*")
        nodes (map #(clojure.string/split %1 #"\s*->\s*") node-strings)]
    (reduce #(assoc %1 (first %2)
                    (frequencies (clojure.string/split (second %2) #",")))
            {} nodes)))

(defn rotate-path [v n]
  (let [i (- (count v) n)]
    (into [] (concat (subvec v i (count v)) (subvec v 0 i)))))

(defn test-eulerian-path
  "tests the solution for a eulerian cycle"
  [gr p]
  (let [start (first p)
        next (second p)
        these-edges (get gr start)
        next-val (get these-edges next)]
    (if (nil? next-val)
      (prn (str "WRONG!! not a valid path from " start " => " next))
      (let [new-edges (if (= 1 next-val)
                        (dissoc these-edges next)
                        (assoc these-edges next (dec next-val)))
            remaining-graph (if (empty? new-edges)
                              (dissoc gr start)
                              (assoc gr start new-edges))]
        (if (= (count p) 2)
          (do (prn (str "VALID!! (so far)"))
              true)
          (if (empty? remaining-graph)
            (do (prn (str "VALID!!") p remaining-graph gr) true)
            (recur remaining-graph (rest p))))))))

(defn eulerian-path-to-str [p]
  (clojure.string/join "->" p))

(defn solve-eulerian-graph-linear-time
  ([start-graph]
    (let [[selected-start next-nodes] (first start-graph)]
      (solve-eulerian-graph-linear-time [selected-start] start-graph [0])))

  ;; i'm assuming that the path through the graph should never cycle completely once
  ;;  - in other words, in the path chosen through the graph, every node shifted before the original starting node will always be closed
  ([start-cycle start-graph open-node-offset-stack]
    (let [selected-node (last start-cycle)
          next-nodes (get start-graph selected-node)
          [next-node node-val] (first next-nodes)
          new-paths (if (= 1 node-val)
                      (dissoc next-nodes next-node)
                      (assoc next-nodes next-node (dec node-val)))
          this-node-closed (empty? new-paths)

          [remaining-graph last-open-node-offset]
          (if this-node-closed
            [(dissoc start-graph selected-node) (inc (last open-node-offset-stack))]
            [(assoc start-graph selected-node new-paths) 1])
          next-node-closed (empty? (get remaining-graph next-node))
          finished? (empty? remaining-graph)

          [next-cycle next-offset-stack]
          (if (and (not finished?) next-node-closed)
            [(rotate-path start-cycle last-open-node-offset)
             (rotate-path open-node-offset-stack last-open-node-offset)]
            [(conj start-cycle next-node)
             (conj open-node-offset-stack last-open-node-offset)])

          next-graph (if (and (not finished?) next-node-closed)
                       (let [new-last-node (last next-cycle)
                             new-first-node (first next-cycle)
                             new-last-edges (get remaining-graph new-last-node {})
                             readd-cycle-edge (assoc new-last-edges new-first-node (get new-last-edges new-first-node 1))]
                         (assoc remaining-graph new-last-node readd-cycle-edge))
                       remaining-graph)]

      (if (empty? next-graph)
        next-cycle
        (recur next-cycle
               next-graph
               next-offset-stack)))))

(defn edge-totals-for-graph [gr]
  (reduce
    (fn [totals node]
      (let [edges (get gr node)
            [ins outs] (reduce
                         (fn [t [e n]]
                           [(assoc (first t) e n) (+ (second t) n)])
                         [{} 0]
                         edges)]
        {:ins  (merge-with + (get totals :ins) ins)
         :outs (merge-with + (get totals :outs) {node outs})}))
    {:ins  {}
     :outs {}}
    (keys gr)))

(defn find-missing-edge-for-eulerian-cycle [gr]
  (let [edge-totals (edge-totals-for-graph gr)
        in-totals (:ins edge-totals)
        out-totals (:outs edge-totals)]
    (reduce
      (fn [missing-edge key]
        (let [ins (get in-totals key)
              ins (if (nil? ins) 0 ins)
              outs (get out-totals key)
              outs (if (nil? outs) 0 outs)]
          (if (not= ins outs)
            (if (= ins (inc outs))
              [key (second missing-edge)]
              (if (= (inc ins) outs)
                [(first missing-edge) key]
                ;; this actually means invalid data
                missing-edge))
            missing-edge)))
      [nil nil]
      (keys gr))))

(defn solve-eulerian-path [start-graph]
  (let [[node edge] (find-missing-edge-for-eulerian-cycle start-graph)
        edges-at-node (get start-graph node)
        new-edges (if (empty? edges-at-node)
                    {node edge}
                    (if (nil? (get edges-at-node edge))
                      (assoc edges-at-node edge 1)
                      (assoc edges-at-node edge (inc (get edges-at-node edge)))))
        new-graph (assoc start-graph node new-edges)
        cycle (solve-eulerian-graph-linear-time [edge] new-graph [0])
        split-index
        ((fn search-cycle [c i]
           (if (> i (count c))
             -1
             (if (and (= node (get c (- (count c) (inc i))))
                      (= edge (get c (- (count c) i))))
               i
               (recur c (inc i)))))
          cycle 0)]
    (if (= split-index (dec (count cycle)))
      (rest cycle)
      (rotate-path (subvec cycle 1 (count cycle)) split-index))))

