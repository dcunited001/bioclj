(ns bioclj.t-assembly
  (:use [midje.sweet]
        [bioclj.assembly]
        [bioclj.string]))

(facts "overlap-graph"
       (let [ans [["GCATG" ["CATGC"]]
                  ["CATGC" ["ATGCG"]]
                  ["AGGCA" ["GGCAT"]]
                  ["GGCAT" ["GCATG"]]]
             og (overlap-graph (split-kmers "ATGCG GCATG CATGC AGGCA GGCAT"))]
         og => ans))


(facts "construct-debruijin-graph"
       (let [k 4
             dna "AAGATTCTCTAAGA"
             ans {"AAG" {"AGA" 2}
                  "AGA" {"GAT" 1}
                  "ATT" {"TTC" 1}
                  "CTA" {"TAA" 1}
                  "CTC" {"TCT" 1}
                  "GAT" {"ATT" 1}
                  "TAA" {"AAG" 1}
                  "TCT" {"CTA" 1 "CTC" 1}
                  "TTC" {"TCT" 1}}
             dbg (construct-debruijin-graph k dna)]
         dbg => ans))

