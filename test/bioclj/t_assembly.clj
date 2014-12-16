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
             kmers (composition k dna)
             ans {"AAG" {"AGA" 2}
                  "AGA" {"GAT" 1}
                  "ATT" {"TTC" 1}
                  "CTA" {"TAA" 1}
                  "CTC" {"TCT" 1}
                  "GAT" {"ATT" 1}
                  "TAA" {"AAG" 1}
                  "TCT" {"CTA" 1 "CTC" 1}
                  "TTC" {"TCT" 1}}
             dbg (construct-debruijin-graph k kmers)]
         dbg => ans))

(facts "parse-graph"
       (fact "works for unbalanced graphs"
             (let [gstr "0 -> 3
                   1 -> 0
                   2 -> 1,6"
                   ans {"0" { "3" 1 }
                        "1" { "0" 1 }
                        "2" { "1" 1 "6" 1 }}]
               (parse-graph gstr) => ans))
       (fact "works for balanced graphs"
             (let [gstr "0 -> 0,1,1
                   1 -> 2
                   2 -> 3
                   3 -> 0,4,0
                   4 -> 3"
                   ans {"0" {"0" 1 "1" 2}
                        "1" {"2" 1}
                        "2" {"3" 1}
                        "3" {"0" 2 "4" 1}
                        "4" {"3" 1}}]
               (parse-graph gstr) => ans)))
