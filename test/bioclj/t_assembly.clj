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
                   ans {"0" {"3" 1}
                        "1" {"0" 1}
                        "2" {"1" 1 "6" 1}}]
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

(facts "solve-eulerian-graph-linear-time"
       (facts "solves a bad example, needing no shifts"
              (let [graph (parse-graph "0 -> 3
                                 1 -> 0
                                 2 -> 1,6
                                 3 -> 2
                                 4 -> 2
                                 5 -> 4
                                 6 -> 5,8
                                 7 -> 9
                                 8 -> 7
                                 9 -> 6")]
                (solve-eulerian-graph-linear-time graph) => ["9" "6" "5" "4" "2" "1" "0" "3" "2" "6" "8" "7" "9"]))

       (facts "solves a better example, needing a few shifts"
              (let [graph (parse-graph "0 -> 3
                                 1 -> 0,7
                                 2 -> 1,6
                                 3 -> 2
                                 4 -> 2
                                 5 -> 4
                                 6 -> 5,8
                                 7 -> 9,1
                                 8 -> 7
                                 9 -> 6")]
                (solve-eulerian-graph-linear-time graph) => ["7" "9" "6" "5" "4" "2" "1" "0" "3" "2" "6" "8" "7" "1" "7"])))

(facts "rotate-path"
       (let [v1 [1 2 3 4 5 6 7 8]]
         (rotate-path v1 2) => [7 8 1 2 3 4 5 6]
         (rotate-path v1 7) => [2 3 4 5 6 7 8 1]))

(facts "edge-totals-for-graph"
       (let [graph (parse-graph "0 -> 3
                                 1 -> 0,7
                                 2 -> 1,6
                                 3 -> 2
                                 4 -> 2
                                 5 -> 4
                                 6 -> 5,8
                                 7 -> 9,1
                                 8 -> 7
                                 9 -> 6")
             totals (edge-totals-for-graph graph)
             ins (:ins totals)
             outs (:outs totals)]
         (get ins "0") => (get outs "0")
         (get ins "7") => 2
         (get outs "7") => 2))

(facts "find-missing-edge-for-eulerian-cycle"
      (let [graph (parse-graph "0 -> 3
                                 1 -> 0,7
                                 2 -> 1,6
                                 3 -> 2
                                 4 -> 2
                                 5 -> 4
                                 6 -> 5,8
                                 7 -> 9
                                 8 -> 7
                                 9 -> 6")
            missing-edge (find-missing-edge-for-eulerian-cycle graph)]
        missing-edge => ["7" "1"])
      (let [graph (parse-graph (slurp "test/data/eulerian-path.txt"))
            missing-edge (find-missing-edge-for-eulerian-cycle graph)]
        missing-edge => ["1587" "1630"])
      (let [kmers (clojure.string/split "CTTA ACCA TACC GGCT GCTT TTAC" #" ")
            graph (construct-debruijin-graph 4 kmers)
            missing-edge (find-missing-edge-for-eulerian-cycle graph)]
        missing-edge => ["CCA" "GGC"]))

(facts "solve-eulerian-path"
       (let [graph (parse-graph "0 -> 3
                                 1 -> 0,7
                                 2 -> 1,6
                                 3 -> 2
                                 4 -> 2
                                 5 -> 4
                                 6 -> 5,8
                                 7 -> 9
                                 8 -> 7
                                 9 -> 6")
             ep (solve-eulerian-path graph)]
         ep => ["1" "0" "3" "2" "1" "7" "9" "6" "5" "4" "2" "6" "8" "7"])
       (let [graph (parse-graph (slurp "test/data/eulerian-path.txt"))
             ep (solve-eulerian-path graph)]
         (first ep) => "1630"
         (last ep) => "1587"))

(facts "string-reconstruction"
       (let [k 4
             kmers (clojure.string/split "CTTA ACCA TACC GGCT GCTT TTAC" #" ")]
         (string-reconstruction k kmers) => "GGCTTACCA"))



