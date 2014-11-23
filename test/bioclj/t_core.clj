(ns bioclj.t-core
  (:use midje.sweet)
  (:require [bioclj.core :as core]))

(facts permute-domain
       (fact "Permuting strings works"
             (core/combinate-domain "AC" 2) => '((\A \A) (\A \C) (\C \A) (\C \C)))

       (fact "Permuting arrays works"
             (core/combinate-domain [\A \B] 2) => '((\A \A) (\A \B) (\B \A) (\B \B))))

(facts kmer-frequency
       (fact "counting kmers works for a 1-mers"
             (let [freq (core/kmer-frequency 1 "ABCDCD")]
               (freq "D") => 2
               (freq "C") => 2
               (freq "A") => 1))

       (fact "counting kmers works for a 2-mers"
             (let [freq (core/kmer-frequency 2 "ABCDCD")]
               (freq "DC") => 1
               (freq "CD") => 2
               (freq "AB") => 1)))

(facts sorted-freq-kmer
       (fact "sorting kmers ma
       s works for short kmers maps"
             (let [dataset "GTACTCAAGGCTAATAGTCGCTAATAGTCGCTAATAGTCGCCTCTGTCCGCTAATAGTCACAGATTCACAGATTCGCTAATAGTCGTACTCAAGACAGATTCGCCTCTGTCCGTACTCAAGGCTAATAGTCACAGATTCGCTAATAGTCGTACTCAAGTCGGAGTAAAGCCTCTGTCCGCTAATAGTCACAGATTCACAGATTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCGTACTCAAGTCGGAGTAAATCGGAGTAAAGCCTCTGTCCTCGGAGTAAAGTACTCAAGACAGATTCGCCTCTGTCCGCCTCTGTCCGTACTCAAGGCCTCTGTCCGCCTCTGTCCTCGGAGTAAAGTACTCAAGTCGGAGTAAAACAGATTCGTACTCAAGACAGATTCTCGGAGTAAAGCCTCTGTCCACAGATTCACAGATTCGTACTCAAGACAGATTCGCCTCTGTCCGCCTCTGTCCACAGATTCGCTAATAGTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCACAGATTCGTACTCAAGGTACTCAAGACAGATTCGTACTCAAGTCGGAGTAAAGCCTCTGTCCTCGGAGTAAAGTACTCAAGTCGGAGTAAAGCTAATAGTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCGTACTCAAGACAGATTCTCGGAGTAAAGTACTCAAGTCGGAGTAAAGTACTCAAGGCCTCTGTCCTCGGAGTAAATCGGAGTAAAGTACTCAAGGCCTCTGTCCGTACTCAAGGCTAATAGTCGCCTCTGTCCACAGATTCTCGGAGTAAAGCTAATAGTCGTACTCAAGACAGATTCGCTAATAGTCGCCTCTGTCCGCCTCTGTCC"
                   freqmap (core/sorted-freq-kmer {"CCGCCTCTGTCCGT" 8, "GAGTAAAACAGATT" 1, "TAAAGTACTCAAGT" 3, "ATAGTCACAGATTC" 4, "CACAGATTCGTACT" 2, "TAAAGCCTCTGTCC" 8})
                   results (core/sorted-freq-kmer (core/kmer-frequency 14 dataset))]

               (first freqmap) => ["TAAAGCCTCTGTCC" 8]
               (second freqmap) => ["CCGCCTCTGTCCGT" 8]
               (first results) => ["TCGGAGTAAAGCCT" 8]
               (second results) => ["TAAAGCCTCTGTCC" 8])))

(facts kmer-indices
       (fact "returns a hashmap with vectors containing the indices substrings are found"
             (let [freq (core/kmer-indices 2 "ABCDCD")]
               (freq "DC") => [3]
               (freq "CD") => [2 4]
               (freq "AB") => [0])))

(facts kmer-clumps
       (fact "returns a list of keys for values in a kmer-indices result which have n+ matches, within L chars"
             (let [freq ()])
             ))

(facts reverse-complement
       (fact "returns the reverse complement of a strand of nucleotides, regardless of whether it's in upper/lower case."
             (let [strand "TCGGAGTAAAGCCT"
                   strand2 "TcgGAgtAaaGcCT"]
               (core/reverse-complement strand) => "AGGCTTTACTCCGA"
               (core/reverse-complement strand2) => "AGGCTTTACTCCGA")))

(facts hamming-distance
       (fact "returns the number of mismatches in two equally sized strings"
             (let [s1 "ACTG"
                   s2 "AATG"
                   s3 "ACTGA"]
               (core/hamming-distance s1 s1) => 0
               (core/hamming-distance s1 s2) => 1
               (core/hamming-distance s1 s1 :d 0) => 0
               (core/hamming-distance s1 s2 :d 0) => 1
               (core/hamming-distance s1 s1 :d 0 :dmax 0) => 1)))

(facts neighborhood-recur
       (let [s1 '(\A \A)
             n1 (sort '("AC" "AT" "AG" "CA" "TA" "GA" "AA"))
             s2 '(\A \A \A)
             n2 (sort '("CAA" "TAA" "GAA" "AAA" "AGA" "ACA" "ATA" "AAG" "AAT" "AAC"))]
         (fact "generates a matching neighborhood of an ACGT string"
               (sort (map (partial apply str) (core/neighborhood-recur s1 1))) => n1
               (sort (map (partial apply str) (core/neighborhood-recur s2 1))) => n2)
         (fact "generates a neighborhood with the appropriate counts"
               (count (core/neighborhood-recur s1 1)) => 7
               (count (core/neighborhood-recur s1 2)) => 16
               (count (core/neighborhood-recur s2 1)) => 10
               (count (core/neighborhood-recur s2 2)) => 37
               (count (core/neighborhood-recur s2 3)) => 64)
         (fact "neighborhoods include the input string itself"
               (.contains (core/neighborhood-recur s1 1) s1) => true
               (.contains (core/neighborhood-recur s2 1) s2) => true)))

(facts transcribe-rna
       (fact "transcribes RNA into a peptide chain"
             (let [rna "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
                   peptide "MAMAPRTEINSTRING*"]
               (core/transcribe-rna rna) => peptide))
       (fact "also works with DNA"
             (let [dna "ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
                   peptide "MAMAPRTEINSTRING*"]
               (core/transcribe-rna dna) => peptide))
       (fact "probably blows up if the RNA chain isn't a multiple of 3"))

(facts rna-encodes-peptide?
       (fact "if the RNA strand encodes the given pepide"
             (let [rna1 "ATGGCC"
                   rna2 "ATGGCC"
                   rna3 "ATGGAT"
                   rna4 "GATATG"
                   peptide "MA"]
               (core/rna-encodes-peptide? rna1 peptide) => true
               (core/rna-encodes-peptide? rna2 peptide) => true
               (core/rna-encodes-peptide? rna3 peptide) => false
               (core/rna-encodes-peptide? rna4 peptide) => false))
       (fact "stops as early as possible"))

(facts find-encoded-peptides
       (fact "finds encoded peptides in both the forward and reverse compliments of an RNA strand"
             (let [rna "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
                   peptide "MA"
                   strands (sort '("ATGGCC" "GGCCAT" "ATGGCC"))
                   rna2 "TTCGAATCCGTATTGCTTTGGTTGGTTGAGAGAGCGACCTACATTGCTAGTCAGAATAGGTGATTCACACAAAGCTTAGCACCTGGGCAGCACCCCGTGATGTAAACCTATGGGAACTAAGGGAGTCCTGCGGTTTTAGCCAGCAAGCGAGCCGGCAGGAACACTCATACCATCGGACGCGTTTGACGCCTCCCCGGAAAGGAAGTATTTGAGCCTCATTATTACGTATTGCCCGTTAGTCGACAAATCAAGCCCTCGTACGCAGCTTATTCGTACGACGTGGAGGCGTTCCCACGGGCCTAACACGATTGGAACACCACCATAGTAGTGTGGTTCAAATACCTCCTTTGGAGATCTAGAGCTTCACTCTGATTCTAGAGGCAACTTTACAATCGCTCTACGAAATTGTATGGACATCATCAACCGGATATTCTGGGGCGGTAGAATTTCTTTTGTTCGAATCGCTCTAGGCCAGGATCAAATTAATTGAATTGCGGACTCAAGGATCGCGATAGCCGACACATCGGACGCTGTAGAAAGCCAGTCTCTGGATTTAATCCACCCTCTATGTTTGACAAAGCACTAAAACGGGATAGTTTCGGGTGGTATAAGTTTCCCAAGACGATTGCATCGCAATTCATCAACAACCATGAACTTACTGTTTTAGTACTTCCACACACCTTGTTAAATTACGCCTTTACTTCATGTTGCGGTGTGTGTTAGATAGTGTGCAGCTACAAGTCTACCGCCATCGCAGCTCGGGATACCGGCAGATGAGATGGTCCTGAGCTCGTACCGGACTCAAACTTTTTCCTTTACTACCTAGGAATCGCCCATGCGAATTTGTCGGACACACACCATTACATTAACGTCACAACAGCTACTGTTAGAATTTTGCTCTTGCAAATCCTGGAAAGAGTTAAAAAAACTCTTCCGCGCGCCAATAGGGTAAATAATAGATAGCCAGACGGCTGTAAGAGGTGATGACATTTGCAACAATCATGCTGTCGCATCTTCCGCAAGTTCATGTCGCGCCTAGGCAATGGATCTGCGAATGGGGGCCACGGGGTATGAACTACGGAATTCTAAGAAAGTTGCCATCCAGAGTTAAGGGTTTGAGGCTAGTTGCATCGCTGGTAACGAACTACCTCATTACTTGGACGCGCAGTGTGACTTCACTCCTGTATAGCGATGATGCCAAGCAGGAATTAGCAAATCTGAAGAGCGTTTCCAAACTGGCCACTTGGACTGACACCTATCGCGGGGGATTTCAGCGCGTGTCGCTCTCACATGAGAGCTGCCGTCAGGAGCGGTAGAGTTTAGAGAGGAATGCGACAAACTCCCTATTCACCTCTCTGGTGATGTAAGGATATTTACGCTTAGTTCTATGCCAGGCTTAGGGCCTCTCGGAACTTTGGTGAGTCCTTATTAATTGATGCTACCTCTCCCTTACCTTCGCCCCAAGTCACGTAGAAGTACTCAATCCTGCTACATGATAATCAAATATTTCCAACGTTGGGAAATCGGTGACATCACATACTAGTTAAGAAACCACTGTCAGTGAACTTATATCCGGGGGAGAAAATCTACTAACTTACATACGCTGTGCGAGCAGTTTTCATTATAAGAAAATATACTCCCGAGGTACCGCATCAAGCACGACATTCCCGGAGAGCATAACATTTCGGTGCACCTGCTTTTGTGCGCTTGCTTGCGGTTATTTATAAACTACGCACAAGGCGCAAACCGCAGTGCGCATGTTTTCTCCGCCTGGCTAGAACTCGACATTCTCGTCAACGCCAATCTATGTGAGAGGATTTAGACCTCTGTGAAAACGAGTCCCTCTATAGAATAAATACCCAGATGCCAATGGGGGTTCTATCCGATGGCAGTGCATGGAGTGGTGGCTCCAGATTAAGGATGAGGAGAGGTAAAGATAACAGTTCGGTCGCCACGACGCGTTGCCAATCGAAATATCAGTACTAAAAGGCCCACCGCTCCGCTTTAGTCCGACTTTACATCCTGTGGAAATTGTCGAACGGAGGCTACATCGGGCTATATGAGTGTGAAAACCTATACTTCTCGCGTCGTTACTCAGTGCCGGTCTCCTGTTTCCCCCAGTCTTACGTACCCTTATTGATATTTGCTTCACGTTGAAACGTCCTAACGCAGCGTAAAGAGGTGTTTGAACCTCATTACTATAAAATCGCGATCGAAGGTAGACTACATACGCAAACGCCGAAACCCTCAGTTGGCCTTGTTGCAAGTATGGAACGTTGTAAAATTTTTCCTAGACGTTGAGCTATCGGTACAAGGTCGTTAGCGTCCTTACCCTTCACTTATATGCCCGACAAAACGCGGGTCCTAGTGCAGTGGTGGGAGCTTGGAATCCCGCAATACAAGGACAACCTGTATCTCGTTCGGCGTTCCGCGATCACTCGATCCCGAACCACTCCAAGCCTGGTTGATCAGCAAAAGCGGAAGGATGGATAAAGGGCTACTGGTTAATGGATGTAAACTTCCAATGATGAAATCCTGGAAACGAGGGATCGGGTTACGGTGGCGAACGGGGTACGGCAACGTGGCTATCTAGAGCCCGACGTTACGACTCATGTACATGCTGCTACGTGGTTGAAGCTGACGTTCAGATGAAGCAGTACTGAGTCCTAGGGCTTACTACTACTCCAATAGGTCTGGCCGGCCAGATACAAAAGTTCGTGGCGGCTCACCCCCTTTCTGGCGGGTGTAGCTTGCTGACCGGTTTGCTCGATAACACAGGCTAGCGAATAGTAATGAGGTTCGAAAACCTCTTTCCAACGACTGAAAGGGTCTACACGAACTATCTACATTTCCCCGCCCATGTCCTTCCGTCTGGTTGCTTCTGGAGATCCTTTCGCATTATACCGCAGCGTAGTGGCTCTGGCATATATGAAAAAATCCTTCTGTGGGTATTTGTGCCATCACTTATTGTTCGTACCGATATGGGATTACAAGTGCGATGTGATAATAAGCGAAGAAGCCAACATGTTACACTGTTCATGCGCTCCGGGTAATGTGCGGGCACCATGCTCAGTTCCCGCTCGCAGTTGTCACTGTCCCTGTTTCGGCACCATAATCAACATTTCCACGGCCACGCTGGTGAATAACCGAGGATACCGAAGTACAGCAAGAATGAGAGCGGGACTCCTCCATCTGCTTGTAATACGCCTTCAAGATAGTCCATAAAACGGTCGGGGTCTTGTTTCGGACTAGCCGCTTTGAAACGGTGCATAGTTGTGTCAAGTGTGGACATTGGCTTTCTATCCTCGTCAGCGATCCTCGGAAAGACTCGGGCAGTCGCCCCGAATCGTAATTAGGTAGTAGTGCGGCTCAAAAACTTCCTTCGACCTAACCGCTATAATGTTCGTAGATATAAATTTCGTTTCAGTATTAACAGGCGCACCGTATATATACGGAATGGTGTCGCCCCATTAGCTGCTCGCCAATATTTATCTAAGACCGCGCGCGTCTAGCGCCTTTAGTAGTTGCACCCGAGTATAGTAATGGGGTTCGAAGACTTCCTTCGCAAGGCTGCCATACTGTATCACAAGTACTGACGGAGCCCCGGAGGAGTGCAGGATACGGCAAAGGAGACCATTACCGGGGCATGAGTCCAAGTTAGCCCGTTAGGTGAAGGACGCTGATACAATAGTGAATCCGTTACTGAAAGGTTTAGAAGACCGGGGGGCTCGCACTAGGTCCAAATATTATGAACCCTACTCCTGCAACTGAATTGGCCGTCCAGGCGATATTTAAAAGGGGTTACTAGCAGGTTCATCGGAGCCCGTACTCCTTCCGGGCATAGTCGTTCGACGGGTAGAAATTCATCCAGTCGTGCCGGATACCCCGAGAATACCCCTATTTTTTGATCCTTCACCATCATCGTCCGCGGACTCATCTAAGTACCTCAGACCGAAACTGTTATCGTAGCGAAGAGCGAACTCGAATGACATCGCTTGTCCAACAGGGAAAATATGTAAAGTATATGCAGATTATTATAGGAAGATCACAAACTCCATCGCGCCTAGGCCAAAGACTTGCCAGAACAACATCTCTTCCAGAGCAAGGAAGTGTTTGAACCTCACTATTATCGAGAGAAGTCCCATGAATTTATAATAGTGAGGCTCAAAAACTTCCTTCATCGTCGGGCGCTGGGGCGAGCTAGGCTTCCCTAGCCGTCATTACTGTACCCACGCCAAATATTATCGAGTATGACTACGAAGCCTTCACAAGGGGCTGCAGTAACAGACTAACTCACGCACACCGCAACTACTATCCTAAATTGAGTAAGAAAACCCTCGACTATAGCCCTTGAATTAAATTCGCTATTTAAGGAAGACCGCGCTTCCGCTCGCCCGCGTATAGTTTACCAGTCCCAAGAAACATGTGTGGCCAGCCTACCTGAAGAGCTGTGTAAAGATAGATTCTGACATCCTCAAAAAGAAGTTTTCGAGCCGCACTACTACGCACGGAGCTCCGTTATTCAAGGCATGTCGAAGTACAGCGTGGTGGCCTCTCCTGTTCTCCACCCCAGCTAAACCCACCTTGTTCGAATTCGCGCAACTGTATATGACATGAACACTTACAGGGGTTAGAAGTTACTGCAACCAAGATAGAGCTCGTCGAAGTAATAGTGCGGTTCAAAAACTTCCTTCAATTGGTCTCATCACTTAAATTTAAGAGCTATTGTGAGTACAGGTACGGATGCGGCTTCAGTGGATCTTCAGCATTCATTCCTTGTAGTAATGGGGTTCGAAGACTTCCTTGCCAGGGTACCAAACAAGTCTTGCGCATCCTCCTCCCTAAGGAGGTATTTGAACCCCACTATTACCCACGATAGAACATGCAGGGTTTGATAGTGGAACACCTTTTAGAATCTGGGGATAAATTCCCAGGACTAATGTATGGCTGTAGTAATGAGGTTCAAAAACCTCCTTTTCAGGTGGATCGCAGGCCGTGCTGCCTCACAAGCTGGGACGCCGTCCACGGTATAGCCGGCGTCGGCAGTTACTGTGAAATAGCGGAAACTCGATCCCAATATATCATCTTACGTTTGGCGCCCAATAGTCGCCCAGTACCCGTTGACAGTTCTTTAACTCGGCTTAGAACTACTAGACAGGTTCAACCGAACCTTGCCCTAGTTCCCACTCCCGTAATTCATTTGGGTTTGCATCTAGAATACTGGAGGGTGCATAGACGAAACGTGTACGTCGGAGAAAACGTAAGAAAATGACCTAGACTCATAGTAGTGAGGTTCGAAGACTTCCTTTCAGTGAAATCGATCCACCACTCGCCGCGAAGAGATAATAGCATAGAGCACAAGTGCGCGAGTAGAGAAAAAGGCTATCCCAACCGGGCACGTCCTTCGTGTTTGGCGTTTACATACGGCACCCCGTTTCTGCACGTTAACCGTCTAGTATCCAACGGTGGATGGCGGACGCTAGACTATAGATATGAGATATCGAGACCTGGAGCTGGGTGTGGCTGCAGCCCGGGTCATTGCGGGCTGTGAATTCAAGGGCATGTAAACAAGCGTATATCGAACAGTGGATGGGCACCTGCAATACTCACGGTAGAGTTAGCTCACAGGATTCACGTTGAGGACTATGAGTCCCTCTTCGCTAGCAGTCTGGGGGGATATGGAGTTTAATAGCTTGACGTAGTAATGGGGCTCAAACACCTCTTTGTGTGAGCACAGCTACTTGCATTAAGAGATTCTAAACAGCGATCATCTCGGCTATTTCGGGCCAGCCTTTTCGGCAGGATGTTATGTAGCATTTCTGGAAGCTTCCCCCTCGAATCTACTAGTGGTGAGAAGATGCCCCACCGATATTACTCTTTAATCTTGAGAAACCTAAAACCGATCTGACCTCAGACGGGCGGCTCCCACCCGAGGATAAACTCGTCAATAATAATGTGGCTCGAACACTTCTTTTCTCACTAGGCTTTTACGACACGCCACATGTATTTAAGCATCTACCTAACTTGTGTCTGCTGATATACAGCGCATTCTACGCCCAACCTACCAATTACTTCAACGTAGTGCGTGGCTAAAATTCAGGGGAGCTTCATCTCTGTCTTAATTTGAAGGTTCTTCCGGGGCGTTTGGGAATCTTCGTGCCTTTTGCGAGGTTAAGGTATCAAAGAAGTTTTCGAACCACATTATTACCGCCTTAAGCACCGGCGCATCCTGCTCGTGACAACTCTACCCTGCCCTGATAAAGGCACTGAACGTTCCAGAGAGTGCATCATTGACACGCGAGCAGGCCACAGTAGCCACAAGACGTATGGGTGATTATAGAATTGGTGGAGGTGTTGTTAACGATCAGGAGGACATTAGTGGGAGTTAGGAAAGACCCTATGTTCTCTCTATCGCGGACTTGTAACTTGACAAGCAAAAGGGTAAGAGAGCTGCACACCGAAGCAGGCCCTTCCTATACCTGTTTTTCCTACGCGTAGAGAGGAATCCAGAAAGGTGATAATTGGCATTCGATGAAAAAACAGTGTGCCACTGACTTAGTTCTATATGTGAAGAGCCTGTTAGCACGTGACGGCGGCCTTGGTATAGAGCCCTTAATGGTCTCCATCGCGTAGTAATGGGGCTCGAAAACCTCCTTACTTGGGATTGCGTGGCCTCCTTGTGAGTCATACACAAGGCTTAGGGCTATGGGGCGATACACTCCTTTTCGCGGCGCATGGGGCGGTGATGCCTACATAGTAGTAGTGACTGCCTTTCTGGGGGGCTATTTGTGGATGACCAACACCTGACCAGCGATGCAATCGCTAGGGGAGGTACACCTCTCATATGTTACAACAATCACCGAATTGTGTTTCGAATTCGAATCAAGTTTGCGGTGTCGACCAGATCTGGTCTTGCTGCCATACCGGGTTCGCCGCCTCCGGTGGATAGAACTGCATCTTAAGACATCTGGACCCAGCGGTAAGTAGCGGGAAGAGTTTAGAGTCATTCGTACAACTACAGGCTAAGGGCTTACTGGGGAGTTGTTGTAGGGCATAAAGATCGCCCCATGACTTTTCGTACTTTCCCCGATAGTTCACTCGCAGCGAGCTGCGGCTGGGCTTCGCCACACGAGTACGGGCAACATTTATCTCCTCTAATCACTGGGCACCGCGCGAGGAAATAGAAAACCCTAATCAGTGCTCATGGGCGCATCTATTGGTCTCCGCATGCACGATGCCGCGGAGTGCTTAGTTGTCCCTGCATAATCTTCGTAGATGTATAAGAGATTACCTATTTATTCGGTTTCGGTTCTAGACGTACCTTGCCGCATGAGTATAGGCTAATGAACTGAGTTGGCGCCAGAGGGAAAGGCATAATAATGCGGCTCGAATACTTCCTTAAGGAAGTATTCGAACCACATTACTAT"
                   peptide2 "KEVFEPHYY"
                   strands2 (sort '("AAGGAAGTATTTGAGCCTCATTATTAC"
                                    "AAAGAGGTGTTTGAACCTCATTACTAT"
                                    "AAGGAGGTATTTGAACCCCACTATTAC"
                                    "AAAGAAGTTTTCGAACCACATTATTAC"
                                    "AAGGAAGTGTTTGAACCTCACTATTAT"
                                    "AAAGAAGTTTTCGAGCCGCACTACTAC"
                                    "AAGGAAGTATTCGAACCACATTACTAT"
                                    "ATAATAATGCGGCTCGAATACTTCCTT"
                                    "GTAGTAATGGGGCTCGAAAACCTCCTT"
                                    "GTAGTAATGAGGTTCAAAAACCTCCTT"
                                    "GTAGTAATGGGGTTCGAAGACTTCCTT"
                                    "ATAATAGTGAGGCTCAAAAACTTCCTT"
                                    "ATAGTAATGGGGTTCGAAGACTTCCTT"
                                    "GTAGTAGTGCGGCTCAAAAACTTCCTT"
                                    "ATAGTAATGAGGTTCGAAAACCTCTTT"
                                    "ATAATAATGTGGCTCGAACACTTCTTT"
                                    "GTAGTAATGGGGCTCAAACACCTCTTT"
                                    "ATAGTAGTGAGGTTCGAAGACTTCCTT"
                                    "GTAATAGTGCGGTTCAAAAACTTCCTT"
                                    "ATAGTAGTGTGGTTCAAATACCTCCTT"))]

               ; $#!@ if i know, that example's wrong!
               ;(sort (core/format-encoded-peptides (core/find-encoded-peptides rna peptide))) => strands
               (sort (core/format-encoded-peptides (core/find-encoded-peptides rna2 peptide2))) => strands2
               )))

(facts linear-subpeptides
       (fact "returns the set of possible subpeptides, given a linear peptide"
             (let [peptide "ABCD"
                   subpeptides '("" "A" "AB" "ABC" "ABCD" "B" "BC" "BCD" "C" "CD" "D")
                   actual (core/linear-subpeptides peptide)]
               (sort actual) => (sort subpeptides)
               (count actual) => (count subpeptides))))

(facts cyclic-subpeptides
       (fact "returns the set of possible subpeptides, given a cyclic peptide"
             ;; TODO: fix so complete protein configurations are only completed once
             (let [peptide "ABCD"
                   subpeptides '("" "A" "AB" "CDA" "ABC" "DA" "ABCD" "B" "BC" "DAB" "BCD" "C" "CD" "D")
                   actual (core/cyclic-subpeptides peptide)]
               (sort actual) => (sort subpeptides)
               (count actual) => (count subpeptides))))

(facts peptide-mass
       (fact "sums each amino acid's mass in daltons"
             (core/peptide-mass "CANCER") => 676))

(facts linear-subpeptide-masses
       (let [peptide "LYFE"
             peptide-masses [113 163 147 129]
             expected (sort '(0 113 129 147 163 276 276 310 423 439 552))
             actual (core/linear-subpeptide-masses peptide)
             actual-max (apply max actual)]
         (fact "maps the mass sums for each linear subpeptide"
               (count actual) => 11
               (sort actual) => expected
               (count (filter #(= %1 actual-max) actual)) => 1)
         (fact "also works when run on a collection of integers, though excludes possibilities with duplicate amino masses"
               (core/linear-subpeptide-masses peptide-masses) => expected)))

(facts cyclic-subpeptide-masses
       (let [peptide "LYFE"
             peptide-masses [113 163 147 129]
             expected (sort '(0 113 129 147 163 242 276 276 310 389 405 423 439 552))
             actual (core/cyclic-subpeptide-masses peptide)
             actual-max (apply max actual)]
         (fact "maps the mass sums for each cyclic subpeptide"
               (count actual) => 14
               (sort actual) => expected)
         (fact "only count the completely peptide configurations once"
               (count (filter #(= %1 actual-max) actual)) => 1) ;(count peptide)
         (fact "also works when run on a collection of integers, though excludes possibilities with duplicate amino masses"
               (core/cyclic-subpeptide-masses peptide-masses) => expected
               )))

(facts expand-peptides
       (fact "expands peptides by adding each possible amino acid mass to the enumeration, returning a list of vectors"
             (let [peptides [[123]]
                   expected (map (fn [i] [123 i]) core/int-mass-values)]
               (core/expand-peptides peptides) => expected
               ))
       (fact "can be expanded using an 'alphabet' of specified int-mass-values'"
             (let [peptides [[123]]
                   alphabet [1 2 3]
                   expected (map (fn [i] [123 i]) alphabet)]
               (core/expand-peptides peptides :alphabet alphabet) => expected)))

(facts consistent-spectra
       (let [actual-spectra-freq (core/spectra-count (core/cyclic-subpeptide-masses "CANCER"))
             spectra1 '(71 129 288)
             spectra2 '(71 129 288 288)
             spectra3 '(71 129 288 288 288)
             spectra4 '(71 129 300)]
         (fact "spectra are consistent which include one each of the masses"
               (core/consistent-spectra spectra1 actual-spectra-freq) => true)
         (fact "spectra are consistent which include up to count of the actual spectra"
               (core/consistent-spectra spectra2 actual-spectra-freq) => true)
         (fact "spectra are inconsistent when they include more than the count of the actual spectra"
               (core/consistent-spectra spectra3 actual-spectra-freq) => false)
         (fact "spectra are inconsistent when they include frequencies not in the actual spectra"
               (core/consistent-spectra spectra4 actual-spectra-freq) => false)))

(facts inverse-num-cyclic-subpeptides
       (fact "returns the inverse of num-cyclic-subpeptides"
             (core/inverse-num-cyclic-subpeptides 2) => 1
             (core/inverse-num-cyclic-subpeptides 4) => 2
             (core/inverse-num-cyclic-subpeptides 8) => 3
             (core/inverse-num-cyclic-subpeptides 14) => 4
             (core/inverse-num-cyclic-subpeptides 22) => 5
             (core/inverse-num-cyclic-subpeptides 32) => 6
             (core/inverse-num-cyclic-subpeptides 44) => 7
             (core/inverse-num-cyclic-subpeptides 58) => 8))

(facts cyclopeptide-sequencing
       (fact "when running against an empty [0] spectrum, an empty vector is returned"
             (core/cyclopeptide-sequencing [0]) => [])
       (fact "when running against a spectrum with one amino acid, the only match is returned"
             (core/cyclopeptide-sequencing [0 71]) => [[71]])
       (fact "when running against longer spectrums, the correct answer is returned"
             (core/cyclopeptide-sequencing [0 113 128 186 241 299 314 427]) => [[186 128 113] [128 186 113] [186 113 128] [113 186 128] [128 113 186] [113 128 186]])
       (fact "when running against a spectrum that has no solution, it doesn't blow up"
             ))

(facts cyclopeptide-score
       (fact "masses in experimental spectra that match to actual spectra contribute +1 to the score."
             (core/spectral-score (core/cyclic-subpeptide-masses "NQEL") [0 99 113 114 128 227 257 299 355 356 370 371 484]) => 11)
       (fact "masses appearing n times in the actual spectra can only add n to the score"
             (core/spectral-score (core/cyclic-subpeptide-masses "NQEL") [0 114 113 114 128 227 257 114 355 356 370 371 484]) => 11))

(facts trim-leaderboard
       (fact "trims leaderboard recursively"
             (let [peptides
                   {1 [15]
                    2 [14]
                    3 [13]
                    4 [9 10 11 12]
                    5 [6 7 8]
                    6 [4 5]
                    7 [3]
                    8 [2]
                    9 [1]}]
               (core/trim-leaderboard 3 [] peptides) => [1 2 3]
               (core/trim-leaderboard 4 [] peptides) => [1 2 3 4 5]
               (core/trim-leaderboard 5 [] peptides) => [1 2 3 4 5]
               (core/trim-leaderboard 6 [] peptides) => [1 2 3 4 5 6 7 8]
               (core/trim-leaderboard 8 [] peptides) => [1 2 3 4 5 6 7 8])))