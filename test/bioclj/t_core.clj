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
               (core/hamming-distance s1 s1 :d 0 :dmax 0) => 1
               )))

(facts neighborhood-recur)
;; why the $#@! is this wrong??!/!/!/!?? lulz

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
                   strands ["ATGGCC" "GGCCAT" "ATGGCC"]]
               (core/find-encoded-peptides rna peptide) => strands)
             ))



