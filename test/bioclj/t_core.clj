(ns bioclj.t-core
  (:use midje.sweet)
  (:require [bioclj.core :as core]))

(facts permute-domain
       (fact "Permuting strings works"
             (core/combinate-domain "AC") => '((\A \A) (\A \C) (\C \A) (\C \C)))

       (fact "Permuting arrays works"
             (core/combinate-domain [\A \B]) => '((\A \A) (\A \B) (\B \A) (\B \B))))

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
       (fact "sorting kmers maps works for short kmers maps"
             (let [dataset "GTACTCAAGGCTAATAGTCGCTAATAGTCGCTAATAGTCGCCTCTGTCCGCTAATAGTCACAGATTCACAGATTCGCTAATAGTCGTACTCAAGACAGATTCGCCTCTGTCCGTACTCAAGGCTAATAGTCACAGATTCGCTAATAGTCGTACTCAAGTCGGAGTAAAGCCTCTGTCCGCTAATAGTCACAGATTCACAGATTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCGTACTCAAGTCGGAGTAAATCGGAGTAAAGCCTCTGTCCTCGGAGTAAAGTACTCAAGACAGATTCGCCTCTGTCCGCCTCTGTCCGTACTCAAGGCCTCTGTCCGCCTCTGTCCTCGGAGTAAAGTACTCAAGTCGGAGTAAAACAGATTCGTACTCAAGACAGATTCTCGGAGTAAAGCCTCTGTCCACAGATTCACAGATTCGTACTCAAGACAGATTCGCCTCTGTCCGCCTCTGTCCACAGATTCGCTAATAGTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCACAGATTCGTACTCAAGGTACTCAAGACAGATTCGTACTCAAGTCGGAGTAAAGCCTCTGTCCTCGGAGTAAAGTACTCAAGTCGGAGTAAAGCTAATAGTCTCGGAGTAAAGCCTCTGTCCGCTAATAGTCGTACTCAAGACAGATTCTCGGAGTAAAGTACTCAAGTCGGAGTAAAGTACTCAAGGCCTCTGTCCTCGGAGTAAATCGGAGTAAAGTACTCAAGGCCTCTGTCCGTACTCAAGGCTAATAGTCGCCTCTGTCCACAGATTCTCGGAGTAAAGCTAATAGTCGTACTCAAGACAGATTCGCTAATAGTCGCCTCTGTCCGCCTCTGTCC"
                   freqmap (core/sorted-freq-kmer {"CCGCCTCTGTCCGT" 8, "GAGTAAAACAGATT" 1, "TAAAGTACTCAAGT" 3, "ATAGTCACAGATTC" 4, "CACAGATTCGTACT" 2, "TAAAGCCTCTGTCC" 8})
                   results (core/sorted-freq-kmer (core/kmer-frequency 14 dataset))]

             (first freqmap) => ["TAAAGCCTCTGTCC" 8]
             (second freqmap) => ["CCGCCTCTGTCCGT" 8]
             (first results) => ["TCGGAGTAAAGCCT" 8]
             (second results) => ["TAAAGCCTCTGTCC" 8])))

(facts reverse-complement
       (fact "returns the reverse complement of a strand of nucleotides, regardless of whether it's in upper/lower case."
             (let [strand "TCGGAGTAAAGCCT"
                   strand2 "TcgGAgtAaaGcCT"]
               (core/reverse-complement strand) => "AGGCTTTACTCCGA"
               (core/reverse-complement strand2) => "AGGCTTTACTCCGA")))

(facts kmers-indices
       (fact "returns a hashmap with vectors containing the indices substrings are found"
             (let [freq (core/kmer-indices 2 "ABCDCD")]
               (freq "DC") => [3]
               (freq "CD") => [2 4]
               (freq "AB") => [0])))