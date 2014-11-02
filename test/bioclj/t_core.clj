(ns bioclj.t-core
  (:use midje.sweet)
  (:require [bioclj.core :as core]))

(facts permute-domain
       (fact "Permuting strings works"
             (core/combinate-domain "AC") => '((\A \A) (\A \C) (\C \A) (\C \C)))

       (fact "Permuting arrays works"
             (core/combinate-domain [\A \B]) => '((\A \A) (\A \B) (\B \A) (\B \B))))

(facts freqent-kmer
       (fact "counting kmers works for a 1-mers"
             (let [freq (core/frequent-kmer 1 "ABCDCD")]
               (freq "D") => 2
               (freq "C") => 2
               (freq "A") => 1))

       (fact "counting kmers works for a 2-mers"
             (let [freq (core/frequent-kmer 2 "ABCDCD")]
               (freq "DC") => 1
               (freq "CD") => 2
               (freq "AB") => 1)))
