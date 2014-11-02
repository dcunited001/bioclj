(ns bioclj.t-core
  (:use midje.sweet)
  (:require [bioclj.core :as core]))

(facts permute-domain
  (fact "Permuting strings works"
    (let [permuted (core/combinate-domain "AC")]
      permuted
      ))

  (fact "Permuting arrays works"
    (let [permuted (core/combinate-domain [\A \B])]
      permuted
      ))

  (fact "Permuting ACTG works"
    (let [permuted (core/combinate-domain "ACTG")]
      permuted
      )))
