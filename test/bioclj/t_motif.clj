(ns bioclj.t-motif
  (:use [midje.sweet]
        [bioclj.motif]
        [bioclj.string]))

(facts "motif-enumeration"
       (fact "finds common neighbors between all seqs in dna"
             (let [dna ["ATTTGGC" "TGCCTTA" "CGGTATC" "GAAAATT"]
                   ans (sort (map acgt-str-to-64b ["ATA" "ATT" "GTT" "TTT"]))]
               (sort (motif-enumeration 3 1 dna)) => ans)
             (let [dna ["TCTGAGCTTGCGTTATTTTTAGACC"
                        "GTTTGACGGGAACCCGACGCCTATA"
                        "TTTTAGATTTCCTCAGTCCACTATA"
                        "CTTACAATTTCGTTATTTATCTAAT"
                        "CAGTAGGAATAGCCACTTTGTTGTA"
                        "AAATCCATTAAGGAAAGACGACCGT"]
                   ans (sort (map acgt-str-to-64b (clojure.string/split "AAACT AAATC AACAC AACAT AACCT AACTA AACTC AACTG AACTT AAGAA AAGCT AAGGT AAGTC AATAC AATAT AATCC AATCT AATGC AATTC AATTG ACAAC ACACA ACACC ACACG ACACT ACAGA ACAGC ACATC ACATG ACCAT ACCCT ACCGT ACCTA ACCTC ACCTG ACCTT ACGAC ACGAG ACGAT ACGCT ACGGT ACGTC ACGTT ACTAA ACTAG ACTAT ACTCA ACTCC ACTCG ACTCT ACTGA ACTGC ACTGT ACTTA ACTTC ACTTT AGAAA AGAAC AGAAG AGAAT AGACA AGACT AGATA AGATC AGCAT AGCCA AGCGT AGCTA AGCTC AGCTG AGCTT AGGAT AGGTA AGGTC AGTAA AGTAC AGTAT AGTCC AGTCG AGTCT AGTGA AGTTG ATAAA ATAAC ATACA ATACC ATAGA ATATA ATATC ATATG ATATT ATCAG ATCCC ATCCG ATCCT ATCGA ATCGC ATCTA ATCTC ATCTG ATGAC ATGAT ATGCA ATGCC ATGGA ATGGC ATGTA ATGTC ATTAA ATTAC ATTAG ATTAT ATTCA ATTCC ATTCG ATTGA ATTGC ATTGG ATTGT ATTTA ATTTC ATTTG ATTTT CAAAG CAACC CAACT CAAGA CAAGC CAATA CAATT CACAC CACAG CACCT CACGT CACTA CACTT CAGAA CAGAC CAGAT CAGGT CAGTA CAGTC CATAA CATAC CATAG CATAT CATCC CATCT CATGA CATGT CATTA CATTG CATTT CCAAG CCATA CCATG CCATT CCCGT CCCTA CCCTT CCGAA CCGAC CCGAT CCGCT CCGGT CCGTA CCGTC CCGTG CCGTT CCTAC CCTAT CCTCA CCTCC CCTTA CCTTC CCTTG CCTTT CGAAA CGAAG CGACA CGACT CGAGT CGATA CGATG CGATT CGCAA CGCAT CGCCA CGCGA CGCTA CGCTC CGCTT CGGAC CGGAT CGGCA CGGTA CGGTC CGGTT CGTAA CGTAC CGTCA CGTCG CGTCT CGTTA CGTTT CTAAC CTAAG CTAAT CTACA CTACC CTACG CTACT CTAGA CTAGC CTAGG CTAGT CTATA CTATC CTATG CTATT CTCAT CTCCG CTCGT CTCTA CTCTT CTGAA CTGAG CTGCA CTGCC CTGTA CTGTT CTTAA CTTAC CTTAG CTTAT CTTCA CTTGA CTTTA CTTTC CTTTG CTTTT GAAAT GAACA GAACT GAAGT GAATG GAATT GACAC GACAT GACCA GACCT GACGT GACTT GAGAA GAGAT GAGCT GATAA GATAC GATAG GATAT GATCA GATCC GATCG GATCT GATGT GATTA GATTC GATTG GATTT GCAAT GCACT GCATC GCATT GCCAT GCCGT GCCTA GCCTT GCGAT GCGGT GCGTC GCGTT GCTAA GCTAC GCTAG GCTAT GCTGA GCTGT GCTTA GCTTT GGAAT GGACA GGATA GGATC GGATT GGCTA GGGAT GGTAC GGTAG GGTAT GGTCA GGTCG GGTTA GTAAA GTAAG GTACA GTACC GTACG GTAGA GTATA GTATC GTATG GTATT GTCAA GTCAG GTCCG GTCCT GTCGA GTCGC GTCGT GTCTA GTCTG GTGAA GTGAG GTGCA GTGCG GTTAA GTTAC GTTAG GTTAT GTTCA GTTCC GTTCG GTTGA GTTTA TAAAC TAAAG TAACA TAACC TAACT TAAGA TAAGC TAATA TAATC TACAC TACAG TACCC TACCG TACCT TACGA TACGC TACGT TACTA TACTC TACTG TAGAA TAGAC TAGAG TAGAT TAGCC TAGCG TAGGA TAGTC TATAA TATAC TATAT TATCA TATCC TATCG TATGA TATGC TATGG TATGT TATTA TATTG TCAAC TCAAT TCACC TCACG TCACT TCAGA TCATA TCATG TCCAA TCCAC TCCAG TCCAT TCCCA TCCCT TCCGA TCCGC TCCGT TCCTA TCCTG TCCTT TCGAA TCGAC TCGAT TCGCC TCGCT TCGGA TCGGC TCGGG TCGGT TCGTC TCTAC TCTAG TCTAT TCTCC TCTCT TCTGG TCTGT TCTTA TCTTT TGAAA TGAAC TGAAT TGACA TGACC TGACT TGAGA TGAGC TGAGT TGATA TGATC TGATG TGATT TGCAA TGCAC TGCAG TGCAT TGCCA TGCCG TGCCT TGCGA TGCGT TGCTT TGGAA TGGAT TGGTA TGTAA TGTAG TGTAT TGTCC TGTCG TGTGG TGTTA TTAAA TTAAC TTAAG TTAAT TTACA TTACC TTACG TTACT TTAGA TTAGC TTAGG TTAGT TTATA TTATC TTATG TTATT TTCAA TTCAC TTCAT TTCCA TTCCC TTCCT TTCGA TTCGG TTCGT TTCTA TTCTG TTGAA TTGAC TTGAG TTGAT TTGCA TTGCG TTGGA TTGGG TTGTG TTTAA TTTAC TTTAG TTTAT TTTCA TTTCC TTTCG TTTGA TTTGG TTTTA TTTTG" #" ")))]
               (sort (motif-enumeration 5 2 dna)) => ans)))

(defn submit-motif-enumeration
  [k d data]
  (let [ans (bioclj.motif/motif-enumeration 5 1 (clojure.string/split data #" "))]
    (clojure.string/join " " (map (comp (partial apply str)
                                        (partial bioclj.string/acgt-64b-to-str 5))
                                  (vec ans)))))


(facts "median-string-distance")
(facts "median-string")

(facts "motif-profile-probability-for-kmer"
       (let [k 5
             kmer (acgt-str-to-64b "TTTAT")
             dna (acgt-get-64b-kmers k "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT")
             profile [0.2 0.2 0.3 0.2 0.3
                      0.4 0.3 0.1 0.5 0.1
                      0.3 0.3 0.5 0.2 0.4
                      0.1 0.2 0.1 0.1 0.2]]
         (motif-probability-for-kmer profile k kmer) => (* 0.1 0.2 0.1 0.2 0.2)))

(facts "motif-profile-most-probable-kmer"
       (let [k 5
             dna (acgt-get-64b-kmers k "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT")
             profile [0.2 0.2 0.3 0.2 0.3
                      0.4 0.3 0.1 0.5 0.1
                      0.3 0.3 0.5 0.2 0.4
                      0.1 0.2 0.1 0.1 0.2]
             ans (motif-profile-most-probable-kmer profile k dna)]
         (apply str (acgt-64b-to-str k (:kmer ans))) => "CCGAG")
       (let [k 6
             dna (acgt-get-64b-kmers k "TGCCCGAGCTATCTTATGCGCATCGCATGCGGACCCTTCCCTAGGCTTGTCGCAAGCCATTATCCTGGGCGCTAGTTGCGCGAGTATTGTCAGACCTGATGACGCTGTAAGCTAGCGTGTTCAGCGGCGCGCAATGAGCGGTTTAGATCACAGAATCCTTTGGCGTATTCCTATCCGTTACATCACCTTCCTCACCCCTA")
             profile [0.364 0.333 0.303 0.212 0.121 0.242
                      0.182 0.182 0.212 0.303 0.182 0.303
                      0.121 0.303 0.182 0.273 0.333 0.303
                      0.333 0.182 0.303 0.212 0.364 0.152]
             ans (motif-profile-most-probable-kmer profile k dna)]
         (apply str (acgt-64b-to-str k (:kmer ans))) => "TGTCGC"))

(facts "motif-profile-consensus-kmer-generator"
       (fact "produces a list of arrays that contain the most probable nucleotides at each index"
             (let [k 5
                   dna (acgt-get-64b-kmers k "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT")
                   profile [0.2 0.2 0.3 0.2 0.3
                            0.4 0.3 0.1 0.5 0.1
                            0.3 0.3 0.5 0.2 0.4
                            0.1 0.2 0.1 0.1 0.2]
                   ans (motif-profile-consensus-kmer-generator dna k profile)]
               (apply str (acgt-64b-to-str k (motif-profile-top-consensus-kmer dna ans))) => "CCGCG")
             (let [k 6
                   dna (acgt-get-64b-kmers k "TGCCCGAGCTATCTTATGCGCATCGCATGCGGACCCTTCCCTAGGCTTGTCGCAAGCCATTATCCTGGGCGCTAGTTGCGCGAGTATTGTCAGACCTGATGACGCTGTAAGCTAGCGTGTTCAGCGGCGCGCAATGAGCGGTTTAGATCACAGAATCCTTTGGCGTATTCCTATCCGTTACATCACCTTCCTCACCCCTA")
                   profile [0.364 0.333 0.303 0.212 0.121 0.242
                            0.182 0.182 0.212 0.303 0.182 0.303
                            0.121 0.303 0.182 0.273 0.333 0.303
                            0.333 0.182 0.303 0.212 0.364 0.152]
                   ans (motif-profile-consensus-kmer-generator dna k profile)]
               (apply str (acgt-64b-to-str k (motif-profile-top-consensus-kmer dna ans))) => "AAACTC")))

(facts "->MotifProfile"
       (let [k 5
             motif-seqs ["AAACT"
                         "AAATC"
                         "AACAC"
                         "AACAT"
                         "AACCT"
                         "AGCTA"
                         "AACTC"
                         "AACTA"
                         "AACTG"
                         "AAGAA"]
             motifs (sort (map acgt-str-to-64b motif-seqs))
             mp (->MotifProfile k motifs)
             m-totals (motif-totals mp)]
         (fact "motif-totals: given a set of motifs of length k, finds the nucleotide counts for each index 1-k, the :max and :max-nucleotides"
               (let [t1 (first m-totals)
                     t3 (nth m-totals 2)
                     t5 (nth m-totals 4)]

                 (:max t1) => 10
                 (:max-nucleotides t1) => [0]
                 (:counts t1) => [10 0 0 0]

                 (:max t3) => 7
                 (:max-nucleotides t3) => [1]
                 (:counts t3) => [2 7 1 0]

                 (:max t5) => 3
                 (sort (:max-nucleotides t5)) => [0 1 3]
                 (:counts t5) => [3 3 1 3]))
         (fact "profile: given the counts of nucleotides for each index 1-k, returns a profile for those counts"
               (let [pro (profile mp)
                     p1 (first pro)
                     p3 (nth pro 2)
                     p5 (nth pro 4)]

                 p1 => [1.0 0.0 0.0 0.0]
                 p3 => [0.2 0.7 0.1 0.0]
                 p5 => [0.3 0.3 0.1 0.3]))

         (fact "consensus-strings: returns all of the possible consensus strings for this profile"
               (let [ans (map acgt-str-to-64b ["AACTA" "AACTC" "AACTT"])
                     con-strings (consensus-strings mp)]
                 (count con-strings) => 3
                 (sort con-strings) => ans))
         (fact "consensus-string: returns the most-likely consensus string for this profile"
               (let [ans (acgt-str-to-64b "AACTA")
                     con (consensus mp)]
                 con => ans))
         (fact "score: given the counts of nucleotides, returns the sum of differences between that nuc"
               (score mp) => 16)))

(facts "greedy-motif-search")
