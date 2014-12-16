(ns bioclj.dynamic-programming)

(defn num-paths-thru-grid [len wid]
  ((fn pascal-triangle [r c]
     (if (or (= r len) (= c wid))
       1
       (+ (pascal-triangle (inc r) c) (pascal-triangle r (inc c)))))
    1 1))
