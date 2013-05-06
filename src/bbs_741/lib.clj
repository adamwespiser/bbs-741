(ns bbs-741.lib
    (:require [clojure.string :as str])
    (:require [clojure.pprint :as pp])
    (:import java.io.File))


(defn help [] 
  (println "USAGE: java -jar [<jar file>] -f [first sequence] -s [second sequence] -m [match score] -i [mismatch score] -g [gap score] ")
  (println "direction maxtrix has the following characters represent multidimensional pointers:
           
            - left 
            \\ diagonal 
            J left, up
            v diagonal, up
            > left, diagonal
            * left, diagonal up
"))

(defn calcMatchScore-with-nums[i j seq1 seq2 array match1 mismatch]
  " determine if two seqs match at the given indexes: if so, return match score, else return mismatch score"
  (let [x (seq1 i)
        y (seq2 j)]
    (if (= x y) match1 mismatch )))

(defn direction-to-char[dir-vec]
  "for a given vector of direction symbol, convert to representative character"
  (let [dir-set (set dir-vec)]
    (cond 
      (= #{:end}  dir-set) \O
      (= #{:left} dir-set) \-
      (= #{:up}    dir-set) \|
      (= #{:diag} dir-set) \\
      (= #{:left :up} dir-set) \J
      (= #{:diag :up} dir-set) \v
      (= #{:left :diag} dir-set) \>
      (= #{:left :diag :up} dir-set) \*
      :else " ")))


(defn getSymbol [i j array]
  "get the direction vector at position (i,j) of an input array,includes index check"
  (if (and (> i -1) (> j -1))
  (get-in array [j  i  1])
   nil))

(defn getSymbol-first [i j array]
  "get the first vector at position (i,j) of an input array, includes index check"
  (if (and (> i -1) (> j -1))
  (get-in array [j  i  1 0])
   nil))


(defn getSymbol-random [i j array]
  "get a random vector at position (i,j) of an input array"
  (if (and (> i -1) (> j -1))
    (let [chs   (get-in array [j i 1])
          rnd  (rand-int (count chs))]
     (nth chs rnd))) )


(defn max-for-cell-with-nums[i j seq1 seq2 array scores]
 "determine the max score for position (i,j) in array using seq1 and seq2 to determine match/mismatch" 
 (let [match1 (scores 0)
       mismatch (scores 1)
       gap    (scores 2)
    cell-val (cond 
    (and (= i 0 ) (= j 0)) [0 [:end]] 
    ;;(= i 0) (vector (+ (get-in array [(dec j) i      0]) gap) [:up])
    (= i 0) (vector 0 [:up])
    (= j 0) (vector (+ (get-in array [ j     (dec i) 0] ) gap) [:left])
    true (let [three {:diag (+ (get-in array [(dec j) (dec i) 0]) (calcMatchScore-with-nums i j seq1 seq2 array match1 mismatch))
                      :left (+ (get-in array [ j (dec i) 0]) gap)
                      :up (+ (get-in array [(dec j) i 0] ) gap )}
		max-val (apply max (vals three))
                ;;direction-vec  (mapv first (filter #(= max-val (second %)) three))]
                direction-vec (filterv #(= max-val (three %)) (keys three))  ]
            (vector max-val direction-vec))) ]
	(assoc-in array [j i] cell-val )))

;; (reduce #(setVal (%2 0) (%2 1) %1 999) a (for [x (range 0 5)]([x x])))


(defn constructArray[seq1 seq2]
 "create the score/direction array from two sequences"
  (vec (repeat (count seq2) (vec (repeat (inc (- (count seq1) 1)) [0 nil])))))
;;(calcMatrix "ATCTGAT" "TGCATA")

(defn calcMatrix[sequence1 sequence2 scores]
"determine the alignmet matrix for two sequences as a vector of rows where each row consists of tuples with [current value [optimal neighbour]]"
  (let [seq1 (vec (cons 0 sequence1))
	seq2 (vec (cons 0 sequence2))
	arr (constructArray seq1 seq2)]
  ;;(println seq1)
  ;;(println seq2)
 ;; (println arr)
  (letfn [(reducer-fn [array [i j]];;(println array "-----") 
  (max-for-cell-with-nums i j seq1 seq2 array scores))]
    (reduce reducer-fn arr 
			     (for [x (range 0 (count seq1))  
				   y (range 0 (count seq2))] [x y])))))


;; (calcMatrix "ATCTGAT" "TGCATA" [5 -4 -7])
;; (map println (map (fn[x](map (fn[y](nth y 0)) x )) (calcMatrix "ATCTGAT" "TGCATA" [5 -4 -7])))
;; (map println (map (fn[x](map (fn[y](nth y 1)) x )) (calcMatrix "ATCTGAT" "TGCATA" [5 -4 -7])))
;; (def scores (map (fn[x](map (fn[y](nth y 0)) x )) (calcMatrix "ATCTGAT" "TGCATA" [5 -4 -7])))
;; (def arrows (map (fn[x](map (fn[y](nth y 1)) x )) (calcMatrix "ATCTGAT" "TGCATA" [5 -4 -7])))

(defn printStack[stack]
  "print the current stack as an alignment in two lines"
  (do
     (println "The alignment is:")
     (println (apply str (cons "first seq:  " (interpose " " (map first stack)))))
     (println (apply str (cons "second seq: " (interpose " " (map second stack)))))))


(defn traceBack-inner-random [seq1 seq2 arr stack i j]
   "helper function for traceback, recurrently calls itself with updated position and alignment args"
    (let [choice (getSymbol-random i j arr)]
      (cond (and (= 0 i) (= 0 j)) (printStack stack)
     (= choice :diag) (traceBack-inner-random seq1 seq2  arr (cons (list (seq1 i) (seq2 j)) stack)  (dec i) (dec j))
     (= choice :up)   (traceBack-inner-random seq1 seq2  arr (cons (list \_       (seq2 j)) stack)       i  (dec j))
     (= choice :left) (traceBack-inner-random seq1 seq2  arr (cons (list (seq1 i)  \_     ) stack)  (dec i)      j))))


(defn traceBack [sequence1 sequence2 combinedArray]
  "for sequence 1 and 2 and combinded array, return the optimal alignment"
  (let [seq1 (vec (cons 0 sequence1))
	seq2  (vec (cons 0 sequence2))
	i (- (count seq1) 1)
        j (- (count seq2) 1)
        stack '() ]
   (traceBack-inner-random seq1 seq2 combinedArray stack i j )))

(defn get-m[sequence1 sequence2 match1 mismatch gap]
  "return the score/direction matrix for an alignment"
  (let [comb (calcMatrix sequence1 sequence2 [match1 mismatch gap])]
   comb))



(defn run[sequence1 sequence2 match1 mismatch gap]
  "determine the score/direction matrix for an alignment and print the results"
  (let [comb (calcMatrix sequence1 sequence2 [match1 mismatch gap])
        v (map (fn[x](map (fn[y](nth y 0)) x )) comb)
        d (map (fn[x](map (fn[y](nth y 1)) x )) comb)
        ]
    ;(pp/pprint (map (fn[x](println (eval (cons 'str  (map direction-to-char x))))) d))
    ;;(pp/pprint (filter #(not (nil? %)) (map (fn[x](println (eval (cons 'str  (map direction-to-char x))))) d)))
    (pp/pprint (filter #(not (nil? %)) (map (fn[x](println (apply str  (map direction-to-char x)))) d)))
    (pp/pprint v)
   ;(pp/pprint (map (fn[x](map (fn[y](nth y 1)) x )) comb ))
  (traceBack sequence1 sequence2 comb)
  println (str "Global Alignment Score: " (last (flatten v)))))
