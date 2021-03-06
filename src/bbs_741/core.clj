 (ns bbs-741.core
  (:gen-class)
  (:use clojure.contrib.core)
  (:use clojure.tools.cli)
   (:require [bbs-741.lib :as lib]))


(defn -main [& args]
  (let [[options args banner] (cli args    
    ["-h" "--help" "Show help" :default false :flag true]
    ["-1" "--first-sequence" "The first sequence to align" :default ""]
    ["-2" "--second-sequence"  "The second sequence to align" :default ""]
    ["-m" "--match" "the match score" :parse-fn #(Integer. %) :default 1]
    ["-i" "--mismatch" "the mismatch score(usually negative)" :parse-fn #(Integer. %) :default 0]
    ["-g" "--gap" "the gap score(usually negative)" :parse-fn #(Integer. %) :default 0])
    pi 3.15
 {:keys [help first-sequence second-sequence match mismatch gap]} options]  ;;end of initital let


(cond 
 (or (or (= true help)(= nil first-sequence)) (= nil second-sequence))  (lib/help) 
   :else (lib/run first-sequence second-sequence match mismatch gap) 
   );;end of cond 
    );;end of let for args
  );;end of main function

                                                                                      
                                                                                      


