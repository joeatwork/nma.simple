(ns nma.simple
  (:require [incanter.core :as incanter]
            [incanter.stats :as stats])
  (:use [incanter.core :only [matrix sel nrow ncol mmult trans]])
  (:import [cern.colt.matrix.tdouble DoubleMatrix2D]))

"POINT BY POINT IMPLIMENTATION OF
 \"Algorithms for Non-negative Matrix Factorization\"
 By Daniel D. Lee and H. Sebastian Seung, 2001"

(defn- random-matrix
  [rows columns] ;; OR rows columns max-range?
  (let [element-count (* rows columns)
        randoms (stats/sample-uniform element-count :min Double/MIN_VALUE :max 2)]
    (matrix randoms columns)))

(defn- sum [ & args]
  (apply reduce + 0 args))

(defn- recip [x] (/ 1 x))

(defn- item-in-product
  "Completely for keeping giant matrices out of memory.
   Expect this to be much slower (but much less swapular)
   than just (sel i j (mmult M N))"
  [i j M N]
  (let [M-row-i (sel M :rows i)
        N-col-j (sel N :cols j)]
    (mmult M-row-i N-col-j)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- euclidean-distance-squared
  [ x y ]
  (let [delta (- x y)]    
    (* delta delta)))

;; For M (n,m) W(n,q) H(q,m)
;; improve functions run in n*n*q*q*m*m time!
;; So choose small q, and go get some coffee...
(defn- improve-H-distance  
  [M W H]

  (let [W-trans (trans W)
        W-trans-W (mmult W-trans W)
        W-trans-M (mmult W-trans M)
        W-trans-W-H (mmult W-trans-W H)
        H-indexes (for [i (range (nrow H)) j (range (ncol H))] [i j])

        ;; Take a closer look at incanter matrix api, you're
        ;; missing stuff here.
        H-element-at
        (fn [[i j]]
          (* (sel H i j)
             (sel W-trans-M i j)
             (recip (sel W-trans-W-H i j))))
        H-elements (map H-element-at H-indexes)]
    (matrix H-elements (ncol H))))

(defn- improve-W-distance
  [M W H]
  (let [H-trans (trans H)
        H-H-trans (mmult H H-trans)
        M-H-trans (mmult M H-trans)
        W-H-H-trans (mmult W H-H-trans)
        W-indexes (for [i (range (nrow W)) j (range (ncol W))] [i j])
        W-element-at
        (fn [[i j]]
          (* (sel W i j)
             (sel M-H-trans i j)
             (recip (sel W-H-H-trans i j))))
        W-elements (map W-element-at W-indexes)]
    (matrix W-elements (ncol W))))

(def distance-costs
  ^{:doc "Cost fn and minimizers associated with euclidean distance between matrices"}
  {:calculate-cost euclidean-distance-squared
   :improve-W improve-W-distance
   :improve-H improve-H-distance})

(defn- cost-of-approximation
  "by-part-cost is a function of two numbers that should
   return a cost, decreasing the two numbers approach one another"
  [cost-function M W H]
  (let [costs-by-part
        (for [i (range (nrow M)) j (range (ncol M))]
          (let [M-i-j (sel M i j)
                WH-i-j (item-in-product i j W H)]
            (cost-function M-i-j WH-i-j)))]
    (sum costs-by-part)))

(defn factorings
  "given an n x m matrix M, and integer r,
   Produces an (EXPENSIVE, INFINITE) LAZY SEQUENCE OF [ cost W H ]
   Where W is an n x r matrix, H is an r x m matrix,
   and cost is a measure of how closely W H approximates M,
   decreasing as the approximation improves."  
  ([M r]
     (factorings distance-costs M r))  
  ([costs M r]
     (let [M (matrix M)
           M-rows (nrow M)
           M-cols (ncol M)
           W (random-matrix M-rows r)
           H (random-matrix r M-cols)
           
           find-cost (partial cost-of-approximation (costs :calculate-cost))
           improve-W (costs :improve-W)
           improve-H (costs :improve-H)
           
           next-factoring
           (fn [ [_ old-W old-H] ]
             (let [new-H (improve-H M old-W old-H)
                   new-W (improve-W M old-W new-H)                
                   new-cost (find-cost M new-W new-H)]
               [new-cost new-W new-H]))]
       (iterate next-factoring [nil W H]))))

(defn good-enough-factoring
  "Will iterate :max-tries, or until cost reduction after iteration
   is less than :threshold-improvement (whichever comes first) and
   return [ tries, improvement, [ cost, W and H matrices ] ]
   as produced by (factorings costs M r)

   threshold-improvement should be expressed as a number in that describes
   the maximum interesting ratio of cost[i]/cost[i + 1]
   - higher numbers will generally run longer and produce better results,
   a value of \"1\" or more will run (theoretically) forever.

   Will impose a :max-tries of 1000 unless a different argument is provided."  
  [factorings &
   {:keys [max-tries threshold-improvement chatty]
    :or {max-tries 1000
         threshold-improvement 1.0
         chatty false}}]

  ;; Yes, this is some (first (reduce (take max-tries factorings))
  ;; cleverness, but I actually think the below is pretty straightforward.

  ;; PICKING IMPROVEMENT
  ;;   improvement as a part of total cost (upside, this is stateless)
  
  (loop [last-cost nil, tries 0, factorings factorings]
    (let [[next-cost _ _ :as next-try] (first factorings)
          improvement (if (and last-cost next-cost)
                        (/ next-cost last-cost)
                        0)
          _ (when chatty
              (println "Searching for best factorization"
                       "last cost" last-cost
                       "this cost" next-cost
                       "try " tries ",improvement " improvement))]
      (if (or (< threshold-improvement improvement)
              (< max-tries tries))
        [tries improvement next-try]
        (recur next-cost (inc tries) (rest factorings))))))


(defn non-negative-factor-approx
  "A simple UI for factoring a given matrix with all postive or zero elements, with sensible defaults"
  [M r & args]
  (let [M (matrix M)]
    (doseq [ i (range (nrow M)) j (range (ncol M))]
      (when-not (<= 0 (sel M i j))
        (throw (Exception. (str "Element " i ", " j " of matrix == " (sel M i j) ", not zero or positive")))))    
    (let [iterations (factorings M r)
          [tries improvement factoring] (apply good-enough-factoring iterations args)
          [cost W H] factoring
          ]
      [W H])))

