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

  ;;; THIS IS THE BOTTLENECK!
  
  ([i j M]
     (sel M i j))
  ([i j M N]
     (let [M-row-i (sel M :rows i)
           N-col-j (sel N :cols j)]
       (mmult M-row-i N-col-j)))
  ([i j M N & rest]     
     (let [M-row-i (sel M :rows i)
           MN-row-i (mmult M-row-i N)]
       (apply item-in-product 0 j MN-row-i rest))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- euclidean-distance-squared
  [ x y ]
  (let [delta (- x y)]
    (* delta delta)))

(defn- improve-H-distance
  ;; For M (n,m) W(n,q) H(q,m)
  ;; This requires
  ;; (n * m) * ((q * n * m) + (q * n * q * m)) operations.
  ;;
  ;; which is a TON.
  ;;
  ;; So keep q small.
  
  [M W H]
  (let [W-trans (trans W)
        new-elements
        (for [i (range (nrow H)) j (range (ncol H))]
          (* (sel H i j)
             (item-in-product i j W-trans M)
             (recip (item-in-product i j W-trans W H))))]
    (matrix new-elements (ncol H))))

(defn- improve-W-distance
  [M W H]
  (let [H-trans (trans H)
        new-elements
        (for [i (range (nrow W)) j (range (ncol W))]
          (* (sel W i j)
             (item-in-product i j M H-trans)
             (recip (item-in-product i j W H H-trans))))]
    (matrix new-elements (ncol W))))

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
  "costs should be relative-entropy-costs or a similar structure

   given an n x m matrix M, and integer r,
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
                       " try " tries ", improvement " improvement))]
      (if (or (< threshold-improvement improvement)
              (< max-tries tries))
        [tries improvement next-try]
        (recur next-cost (inc tries) (rest factorings))))))


(defn non-negative-factor-approx
  "A simple UI for factoring a given matrix, with sensible defaults"
  [M r & args]
  (let [iterations (factorings M r)]
    (apply good-enough-factoring iterations args)))