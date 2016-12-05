;; gorilla-repl.fileformat = 1

;; **
;;; # Hidden Markov Model with Unknown Number of States
;; **

;; @@
(ns worksheets.hmm
    (:require [gorilla-plot.core :as plot]
           [anglican.core :refer [doquery]]
           [bopp.core :refer :all]
           [bopp.helper-functions :refer [argmax]]
           [clojure-csv.core :refer :all]
           [clojure.data.csv :as csv]
           [clojure.java.io :as io]
           [clojure.core.matrix :as m]
           [clojure.data.json :as json])
  (:use
    clojure.repl
        [anglican
          runtime
          emit
          smc
          stat
          [state :only [get-result get-log-weight]]
          [inference :only [log-marginal rand-roulette]]]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; ## The Problem
;;; 
;;; We consider a hidden Markov model (HMM) with an unknown number of states.  This example demonstrates how BOPP can be applied to models which conceptually have an unknown number of variables, by generating all possible variables that might be needed, but then leaving some variables unused for some execution traces.  This avoids problems of varying base measures so that the MMAP problem is well defined  and provides a function with a fixed number of inputs as required by the BO scheme.  From the BO perspective, the target function is simply constant for variations in an unused variable.
;;; 
;;; HMMs are Markovian state space models with discrete latent variables.  Each latent state @@x\_t \in\{1,\dots,K\}, t=1,\dots,T@@ is defined conditionally on @@x\_{t-1}@@ through a set of discrete transition probabilities, whilst each output @@y\_t\in\mathbb{R}@@ is considered to be generated i.i.d. given @@x\_t@@.  We consider the following HMM, in which the number of states @@K@@, is also a random variable: 
;;; 
;;; @@\begin{align}
;;; K & \sim \text{Discrete}\{1,2,3,4,5\} \\\\
;;; T\_k &\sim \text{Dirichlet}\\{{1}\_{1:K}\\}, \quad \forall k=1,\dots,K \\\\
;;; \phi\_k &\sim \text{Uniform}[0,1], \quad \forall k=1,\dots,K \\\\
;;; \mu\_0 &\leftarrow \min \\{y\_{1:T}\\} \\\\
;;; \mu\_k &\leftarrow \mu\_{k-1}+\phi\_k \cdot (\max \\{y\_{1:T}\\} -\mu\_{k-1}), \quad \forall k=1,\dots,K \\\\
;;; x\_1 &\leftarrow 1 \\\\
;;; x\_t | x\_{t-1} &\sim \text{Discrete}\\{T\_{x\_{t-1}}\\} \\\\
;;; y\_t | x\_t &\sim\mathcal{N}(\mu(x\_{t-1}),0.2).
;;; \end{align}@@
;;; 
;;; Our experiment is based on applying BOPP to the above model to do MMAP estimation with a single synthetic dataset, generated using @@K=3, \;\mu\_1 = -1, \;\mu\_2 = 0, \;\mu\_3 = 4, \;T\_1 = [0.9,0.1,0], \;T\_2=[0.2,0.75,0.05]@@ and @@T\_3=[0.1,0.2,0.7]@@.  Lets first first load the data and set the known parameters.
;;; 
;; **

;; @@
(defn hmm-data [n]
  (let [ground-truth-x (mapv read-string
                             (into []
                                   (flatten
                                    (csv/read-csv
                                     (slurp
                                      (io/reader
                                       (io/resource (str "data/hmm/gt_x" n ".csv"))))))))
        y (mapv read-string
                (into []
                      (flatten
                       (csv/read-csv
                        (slurp
                         (io/reader
                          (io/resource (str "data/hmm/y" n ".csv"))))))))]
    [ground-truth-x y]))

(def T 500)

(def observations (take T (second (hmm-data 1))))

(defn index->ind
  "converts a collection of indices to a matrix of indicator vectors"
  [values]
  (let [max-v (reduce max values)
    	zero-vec (into [] (repeat (inc max-v) 0))]
    (m/matrix (map #(assoc zero-vec % 1) values))))

(defn square
  [x]
  (m/mul x x))

(def sig 0.2)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.hmm/sig</span>","value":"#'worksheets.hmm/sig"}
;; <=

;; **
;;; ## Solution using BOPP
;;; 
;;; We now use BOPP to optimize both the number of states @@K@@ and the stick-breaking parameters @@\phi_k@@, with full inference performed on the other parameters.  BOPP therefore aims to maximize
;;; 
;;; @@\begin{align}
;;; p(K,\phi\_{k=1:5}|y\_{t=1:T}) = \iint p(K,\phi\_{k=1:5},x\_{t=1:T},T\_{k=1:K}|y\_{t=1:T}) \mathrm{d}x\_{t=1:T} \mathrm{d}T\_{k=1:K}.
;;; \end{align}@@
;;; 
;;; First we define our model using defopt
;; **

;; @@
(defopt hmm-simple-opt
  []
  [n-states phi1 phi2 phi3 phi4 phi5]
  (let [opt-min (apply min observations)
        opt-max (apply max observations)
        n-states (sample (uniform-discrete 1 6))
        init-dist (discrete (apply conj [1] (repeat (dec  n-states) 0)))
        trans-dist-mem (mem (fn [n] (discrete (sample (dirichlet (into [] (repeat n-states 1)))))))

        phi1 (sample (uniform-continuous 0 1))
        phi2 (sample (uniform-continuous 0 1))
        phi3 (sample (uniform-continuous 0 1))
        phi4 (sample (uniform-continuous 0 1))
        phi5 (sample (uniform-continuous 0 1))

        mus (let [left (- 1 phi1)
                  mu1 phi1]
             (if (< n-states 2)
               [mu1]
               (let [diff2 (* phi2 left)
                     left (- left diff2)
                     mu2 (- 1 left)]
                 (if (< n-states 3)
                   [mu1 mu2]
                   (let [diff3 (* phi3 left)
                         left (- left diff3)
                         mu3 (- 1 left)]
                     (if (< n-states 4)
                       [mu1 mu2 mu3]
                       (let [diff4 (* phi4 left)
                             left (- left diff4)
                             mu4 (- 1 left)]
                         (if (< n-states 5)
                           [mu1 mu2 mu3 mu4]
                           (let [diff5 (* phi5 left)
                                 left (- left diff5)
                                 mu5 (- 1 left)]
                             [mu1 mu2 mu3 mu4 mu5])))))))))
        mus (map #(+ opt-min (* % (- opt-max opt-min))) mus)]
   ;; Return states, n-states, mus and transition distribution
  {:states
    (reduce
      (fn [states obs]
        (let [state (sample (trans-dist-mem (peek states)))]
          (observe (normal (nth mus state) sig) obs)
          (conj states state)))
      [(sample init-dist)]
      observations)
   :n-states n-states
   :mus mus
   :transition-dist {0 (trans-dist-mem 0)
                     1 (trans-dist-mem 1)
                     2 (trans-dist-mem 2)
                     3 (trans-dist-mem 3)
                     4 (trans-dist-mem 4)}}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.hmm/hmm-simple-opt</span>","value":"#'worksheets.hmm/hmm-simple-opt"}
;; <=

;; **
;;; Now we carry out the MMAP estimation by calling BOPP
;; **

;; @@
(def samples (->> (doopt :smc 
                         hmm-simple-opt
                         []
                         200 ;; Number of particles
                         :bo-options {:verbose true})
                (take 100) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@

;; @@
samples
;; @@

;; @@

;; @@
