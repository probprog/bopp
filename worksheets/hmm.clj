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
                         :bo-options {:verbose true}) ;; So that BOPP spits out some information to follow its progress
                (take 50) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[2 0.38365382681266014 0.4823430869201011 0.4086743943587312 0.608895907898424 0.6409658550809814] [5 0.03449142841949904 0.8045265641381949 0.5387556639684521 0.49389524928254125 0.3738873607520308] [3 0.8404384066870714 0.352178964968338 0.6369563210635534 0.7453552843779148 0.6931280395832291] [2 0.8206756712957628 0.017563110547761163 0.34978755780360826 0.4229099463109336 0.512254896542613] [3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]]
;;; :initial-log-Zs [-51928.50503927859 -51427.89603971816 -116985.61563774954 -103473.76425664629 -22107.542312660145]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.02450445789016542]
;;; :theta-best [3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]     :log-Z-theta-best -22107.542312660145     :mean-theta-best -22369.775924801885     :std-dev-theta-best 2339.919702536459     :i-best 4
;;; :theta-next [3 0.12593448768732027 0.4050690937553105 0.8136403627433 0.6471216405525084 0.776727666950777]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -25083.89993128104
;;; :log-Z-i-best -22107.542312660145
;;; :theta-mean-best ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22369.775924801885)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.05599812569770439]
;;; :theta-best [3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]     :log-Z-theta-best -22107.542312660145     :mean-theta-best -22276.725984139674     :std-dev-theta-best 1474.4705129100885     :i-best 5
;;; :theta-next [3 0.4412012231020639 0.29840084655229204 0.8131458170575449 0.6264601647735567 0.7019277347625381]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -50813.12216720674
;;; :log-Z-i-best -22107.542312660145
;;; :theta-mean-best ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22276.725984139674)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.023545904771284605]
;;; :theta-best [3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]     :log-Z-theta-best -22107.542312660145     :mean-theta-best -22107.687852866744     :std-dev-theta-best 49.11475680910615     :i-best 6
;;; :theta-next [1 0.4690458668069968 0.5154632316910503 0.6626988804118085 0.7589733086072618 0.7842030637829932]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27148.236208382405
;;; :log-Z-i-best -22107.542312660145
;;; :theta-mean-best ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22107.687852866744)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03712842450428576]
;;; :theta-best [3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]     :log-Z-theta-best -22107.542312660145     :mean-theta-best -22366.528599528843     :std-dev-theta-best 1978.858489557478     :i-best 7
;;; :theta-next [2 0.44759800890789925 0.37252482228192596 0.8729057718488862 0.6708578513387761 0.7965583070599128]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30586.239835756063
;;; :log-Z-i-best -22107.542312660145
;;; :theta-mean-best ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22366.528599528843)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.hmm/samples</span>","value":"#'worksheets.hmm/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-double'>0.4765634418787814</span>","value":"0.4765634418787814"},{"type":"html","content":"<span class='clj-double'>0.31581873014572337</span>","value":"0.31581873014572337"},{"type":"html","content":"<span class='clj-double'>0.9291668739422139</span>","value":"0.9291668739422139"},{"type":"html","content":"<span class='clj-double'>0.7101523845173159</span>","value":"0.7101523845173159"},{"type":"html","content":"<span class='clj-double'>0.6699824433495594</span>","value":"0.6699824433495594"}],"value":"[3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]"},{"type":"html","content":"<span class='clj-double'>-22369.775924801885</span>","value":"-22369.775924801885"}],"value":"([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22369.775924801885)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-double'>0.4765634418787814</span>","value":"0.4765634418787814"},{"type":"html","content":"<span class='clj-double'>0.31581873014572337</span>","value":"0.31581873014572337"},{"type":"html","content":"<span class='clj-double'>0.9291668739422139</span>","value":"0.9291668739422139"},{"type":"html","content":"<span class='clj-double'>0.7101523845173159</span>","value":"0.7101523845173159"},{"type":"html","content":"<span class='clj-double'>0.6699824433495594</span>","value":"0.6699824433495594"}],"value":"[3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]"},{"type":"html","content":"<span class='clj-double'>-22276.725984139674</span>","value":"-22276.725984139674"}],"value":"([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22276.725984139674)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-double'>0.4765634418787814</span>","value":"0.4765634418787814"},{"type":"html","content":"<span class='clj-double'>0.31581873014572337</span>","value":"0.31581873014572337"},{"type":"html","content":"<span class='clj-double'>0.9291668739422139</span>","value":"0.9291668739422139"},{"type":"html","content":"<span class='clj-double'>0.7101523845173159</span>","value":"0.7101523845173159"},{"type":"html","content":"<span class='clj-double'>0.6699824433495594</span>","value":"0.6699824433495594"}],"value":"[3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]"},{"type":"html","content":"<span class='clj-double'>-22107.687852866744</span>","value":"-22107.687852866744"}],"value":"([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22107.687852866744)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-double'>0.4765634418787814</span>","value":"0.4765634418787814"},{"type":"html","content":"<span class='clj-double'>0.31581873014572337</span>","value":"0.31581873014572337"},{"type":"html","content":"<span class='clj-double'>0.9291668739422139</span>","value":"0.9291668739422139"},{"type":"html","content":"<span class='clj-double'>0.7101523845173159</span>","value":"0.7101523845173159"},{"type":"html","content":"<span class='clj-double'>0.6699824433495594</span>","value":"0.6699824433495594"}],"value":"[3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594]"},{"type":"html","content":"<span class='clj-double'>-22366.528599528843</span>","value":"-22366.528599528843"}],"value":"([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22366.528599528843)"}],"value":"[([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22369.775924801885) ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22276.725984139674) ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22107.687852866744) ([3 0.4765634418787814 0.31581873014572337 0.9291668739422139 0.7101523845173159 0.6699824433495594] -22366.528599528843)]"}
;; <=

;; @@

;; @@
