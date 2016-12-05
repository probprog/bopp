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
                         :bo-verbose true) ;; So that BOPP spits out some information to follow its progress
                (take 50) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[1 0.6136355396199027 0.014785302537952516 0.15131100287429766 0.17717556389642186 0.8869910051580678] [3 0.5010145893316937 0.38999005287040234 0.3365759454708628 0.8373857310201516 0.5845293966635656] [5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592] [2 0.9978444454588267 0.42204699210863916 0.628273637445474 0.7931088980723362 0.13063862476668198] [1 0.8440847942620262 0.3527892131995565 0.3260118979546991 0.8701021527511927 0.43806542114866165]]
;;; :initial-log-Zs [-47575.23507188158 -23273.086860314077 -9140.82943270725 -176203.50166435775 -111757.51334959494]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.057961216893897084]
;;; :theta-best [5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592]     :log-Z-theta-best -9140.82943270725     :mean-theta-best -9141.121530816017     :std-dev-theta-best 92.01108917278682     :i-best 2
;;; :theta-next [5 0.42264856253585 0.6101409599625166 0.9889480865713883 0.9431647286375652 0.881063840630045]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -12670.878737500083
;;; :log-Z-i-best -9140.82943270725
;;; :theta-mean-best ([5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592] -9141.121530816017)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03152520140477911]
;;; :theta-best [5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592]     :log-Z-theta-best -9140.82943270725     :mean-theta-best -9142.65599985933     :std-dev-theta-best 393.8365613221697     :i-best 3
;;; :theta-next [2 0.5152129368816105 0.5093737822801269 0.44813799469583054 0.9255435594472519 0.5313593085368673]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -25651.68278298364
;;; :log-Z-i-best -9140.82943270725
;;; :theta-mean-best ([5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592] -9142.65599985933)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.028431555274355192]
;;; :theta-best [5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592]     :log-Z-theta-best -9140.82943270725     :mean-theta-best -9169.889664894989     :std-dev-theta-best 1090.5506043395435     :i-best 4
;;; :theta-next [3 0.5574660210139954 0.2581399437107865 0.3196040794591503 0.8453437576546695 0.9410976281695481]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33055.12415126689
;;; :log-Z-i-best -9140.82943270725
;;; :theta-mean-best ([5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592] -9169.889664894989)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03491231946448294]
;;; :theta-best [5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592]     :log-Z-theta-best -9140.82943270725     :mean-theta-best -9141.318323084968     :std-dev-theta-best 140.3402251462173     :i-best 5
;;; :theta-next [3 0.3898801292858136 0.34927144892070144 0.24128350127010978 0.8197571318316361 0.6948807582571678]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -10568.538705099825
;;; :log-Z-i-best -9140.82943270725
;;; :theta-mean-best ([5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592] -9141.318323084968)
;;; :BO-Iteration 4
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.10292089309996387]
;;; :theta-best [5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592]     :log-Z-theta-best -9140.82943270725     :mean-theta-best -18508.218868094147     :std-dev-theta-best 23526.556502157884     :i-best 6
;;; :theta-next [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -7028.214000208321
;;; :log-Z-i-best -9140.82943270725
;;; :theta-mean-best ([5 0.38356883578290435 0.5722879480986285 0.9753750157793462 0.8664941875207659 0.869746117375592] -18508.218868094147)
;;; :BO-Iteration 5
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.021936095286129508]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -8927.874895775516     :std-dev-theta-best 11630.56741326268     :i-best 0
;;; :theta-next [1 0.4034022195900888 0.5048169512938411 0.5231339034232568 0.6557001295876259 0.5543520205069792]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -22924.241379455212
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -8927.874895775516)
;;; :BO-Iteration 6
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.0372942488355033]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -7193.511324892286     :std-dev-theta-best 3043.387237467751     :i-best 1
;;; :theta-next [3 0.48056098717326634 0.20108181595746902 0.4871818133915819 0.8732462678467087 0.7065833197872079]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -20214.913177158054
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -7193.511324892286)
;;; :BO-Iteration 7
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07562659609528377]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -7160.154427188827     :std-dev-theta-best 2375.7929223981405     :i-best 2
;;; :theta-next [3 0.35954380020331067 0.3801554745610441 0.3255227832690967 0.8069477885894158 0.647911869198202]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -7713.424785341964
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -7160.154427188827)
;;; :BO-Iteration 8
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.061706163679351736]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -7067.555746898201     :std-dev-theta-best 1344.411872138923     :i-best 3
;;; :theta-next [3 0.4558285486483875 0.31344266193991244 0.5122308164247007 0.7971621040066212 0.7457848673987783]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -16697.340770827694
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -7067.555746898201)
;;; :BO-Iteration 9
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03217717347541454]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -9151.506185239443     :std-dev-theta-best 9642.010338327766     :i-best 4
;;; :theta-next [3 0.6410069373057643 0.2608207586738971 0.2616634580976782 0.9563654153727258 0.5963128266163921]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -50915.3207161802
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -9151.506185239443)
;;; :BO-Iteration 10
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.043857546698378234]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -7170.449640068895     :std-dev-theta-best 2237.1201506665666     :i-best 5
;;; :theta-next [3 0.5120159810469371 0.3599851176662992 0.5503903790118868 0.9319178776407794 0.7477528662253953]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -24730.545853241707
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -7170.449640068895)
;;; :BO-Iteration 11
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.057891900632682534]
;;; :theta-best [3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198]     :log-Z-theta-best -7028.214000208321     :mean-theta-best -7069.9377138600685     :std-dev-theta-best 1157.1185483136662     :i-best 6
;;; :theta-next [3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -5151.474550708618
;;; :log-Z-i-best -7028.214000208321
;;; :theta-mean-best ([3 0.3519498466952635 0.42849915634782576 0.3101492018378686 0.8179689386369449 0.6279088368953198] -7069.9377138600685)
;;; :BO-Iteration 12
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.06579017531530418]
;;; :theta-best [3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135]     :log-Z-theta-best -5151.474550708618     :mean-theta-best -5160.823387309763     :std-dev-theta-best 757.1266586400836     :i-best 0
;;; :theta-next [2 0.4272261690117951 0.4625155606401155 0.345516900448912 0.9559562533089512 0.6668516634468601]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -14522.37779293811
;;; :log-Z-i-best -5151.474550708618
;;; :theta-mean-best ([3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135] -5160.823387309763)
;;; :BO-Iteration 13
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.06552630727177405]
;;; :theta-best [3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135]     :log-Z-theta-best -5151.474550708618     :mean-theta-best -5803.089405684266     :std-dev-theta-best 4664.900698489574     :i-best 1
;;; :theta-next [3 0.4206730144273494 0.3327047669234815 0.5495083892393071 0.943448221535512 0.6568033801571028]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -12390.908571233938
;;; :log-Z-i-best -5151.474550708618
;;; :theta-mean-best ([3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135] -5803.089405684266)
;;; :BO-Iteration 14
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.09798617698442373]
;;; :theta-best [3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135]     :log-Z-theta-best -5151.474550708618     :mean-theta-best -5163.176595782454     :std-dev-theta-best 849.7110792861748     :i-best 2
;;; :theta-next [3 0.4234243066828023 0.49657153358531897 0.469289500205289 0.7105042803298833 0.6487744278238664]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -12620.328943290007
;;; :log-Z-i-best -5151.474550708618
;;; :theta-mean-best ([3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135] -5163.176595782454)
;;; :BO-Iteration 15
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04141729984269055]
;;; :theta-best [3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135]     :log-Z-theta-best -5151.474550708618     :mean-theta-best -5161.927697919426     :std-dev-theta-best 1159.306047993774     :i-best 3
;;; :theta-next [2 0.25993694098037046 0.21661112568804858 0.5536458341318001 0.8070202384756248 0.6718501587268982]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -13726.08641247073
;;; :log-Z-i-best -5151.474550708618
;;; :theta-mean-best ([3 0.33254733848443774 0.4803443247227253 0.47597907250766774 0.8778208797595244 0.6047079293873135] -5161.927697919426)
;;; :BO-Iteration 16
;;; :n-gps-in-acq-function 20
;;; 
;; <-

;; @@
samples
;; @@

;; @@

;; @@
