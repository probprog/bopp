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

(def observations (second (hmm-data 1)))

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
                         500 ;; Number of particles
                         :bo-verbose true) ;; So that BOPP spits out some information to follow its progress
                (take 4) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[2 0.6202982032705744 0.7934506624341537 0.8834522606547959 0.905297589219096 0.2911872327070131] [4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] [5 0.867078570613671 0.6677387701790267 0.4574771789581269 0.3641219305836285 0.6693353263615105] [2 0.9682379285927254 0.7893721338796724 0.6806804920437028 0.569917507168052 0.8626491980918682] [1 0.34293186296595524 0.5630892385671835 0.7872748541653938 0.621132084191365 0.393776577614255]]
;;; :initial-log-Zs [-89817.36423377706 -36323.35241817959 -232972.17650582184 -314715.71565813996 -47646.12309537168]
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.04075740651076098]
;;; :theta-best [4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269]     :log-Z-theta-best -36323.35241817959     :mean-theta-best -36327.96946444863     :std-dev-theta-best 508.8141459704806     :i-best 1
;;; :theta-next [4.0 0.5864806913330755 0.6679956343339919 0.8891072623508182 0.5407780925313492 0.21710235797024158]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -75183.92889923378
;;; :log-Z-i-best -36323.35241817959
;;; :theta-mean-best ([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36327.96946444863)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.07698396815517757]
;;; :theta-best [4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269]     :log-Z-theta-best -36323.35241817959     :mean-theta-best -36324.807734202244     :std-dev-theta-best 312.59158453817383     :i-best 2
;;; :theta-next [4.0 0.49437193363581383 0.6197707903047374 0.9600924190588015 0.5296339868063966 0.331819475995425]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -43058.961562189295
;;; :log-Z-i-best -36323.35241817959
;;; :theta-mean-best ([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36324.807734202244)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.07983816598747852]
;;; :theta-best [4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269]     :log-Z-theta-best -36323.35241817959     :mean-theta-best -36324.55161192414     :std-dev-theta-best 304.17862962368105     :i-best 3
;;; :theta-next [4.0 0.4331945091301894 0.6557702102376561 0.8846866753363861 0.5017591055911628 0.32889711468883204]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27101.446504666383
;;; :log-Z-i-best -36323.35241817959
;;; :theta-mean-best ([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36324.55161192414)
;;; :BO-Iteration 4
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.05866916148745695]
;;; :theta-best [4.0 0.4331945091301894 0.6557702102376561 0.8846866753363861 0.5017591055911628 0.32889711468883204]     :log-Z-theta-best -27101.446504666383     :mean-theta-best -27102.552950098703     :std-dev-theta-best 318.8764242820278     :i-best 0
;;; :theta-next [4.0 0.48156346729440574 0.6543270605109099 0.7877920542443401 0.6065929271953667 0.36517926737390805]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -39085.5599608049
;;; :log-Z-i-best -27101.446504666383
;;; :theta-mean-best ([4.0 0.4331945091301894 0.6557702102376561 0.8846866753363861 0.5017591055911628 0.32889711468883204] -27102.552950098703)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.hmm/samples</span>","value":"#'worksheets.hmm/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-double'>0.47131746688572873</span>","value":"0.47131746688572873"},{"type":"html","content":"<span class='clj-double'>0.614396458990947</span>","value":"0.614396458990947"},{"type":"html","content":"<span class='clj-double'>0.8806129326987864</span>","value":"0.8806129326987864"},{"type":"html","content":"<span class='clj-double'>0.45397650897668074</span>","value":"0.45397650897668074"},{"type":"html","content":"<span class='clj-double'>0.3331249262825269</span>","value":"0.3331249262825269"}],"value":"[4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269]"},{"type":"html","content":"<span class='clj-double'>-36327.96946444863</span>","value":"-36327.96946444863"}],"value":"([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36327.96946444863)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-double'>0.47131746688572873</span>","value":"0.47131746688572873"},{"type":"html","content":"<span class='clj-double'>0.614396458990947</span>","value":"0.614396458990947"},{"type":"html","content":"<span class='clj-double'>0.8806129326987864</span>","value":"0.8806129326987864"},{"type":"html","content":"<span class='clj-double'>0.45397650897668074</span>","value":"0.45397650897668074"},{"type":"html","content":"<span class='clj-double'>0.3331249262825269</span>","value":"0.3331249262825269"}],"value":"[4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269]"},{"type":"html","content":"<span class='clj-double'>-36324.807734202244</span>","value":"-36324.807734202244"}],"value":"([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36324.807734202244)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-double'>0.47131746688572873</span>","value":"0.47131746688572873"},{"type":"html","content":"<span class='clj-double'>0.614396458990947</span>","value":"0.614396458990947"},{"type":"html","content":"<span class='clj-double'>0.8806129326987864</span>","value":"0.8806129326987864"},{"type":"html","content":"<span class='clj-double'>0.45397650897668074</span>","value":"0.45397650897668074"},{"type":"html","content":"<span class='clj-double'>0.3331249262825269</span>","value":"0.3331249262825269"}],"value":"[4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269]"},{"type":"html","content":"<span class='clj-double'>-36324.55161192414</span>","value":"-36324.55161192414"}],"value":"([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36324.55161192414)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>4.0</span>","value":"4.0"},{"type":"html","content":"<span class='clj-double'>0.4331945091301894</span>","value":"0.4331945091301894"},{"type":"html","content":"<span class='clj-double'>0.6557702102376561</span>","value":"0.6557702102376561"},{"type":"html","content":"<span class='clj-double'>0.8846866753363861</span>","value":"0.8846866753363861"},{"type":"html","content":"<span class='clj-double'>0.5017591055911628</span>","value":"0.5017591055911628"},{"type":"html","content":"<span class='clj-double'>0.32889711468883204</span>","value":"0.32889711468883204"}],"value":"[4.0 0.4331945091301894 0.6557702102376561 0.8846866753363861 0.5017591055911628 0.32889711468883204]"},{"type":"html","content":"<span class='clj-double'>-27102.552950098703</span>","value":"-27102.552950098703"}],"value":"([4.0 0.4331945091301894 0.6557702102376561 0.8846866753363861 0.5017591055911628 0.32889711468883204] -27102.552950098703)"}],"value":"[([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36327.96946444863) ([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36324.807734202244) ([4 0.47131746688572873 0.614396458990947 0.8806129326987864 0.45397650897668074 0.3331249262825269] -36324.55161192414) ([4.0 0.4331945091301894 0.6557702102376561 0.8846866753363861 0.5017591055911628 0.32889711468883204] -27102.552950098703)]"}
;; <=

;; @@

;; @@
