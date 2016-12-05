;; gorilla-repl.fileformat = 1

;; **
;;; # Simple Bimodal Modal
;; **

;; @@
(ns examples.simple-bimodal
  (require [bopp.core :refer :all]
           [anglican.runtime :refer :all]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; In this simple example, we demonstrate the ability of BOPP to carry out unbounded optimization on a 1D problem with a significant prior-posterior mismatch, namely @@p \left(\theta\right)={\rm Normal}(0, 0.5)@@ and @@p \left(Y|\theta\right)={\rm Normal}(5-\left|\theta\right|,0.5)@@ where we want to optimize @@p(\theta | Y)@@ with respec to @@\theta@@.  This requires BOPP to adapt to the target and effectively establish maxima in the presence of multiple modes.  It should be noted that this is a deterministic optimization problem.
;; **

;; @@
(defopt simple-bimodal [y] [theta]
  (let [theta (sample (normal 0 0.5))]
    (observe
     (normal (sqrt (* theta theta)) 0.5) y)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;examples.simple-bimodal/simple-bimodal</span>","value":"#'examples.simple-bimodal/simple-bimodal"}
;; <=

;; @@
(def samples (->> (doopt :importance 
                         simple-bimodal
                         [5]
                         1 ;; Number of particles
                         :bo-options {:verbose true}) ;; So that BOPP spits out some information to follow its progress
                (take 50) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[0.046916593374184544] [-0.37459067270322316] [0.22628511572928373] [0.7255113943155673] [0.26679095290587657]]
;;; :initial-log-Zs [-49.52205550474111 -43.52104193953 -46.13070020510624 -38.04682195210498 -45.40047329738163]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.09868266787245222]
;;; :theta-best [0.7255113943155673]     :log-Z-theta-best -38.04682195210498     :mean-theta-best -38.04752861693298     :std-dev-theta-best 0.026301533605594535     :i-best 3
;;; :theta-next [-0.9700824126765395]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -34.81417400129601
;;; :log-Z-i-best -38.04682195210498
;;; :theta-mean-best ([0.7255113943155673] -38.04752861693298)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04205909519117866]
;;; :theta-best [-0.9700824126765395]     :log-Z-theta-best -34.81417400129601     :mean-theta-best -34.81922958920508     :std-dev-theta-best 0.13905131405722734     :i-best 0
;;; :theta-next [-1.053716348842142]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.818528303715475
;;; :log-Z-i-best -34.81417400129601
;;; :theta-mean-best ([-0.9700824126765395] -34.81922958920508)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.08904278566434207]
;;; :theta-best [-1.053716348842142]     :log-Z-theta-best -33.818528303715475     :mean-theta-best -34.17394342572139     :std-dev-theta-best 1.0131396360302964     :i-best 0
;;; :theta-next [-1.2386245712660944]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.815854594144028
;;; :log-Z-i-best -33.818528303715475
;;; :theta-mean-best ([-1.053716348842142] -34.17394342572139)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.06039336904097954]
;;; :theta-best [-1.2386245712660944]     :log-Z-theta-best -31.815854594144028     :mean-theta-best -31.878051353594497     :std-dev-theta-best 0.2321974931366066     :i-best 0
;;; :theta-next [-1.3258726273901784]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.965883053736025
;;; :log-Z-i-best -31.815854594144028
;;; :theta-mean-best ([-1.2386245712660944] -31.878051353594497)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;examples.simple-bimodal/samples</span>","value":"#'examples.simple-bimodal/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>0.7255113943155673</span>","value":"0.7255113943155673"}],"value":"[0.7255113943155673]"},{"type":"html","content":"<span class='clj-double'>-38.04752861693298</span>","value":"-38.04752861693298"}],"value":"([0.7255113943155673] -38.04752861693298)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-0.9700824126765395</span>","value":"-0.9700824126765395"}],"value":"[-0.9700824126765395]"},{"type":"html","content":"<span class='clj-double'>-34.81922958920508</span>","value":"-34.81922958920508"}],"value":"([-0.9700824126765395] -34.81922958920508)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.053716348842142</span>","value":"-1.053716348842142"}],"value":"[-1.053716348842142]"},{"type":"html","content":"<span class='clj-double'>-34.17394342572139</span>","value":"-34.17394342572139"}],"value":"([-1.053716348842142] -34.17394342572139)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.2386245712660944</span>","value":"-1.2386245712660944"}],"value":"[-1.2386245712660944]"},{"type":"html","content":"<span class='clj-double'>-31.878051353594497</span>","value":"-31.878051353594497"}],"value":"([-1.2386245712660944] -31.878051353594497)"}],"value":"[([0.7255113943155673] -38.04752861693298) ([-0.9700824126765395] -34.81922958920508) ([-1.053716348842142] -34.17394342572139) ([-1.2386245712660944] -31.878051353594497)]"}
;; <=

;; @@

;; @@
