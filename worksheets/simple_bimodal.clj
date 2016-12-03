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
                         :bo-verbose true) ;; So that BOPP spits out some information to follow its progress
                (take 20) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[0.5088777158653182] [0.9883271589327601] [0.4057456579509261] [-0.4331005508065856] [-0.5431494322427227]]
;;; :initial-log-Zs [-41.309854506800306 -34.59220181897066 -42.99518770205505 -42.539876037593615 -40.768639283417365]
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.08764013619110828]
;;; :theta-best [0.9883271589327601]     :log-Z-theta-best -34.59220181897066     :mean-theta-best -34.59223421262662     :std-dev-theta-best 0.007538119008210946     :i-best 1
;;; :theta-next [1.1096450569446823]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.183930176002875
;;; :log-Z-i-best -34.59220181897066
;;; :theta-mean-best ([0.9883271589327601] -34.59223421262662)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.111533241046195]
;;; :theta-best [1.1096450569446823]     :log-Z-theta-best -33.183930176002875     :mean-theta-best -33.184054753326514     :std-dev-theta-best 0.009108275657186976     :i-best 0
;;; :theta-next [1.2495093242553224]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.706490425786974
;;; :log-Z-i-best -33.183930176002875
;;; :theta-mean-best ([1.1096450569446823] -33.184054753326514)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.105045142484742]
;;; :theta-best [1.2495093242553224]     :log-Z-theta-best -31.706490425786974     :mean-theta-best -31.706676295631592     :std-dev-theta-best 0.012792102696741846     :i-best 0
;;; :theta-next [1.384337254609656]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.430396151097135
;;; :log-Z-i-best -31.706490425786974
;;; :theta-mean-best ([1.2495093242553224] -31.706676295631592)
;;; :BO-Iteration 4
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.021249451545865553]
;;; :theta-best [1.384337254609656]     :log-Z-theta-best -30.430396151097135     :mean-theta-best -30.430632369622735     :std-dev-theta-best 0.014137293442199215     :i-best 0
;;; :theta-next [1.4042031727364908]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.254665451852546
;;; :log-Z-i-best -30.430396151097135
;;; :theta-mean-best ([1.384337254609656] -30.430632369622735)
;;; :BO-Iteration 5
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.11345675842589817]
;;; :theta-best [1.4042031727364908]     :log-Z-theta-best -30.254665451852546     :mean-theta-best -30.255684024212666     :std-dev-theta-best 0.01432404660503667     :i-best 0
;;; :theta-next [1.5451908190930657]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.098224993066136
;;; :log-Z-i-best -30.254665451852546
;;; :theta-mean-best ([1.4042031727364908] -30.255684024212666)
;;; :BO-Iteration 6
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.02368554650618615]
;;; :theta-best [1.5451908190930657]     :log-Z-theta-best -29.098224993066136     :mean-theta-best -29.098365492667675     :std-dev-theta-best 0.010233288411887347     :i-best 0
;;; :theta-next [1.5730511109701877]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.88851967678387
;;; :log-Z-i-best -29.098224993066136
;;; :theta-mean-best ([1.5451908190930657] -29.098365492667675)
;;; :BO-Iteration 7
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.042997648872116885]
;;; :theta-best [1.5730511109701877]     :log-Z-theta-best -28.88851967678387     :mean-theta-best -28.889308705509976     :std-dev-theta-best 0.015418713409837132     :i-best 0
;;; :theta-next [1.6256620189788173]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.509450325514248
;;; :log-Z-i-best -28.88851967678387
;;; :theta-mean-best ([1.5730511109701877] -28.889308705509976)
;;; :BO-Iteration 8
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.02389315144200711]
;;; :theta-best [1.6256620189788173]     :log-Z-theta-best -28.509450325514248     :mean-theta-best -28.509629301954966     :std-dev-theta-best 0.009424395234034136     :i-best 0
;;; :theta-next [1.9682738058088334]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -26.58251368764554
;;; :log-Z-i-best -28.509450325514248
;;; :theta-mean-best ([1.6256620189788173] -28.509629301954966)
;;; :BO-Iteration 9
;;; :n-gps-in-acq-function 50
;;; :acq-opt [0.0010197060436863631]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.58256211926715     :std-dev-theta-best 0.009436661819465671     :i-best 0
;;; :theta-next [-1.1968316351232733]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -32.244573854151376
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.58256211926715)
;;; :BO-Iteration 10
;;; :n-gps-in-acq-function 50
;;; :acq-opt [2.6087752982938774E-7]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.58254440777882     :std-dev-theta-best 0.007038375412774474     :i-best 1
;;; :theta-next [-0.9568657695980969]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -34.97663571744175
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.58254440777882)
;;; :BO-Iteration 11
;;; :n-gps-in-acq-function 50
;;; :acq-opt [2.948817359614054E-4]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.582582603882503     :std-dev-theta-best 0.011009876676756371     :i-best 2
;;; :theta-next [-1.6731916567132488]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.186030851403782
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.582582603882503)
;;; :BO-Iteration 12
;;; :n-gps-in-acq-function 50
;;; :acq-opt [1.462977215451319E-4]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.5825829589913     :std-dev-theta-best 0.01017086561638003     :i-best 3
;;; :theta-next [-1.7845953870926892]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.498797745965696
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.5825829589913)
;;; :BO-Iteration 13
;;; :n-gps-in-acq-function 50
;;; :acq-opt [6.353570063715059E-8]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.582545361535463     :std-dev-theta-best 0.006879222794722128     :i-best 4
;;; :theta-next [-1.8442086418157724]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.17183192716591
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.582545361535463)
;;; :BO-Iteration 14
;;; :n-gps-in-acq-function 50
;;; :acq-opt [5.07617468365555E-4]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.582536086071155     :std-dev-theta-best 0.005965297553499494     :i-best 5
;;; :theta-next [1.8181197709844346]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.311425292178733
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.582536086071155)
;;; :BO-Iteration 15
;;; :n-gps-in-acq-function 50
;;; :acq-opt [2.326240901665793E-7]
;;; :theta-best [1.9682738058088334]     :log-Z-theta-best -26.58251368764554     :mean-theta-best -26.582572783005162     :std-dev-theta-best 0.007210328897662849     :i-best 6
;;; :theta-next [-1.448870998436327]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.871071417002433
;;; :log-Z-i-best -26.58251368764554
;;; :theta-mean-best ([1.9682738058088334] -26.582572783005162)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;examples.simple-bimodal/samples</span>","value":"#'examples.simple-bimodal/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>0.9883271589327601</span>","value":"0.9883271589327601"}],"value":"[0.9883271589327601]"},{"type":"html","content":"<span class='clj-double'>-34.59223421262662</span>","value":"-34.59223421262662"}],"value":"([0.9883271589327601] -34.59223421262662)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.1096450569446823</span>","value":"1.1096450569446823"}],"value":"[1.1096450569446823]"},{"type":"html","content":"<span class='clj-double'>-33.184054753326514</span>","value":"-33.184054753326514"}],"value":"([1.1096450569446823] -33.184054753326514)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.2495093242553224</span>","value":"1.2495093242553224"}],"value":"[1.2495093242553224]"},{"type":"html","content":"<span class='clj-double'>-31.706676295631592</span>","value":"-31.706676295631592"}],"value":"([1.2495093242553224] -31.706676295631592)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.384337254609656</span>","value":"1.384337254609656"}],"value":"[1.384337254609656]"},{"type":"html","content":"<span class='clj-double'>-30.430632369622735</span>","value":"-30.430632369622735"}],"value":"([1.384337254609656] -30.430632369622735)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.4042031727364908</span>","value":"1.4042031727364908"}],"value":"[1.4042031727364908]"},{"type":"html","content":"<span class='clj-double'>-30.255684024212666</span>","value":"-30.255684024212666"}],"value":"([1.4042031727364908] -30.255684024212666)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.5451908190930657</span>","value":"1.5451908190930657"}],"value":"[1.5451908190930657]"},{"type":"html","content":"<span class='clj-double'>-29.098365492667675</span>","value":"-29.098365492667675"}],"value":"([1.5451908190930657] -29.098365492667675)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.5730511109701877</span>","value":"1.5730511109701877"}],"value":"[1.5730511109701877]"},{"type":"html","content":"<span class='clj-double'>-28.889308705509976</span>","value":"-28.889308705509976"}],"value":"([1.5730511109701877] -28.889308705509976)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.6256620189788173</span>","value":"1.6256620189788173"}],"value":"[1.6256620189788173]"},{"type":"html","content":"<span class='clj-double'>-28.509629301954966</span>","value":"-28.509629301954966"}],"value":"([1.6256620189788173] -28.509629301954966)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.58256211926715</span>","value":"-26.58256211926715"}],"value":"([1.9682738058088334] -26.58256211926715)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.58254440777882</span>","value":"-26.58254440777882"}],"value":"([1.9682738058088334] -26.58254440777882)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.582582603882503</span>","value":"-26.582582603882503"}],"value":"([1.9682738058088334] -26.582582603882503)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.5825829589913</span>","value":"-26.5825829589913"}],"value":"([1.9682738058088334] -26.5825829589913)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.582545361535463</span>","value":"-26.582545361535463"}],"value":"([1.9682738058088334] -26.582545361535463)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.582536086071155</span>","value":"-26.582536086071155"}],"value":"([1.9682738058088334] -26.582536086071155)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>1.9682738058088334</span>","value":"1.9682738058088334"}],"value":"[1.9682738058088334]"},{"type":"html","content":"<span class='clj-double'>-26.582572783005162</span>","value":"-26.582572783005162"}],"value":"([1.9682738058088334] -26.582572783005162)"}],"value":"[([0.9883271589327601] -34.59223421262662) ([1.1096450569446823] -33.184054753326514) ([1.2495093242553224] -31.706676295631592) ([1.384337254609656] -30.430632369622735) ([1.4042031727364908] -30.255684024212666) ([1.5451908190930657] -29.098365492667675) ([1.5730511109701877] -28.889308705509976) ([1.6256620189788173] -28.509629301954966) ([1.9682738058088334] -26.58256211926715) ([1.9682738058088334] -26.58254440777882) ([1.9682738058088334] -26.582582603882503) ([1.9682738058088334] -26.5825829589913) ([1.9682738058088334] -26.582545361535463) ([1.9682738058088334] -26.582536086071155) ([1.9682738058088334] -26.582572783005162)]"}
;; <=

;; @@

;; @@
