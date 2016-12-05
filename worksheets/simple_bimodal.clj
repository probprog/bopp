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
                         :bo-options {:verbose true})
                (take 50) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@

;; @@
samples
;; @@
