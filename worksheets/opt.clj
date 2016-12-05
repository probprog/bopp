;; gorilla-repl.fileformat = 1

;; **
;;; # Optimization Benchmarks
;; **

;; @@
(ns worksheets.opt
  (:require [gorilla-plot.core :as plot]
            [clojure.core.matrix :as mat]
            [bopp.core :refer :all])
  (:use [anglican runtime emit]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; In this worksheet we consider simply using BOPP as an optimizer on some classic Bayesian optimization benchmark problems.  Results for some popular packages on the same problems are presented in the paper.
;; **

;; @@
(defdist factor
  []
  (sample* [this] 0)
  (observe* [this value] value))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>#multifn[print-method 0x27933958]</span>","value":"#multifn[print-method 0x27933958]"}
;; <=

;; **
;;; ## Branin
;; **

;; **
;;; @@f(x) = \left(x\_2-\frac{5.1}{4\pi^2}x\_1^2+\frac{5}{\pi}x\_1-6\right)^2+10(1-\frac{1}{8\pi}\cos(x\_1))+10@@
;;; 
;;; Bounds @@-5 \le x\_1 \le 10@@ and @@0 \le x\_2 \le 15@@.
;;; 
;;; There are three global minimum, @@f(x^\*) = 0.397887@@ at @@x^\* = (-\pi,12,1275), (\pi,2.275) ~\rm{and}~ (9.42478,2.475)@@.
;; **

;; @@
(with-primitive-procedures [factor]
  (defopt branin-opt [] [x1 x2]
   (let [x1 (sample (uniform-continuous -5 10))
         x2 (sample (uniform-continuous 0 15))
         t1 (- (+ (- x2 (* (pow x1 2) (/ 5.1 (* 4 (pow Math/PI 2)))))
                  (* x1 (/ 5 Math/PI)))
               6)
         t2 (* 10
               (- 1 (/ 1 (* 8 Math/PI)))
               (cos x1))
         branin (- 0 (+ (pow t1 2) t2 10))
         prior-correction (+ (log 15) (log 15))]
     (observe (factor) (+ branin prior-correction)))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.opt/branin-opt</span>","value":"#'worksheets.opt/branin-opt"}
;; <=

;; @@
(def branin-samples (->> (doopt :importance branin-opt [] 1 :bo-options {:verbose true})
                		 (take 50) ;; Number of optimization iterations to do
                		 doall
                		 (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[-3.793932187108111 4.95087634462598] [-0.21853856501694713 8.03366707915622] [2.3634544645354287 5.775270545231232] [2.9647360458253846 10.096026043663473] [3.899983271831735 14.627279346693625]]
;;; :initial-log-Zs [-82.41522427305625 -22.19506262585381 -11.086581355853344 -59.51566008695524 -168.65127332062266]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.11617931247633569]
;;; :theta-best [2.3634544645354287 5.775270545231232]     :log-Z-theta-best -11.086581355853344     :mean-theta-best -11.095367896581308     :std-dev-theta-best 0.6674155200461411     :i-best 2
;;; :theta-next [1.7482032213544114 7.1027203658319324]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -20.487332027288623
;;; :log-Z-i-best -11.086581355853344
;;; :theta-mean-best ([2.3634544645354287 5.775270545231232] -11.095367896581308)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.05967694301634994]
;;; :theta-best [2.3634544645354287 5.775270545231232]     :log-Z-theta-best -11.086581355853344     :mean-theta-best -11.096611631764716     :std-dev-theta-best 0.3775621817919005     :i-best 3
;;; :theta-next [4.254583958653293 4.56692785346628]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -14.755318850543087
;;; :log-Z-i-best -11.086581355853344
;;; :theta-mean-best ([2.3634544645354287 5.775270545231232] -11.096611631764716)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.08057700135500676]
;;; :theta-best [2.3634544645354287 5.775270545231232]     :log-Z-theta-best -11.086581355853344     :mean-theta-best -11.5654758366708     :std-dev-theta-best 4.067199680893291     :i-best 4
;;; :theta-next [3.218340263717427 4.856522781795318]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -7.398995947616045
;;; :log-Z-i-best -11.086581355853344
;;; :theta-mean-best ([2.3634544645354287 5.775270545231232] -11.5654758366708)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04817424879307838]
;;; :theta-best [3.218340263717427 4.856522781795318]     :log-Z-theta-best -7.398995947616045     :mean-theta-best -7.4828549862019145     :std-dev-theta-best 3.6038229971496283     :i-best 0
;;; :theta-next [1.1071721967128694 3.585266947980962]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -14.951669345088717
;;; :log-Z-i-best -7.398995947616045
;;; :theta-mean-best ([3.218340263717427 4.856522781795318] -7.4828549862019145)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.opt/branin-samples</span>","value":"#'worksheets.opt/branin-samples"}
;; <=

;; @@
branin-samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.3634544645354287</span>","value":"2.3634544645354287"},{"type":"html","content":"<span class='clj-double'>5.775270545231232</span>","value":"5.775270545231232"}],"value":"[2.3634544645354287 5.775270545231232]"},{"type":"html","content":"<span class='clj-double'>-11.095367896581308</span>","value":"-11.095367896581308"}],"value":"([2.3634544645354287 5.775270545231232] -11.095367896581308)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.3634544645354287</span>","value":"2.3634544645354287"},{"type":"html","content":"<span class='clj-double'>5.775270545231232</span>","value":"5.775270545231232"}],"value":"[2.3634544645354287 5.775270545231232]"},{"type":"html","content":"<span class='clj-double'>-11.096611631764716</span>","value":"-11.096611631764716"}],"value":"([2.3634544645354287 5.775270545231232] -11.096611631764716)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.3634544645354287</span>","value":"2.3634544645354287"},{"type":"html","content":"<span class='clj-double'>5.775270545231232</span>","value":"5.775270545231232"}],"value":"[2.3634544645354287 5.775270545231232]"},{"type":"html","content":"<span class='clj-double'>-11.5654758366708</span>","value":"-11.5654758366708"}],"value":"([2.3634544645354287 5.775270545231232] -11.5654758366708)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>3.218340263717427</span>","value":"3.218340263717427"},{"type":"html","content":"<span class='clj-double'>4.856522781795318</span>","value":"4.856522781795318"}],"value":"[3.218340263717427 4.856522781795318]"},{"type":"html","content":"<span class='clj-double'>-7.4828549862019145</span>","value":"-7.4828549862019145"}],"value":"([3.218340263717427 4.856522781795318] -7.4828549862019145)"}],"value":"[([2.3634544645354287 5.775270545231232] -11.095367896581308) ([2.3634544645354287 5.775270545231232] -11.096611631764716) ([2.3634544645354287 5.775270545231232] -11.5654758366708) ([3.218340263717427 4.856522781795318] -7.4828549862019145)]"}
;; <=

;; **
;;; ## Hartmann 6D
;; **

;; **
;;; @@\begin{align}
;;; f(x) &= -\sum\_{i=1}^{4} \alpha\_i \exp \left(-\sum\_{j=1}^{6} A\_{i,j}(x\_j-P\_{i,j})^2\right), \quad \rm{where} \\\\
;;; \alpha &= (1,1.2,3,3.2)^T \\\\
;;; A &= \left(\begin{matrix}
;;; 10 & 3 & 17 & 3.5 & 1.6 & 8 \\\\
;;; 0.05 & 10 & 17 & 0.1 & 8 & 14 \\\\
;;; 3 & 3.5 & 1.7 & 10 & 17 & 8 \\\\
;;; 17 & 8 & 0.05 & 10 & 0.1 & 14
;;; \end{matrix}\right) \\\\
;;; P &= 10^{-4} \left(
;;; \begin{matrix}
;;; 1312 & 1696 & 5569 & 124 & 8283 & 5886\\\\
;;; 2329 & 4135 & 8307 & 3736 & 1004 & 9991 \\\\
;;; 2358 & 1451 & 3522 & 2883 & 3047 & 6650 \\\\
;;; 4047 & 8828 & 8732 & 5743 & 1091 & 381
;;; \end{matrix}
;;; \right)
;;; \end{align}@@
;;; 
;;; Bounds @@0 \le x\_i \le 1, \; \forall i \in \\{1,\dots,6\\}@@
;;; 
;;; The function has 6 local minima and one global minima at
;;; @@\begin{align}
;;; x^\* &= (0.20169,0.150011,0.476874,0.275332,0.311652,0.6573),
;;; \end{align}@@
;;; with @@f(x^\*) = -3.32237@@.  
;; **

;; @@
(def A
  [[10 3 17 3.5 1.7 8]
   [0.05 10 17 0.1 8 14]
   [3 3.5 1.7 10 17 8]
   [17 8 0.05 10 0.1 14]])

(def P
  [[0.1312 0.1696 0.5569 0.0124 0.8283 0.5886]
   [0.2329 0.4135 0.8307 0.3736 0.1004 0.9991]
   [0.2348 0.1451 0.3522 0.2883 0.3047 0.6650]
   [0.4047 0.8828 0.8732 0.5743 0.1091 0.0381]])

(def alpha
  [1 1.2 3 3.2])

(defn hartman-6d [x1 x2 x3 x4 x5 x6]
  (let [x [x1 x2 x3 x4 x5 x6]
        dxP (mat/pow
             (mat/sub P x)
             2)
        terms (mat/transpose
               (mat/mul A dxP))
        terms-reduced (reduce mat/add terms)
        exp-terms (mat/exp
                   (mat/sub 0 terms-reduced))
        f (- 0
             (reduce + (mat/mul alpha exp-terms)))]
    (- f)))

(with-primitive-procedures [factor hartman-6d]
 (defopt h6 [] [x1 x2 x3 x4 x5 x6]
   (let [x1 (sample (uniform-continuous 0 1))
         x2 (sample (uniform-continuous 0 1))
         x3 (sample (uniform-continuous 0 1))
         x4 (sample (uniform-continuous 0 1))
         x5 (sample (uniform-continuous 0 1))
         x6 (sample (uniform-continuous 0 1))
         f (hartman-6d x1 x2 x3 x4 x5 x6)
         prior-correction 0]
     (observe (factor)
              (+ f prior-correction)))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.opt/h6</span>","value":"#'worksheets.opt/h6"}
;; <=

;; @@

;; @@
