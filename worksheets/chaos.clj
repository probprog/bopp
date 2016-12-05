;; gorilla-repl.fileformat = 1

;; **
;;; # Chaos
;; **

;; @@
(ns worksheets.chaos
    (:require [gorilla-plot.core :as plot]
           [anglican.core :refer [doquery]]
           [bopp.core :refer :all]
           [bopp.helper-functions :refer [argmax]]
           [clojure-csv.core :refer :all]
           [clojure.data.csv :as csv]
           [clojure.java.io :as io]
           [clojure.core.matrix :as m
            :refer
            [matrix identity-matrix zero-vector
             shape inverse transpose
             mul mmul add sub div]]
           [clojure.data.json :as json])
  (:use
    clojure.repl
        [anglican
          runtime
          emit
          smc
          stat
          [state :only [get-predicts get-log-weight]]
          [inference :only [log-marginal rand-roulette]]]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; ## The Problem
;;; 
;;; As an example application we consider the case of optimizing the transition function parameters of an extended Kalman filter for the tracking of a chaotic attractor.  Chaotic attractors present an interesting case for tracking problems as, although their underlying dynamics are strictly deterministic with bounded trajectories, neighbouring trajectories diverge exponentially. Therefore regardless of the available precision, a trajectory cannot be indefinitely extrapolated to within a given accuracy and probabilistic methods such as the extended Kalman filter must be incorporated. From an empirical perspective, this forms a challenging optimization problem as the target transpires to be multi-modal, has variations at different length scales and has local minima close to the global maximum.
;; **

;; **
;;; Suppose we observe a noisy signal @@y\_t \in \mathbb{R}^{K}, \; t = 1,2,\dots@@ in some @@K@@ dimensional observation space which we believe has a lower dimensional latent space @@x\_t \in \mathbb{R}^{D}@@ corresponding to a chaotic attractor of known type but with unknown parameters.  Given observations up to some time @@T@@, we wish to performance inference over the latent space using an extended Kalman filter as defined by
;;; 
;;; @@\begin{align}
;;; x\_0 \sim & \mathcal{N} \left(\mu\_0, \sigma\_0 I\right) \\\\
;;; x\_t = & A \left(x\_{t-1}, \theta\right)+\delta\_{t-1}, \quad \delta\_{t-1} \sim \mathcal{N} \left(0, \sigma\_q I\right) \\\\
;;; y\_t = & C x\_{t}+\varepsilon\_{t}, \quad \varepsilon\_{t} \sim \mathcal{N} \left(0, \sigma\_y I\right)
;;; \end{align}@@
;;; 
;;; where @@I@@ is the identity matrix, @@C@@ is a known @@K \times D@@ matrix,  @@\mu\_0@@ is the expected starting position, and @@\sigma\_0, \sigma\_q@@ and @@\sigma\_y@@ are all scalars which are assumed to be known.  The transition function @@A \left(\cdot,\cdot\right)@@ is
;;; 
;;; @@\begin{align}
;;; 	x\_{t,1} = & \; \sin \left(\beta x\_{t-1,2}\right)-\cos\left(\frac{5x\_{t-1,1}}{2}\right)x\_{t-1,3}  \\\\
;;;     	x\_{t,2} = & \;-\sin \left(\frac{3x\_{t-1,1}}{2}\right)x\_{t-1,3}-\cos\left(\eta x\_{t-1,2}\right) \\\\
;;; x\_{t,3} = & \; \sin \left(x\_{t-1,1}\right)
;;; \end{align}@@
;;; 	
;;; 
;;; corresponding to a type of Pickover attractor with unknown parameters @@\theta = \\{\beta,\eta\\}@@ which we wish to optimize.  Note that @@\eta@@ and @@-\eta@@ will give the same behaviour. Lets start by defining some helper functions, defining the variables, and importing the data.                         
;; **

;; @@
(defn get-chaos-data [n-points]
  (->> (str "data/chaos/y1.csv")
       io/resource
       io/reader
       slurp
       csv/read-csv
       (into [])
       (mapv #(mapv read-string %))
       (take n-points)
       (into [])))

(def T 100)
(def observations (into [] (get-chaos-data T)))

(defn A
  [beta nu x y z]
  "Chaotic transition according to the pickover attractor"
  [(- (sin (* beta y)) (* z (cos (* 5/2 x))))
   (- (* z (sin (* -3/2 x))) (cos (* nu y)))
   (sin x)])

;; Define the other parameters

(def C
  (transpose (matrix
   [[0.0243087960714491	0.0168113871672383	0.0691772445397626	1.57387136968515e-05	0.303732278352682	0.0618563110870659	4.95623561031004e-06	1.89039711478951e-08	0.00113995322996380	0.00515810614718469	0.00341658550834822	7.85881265878172e-07	1.30076659987469e-16	0.291018005156189	3.67275441629293e-10	3.20483237445295e-05	1.63410289930715e-18	1.65003111575406e-05	0.211775125663058	0.0115361583403371]
    [2.03156858097514e-14	1.61789389343840e-08	0.0518864591774790	3.77174213741851e-10	0.000580831681939922	7.97405604743620e-30	0.0660259803416299	6.29137955059571e-12	0.000500238717434247	2.78157203235329e-05	3.67744070950491e-05	0.706675337480128	0.00104581188589214	7.94060469213274e-06	0.0634211844981511	7.77022364283257e-19	0.0558005175443096	0.0539415131008550	3.43420457005632e-05	1.52362319447776e-05]
    [2.02934951680463e-05	0.0819514426812765	0.00587746530272417	1.12942910915840e-08	7.07682025375137e-07	0.000219700599512478	0.229511658185978	0.209524004136316	2.35546751056699e-19	0.286199014662333	0.00126960750520102	0.0892967894625532	0.00508817870175021	2.91561618504951e-24	0.0895862772514099	6.47343828059709e-06	5.67384242548235e-09	0.000224094066002085	0.00122301803539802	1.25782593868869e-06]])))

(def mu0 [0 0 0])
(def sig-0 1)
(def sig-q 0.01)
(def sig-y 0.2)

;; For speed, define a custom observation distribution
(defdist obs-dist
  [x] ; distribution parameters
  [dist (normal 0 sig-y)]        ; auxiliary bindings
  (sample* [this] (mapv #(+ % (sample* dist)) (mmul C x)))
  (observe* [this value]
           (reduce + (map #(observe* dist (- %1 %2))
                          (mmul C x)
                          value))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>#multifn[print-method 0x27933958]</span>","value":"#multifn[print-method 0x27933958]"}
;; <=

;; **
;;; ## Solution using BOPP
;;; 
;;; We can solve this problem using BOPP by coding the problem and running automatic MMAP estimation to optimize @@\theta@@.
;; **

;; @@
(with-primitive-procedures [matrix add obs-dist A]
   (defopt kalman-chaos-opt
    [observations] ;; Fixed inputs
    [beta nu] ;; Parameters to optimize
    (let [beta (sample (uniform-continuous -3 3)) ;; These effectly represent the bounds on our optimization
          nu (sample (uniform-continuous 0 3))
          delta-dist (normal 0 sig-q)
          trans-sample (fn [x] (into [] (map #(+ % (sample delta-dist)) ;; mapv not cps transformed by default
                       		  		 		  (A beta nu (first x) (second x) (nth x 2)))))
          initial-vals (add mu0 (into [] (repeatedly 3 #(sample (normal 0 sig-y)))))]
      (reduce (fn [states obs]
                (let [prev-state (peek states);; sample next state
                      state (if prev-state
                                (trans-sample prev-state)
                                initial-vals)]
                  (observe (count states) (obs-dist state) obs) ;; observe next data point (when available)
                  (conj states state))) ;; append state to sequence and continue with next obs       
         [] ;; start with empty sequence         
         observations)))) ;; loop over data, return states
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.chaos/kalman-chaos-opt</span>","value":"#'worksheets.chaos/kalman-chaos-opt"}
;; <=

;; **
;;; Next we can call BOPP directly on this to optimize beta and nu
;; **

;; @@
(def samples (->> (doopt :smc 
                         kalman-chaos-opt
                         [observations]
                         100 ;; Number of particles
                         :bo-options {:verbose true}) ;; So that BOPP spits out some information to follow its progress
                (take 100) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[2.898605743257468 2.976789838981034] [0.9321666449902151 1.6530519322024295] [0.9232603298954776 1.4435217272425034] [1.7401972121163336 0.640260653667087] [-1.6115486118783764 0.4938480046732088]]
;;; :initial-log-Zs [46.159521525763104 466.09240763528584 530.9648572195993 -403.9185822528781 577.7623794845746]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07935614898660472]
;;; :theta-best [-1.6115486118783764 0.4938480046732088]     :log-Z-theta-best 577.7623794845746     :mean-theta-best 577.4264344509007     :std-dev-theta-best 8.850997871015775     :i-best 4
;;; :theta-next [-1.770938791817185 0.7181611862569413]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 753.3775596333159
;;; :log-Z-i-best 577.7623794845746
;;; :theta-mean-best ([-1.6115486118783764 0.4938480046732088] 577.4264344509007)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.06400354036441255]
;;; :theta-best [-1.770938791817185 0.7181611862569413]     :log-Z-theta-best 753.3775596333159     :mean-theta-best 752.0117491094121     :std-dev-theta-best 20.756508929746364     :i-best 0
;;; :theta-next [-2.5306307968431003 0.8879770358270882]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 803.1047038159323
;;; :log-Z-i-best 753.3775596333159
;;; :theta-mean-best ([-1.770938791817185 0.7181611862569413] 752.0117491094121)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07341411258356674]
;;; :theta-best [-2.5306307968431003 0.8879770358270882]     :log-Z-theta-best 803.1047038159323     :mean-theta-best 803.0772263468535     :std-dev-theta-best 2.741344054630249     :i-best 0
;;; :theta-next [-1.6884899274367589 1.051959328053588]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 862.6707080256344
;;; :log-Z-i-best 803.1047038159323
;;; :theta-mean-best ([-2.5306307968431003 0.8879770358270882] 803.0772263468535)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.049082563718219584]
;;; :theta-best [-1.6884899274367589 1.051959328053588]     :log-Z-theta-best 862.6707080256344     :mean-theta-best 862.3037799355712     :std-dev-theta-best 10.916587178383756     :i-best 0
;;; :theta-next [-1.1723829583242666 1.3942295005824406]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 765.2437725022246
;;; :log-Z-i-best 862.6707080256344
;;; :theta-mean-best ([-1.6884899274367589 1.051959328053588] 862.3037799355712)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.chaos/samples</span>","value":"#'worksheets.chaos/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.6115486118783764</span>","value":"-1.6115486118783764"},{"type":"html","content":"<span class='clj-double'>0.4938480046732088</span>","value":"0.4938480046732088"}],"value":"[-1.6115486118783764 0.4938480046732088]"},{"type":"html","content":"<span class='clj-double'>577.4264344509007</span>","value":"577.4264344509007"}],"value":"([-1.6115486118783764 0.4938480046732088] 577.4264344509007)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.770938791817185</span>","value":"-1.770938791817185"},{"type":"html","content":"<span class='clj-double'>0.7181611862569413</span>","value":"0.7181611862569413"}],"value":"[-1.770938791817185 0.7181611862569413]"},{"type":"html","content":"<span class='clj-double'>752.0117491094121</span>","value":"752.0117491094121"}],"value":"([-1.770938791817185 0.7181611862569413] 752.0117491094121)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.5306307968431003</span>","value":"-2.5306307968431003"},{"type":"html","content":"<span class='clj-double'>0.8879770358270882</span>","value":"0.8879770358270882"}],"value":"[-2.5306307968431003 0.8879770358270882]"},{"type":"html","content":"<span class='clj-double'>803.0772263468535</span>","value":"803.0772263468535"}],"value":"([-2.5306307968431003 0.8879770358270882] 803.0772263468535)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.6884899274367589</span>","value":"-1.6884899274367589"},{"type":"html","content":"<span class='clj-double'>1.051959328053588</span>","value":"1.051959328053588"}],"value":"[-1.6884899274367589 1.051959328053588]"},{"type":"html","content":"<span class='clj-double'>862.3037799355712</span>","value":"862.3037799355712"}],"value":"([-1.6884899274367589 1.051959328053588] 862.3037799355712)"}],"value":"[([-1.6115486118783764 0.4938480046732088] 577.4264344509007) ([-1.770938791817185 0.7181611862569413] 752.0117491094121) ([-2.5306307968431003 0.8879770358270882] 803.0772263468535) ([-1.6884899274367589 1.051959328053588] 862.3037799355712)]"}
;; <=

;; @@

;; @@
