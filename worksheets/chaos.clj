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
;;; {"type":"html","content":"<span class='clj-unkown'>#multifn[print-method 0x601dc1b1]</span>","value":"#multifn[print-method 0x601dc1b1]"}
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
                         :bo-verbose true) ;; So that BOPP spits out some information to follow its progress
                (take 100) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[-1.3046647236636435 0.12992312058833022] [0.5067490852423342 1.879969151250252] [-1.8531507208526143 1.763583387232576] [0.44867877614549867 2.0926295491474285] [-1.66078143971 0.19133833871806782]]
;;; :initial-log-Zs [366.3850490004529 635.9332367480033 718.9111067851228 605.9928083450349 619.0319233803259]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.13504592686852585]
;;; :theta-best [-1.8531507208526143 1.763583387232576]     :log-Z-theta-best 718.9111067851228     :mean-theta-best 715.622733900344     :std-dev-theta-best 14.818472732810097     :i-best 2
;;; :theta-next [-2.398537226969978 0.566750291398498]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 833.0858274604582
;;; :log-Z-i-best 718.9111067851228
;;; :theta-mean-best ([-1.8531507208526143 1.763583387232576] 715.622733900344)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.10359537658810891]
;;; :theta-best [-2.398537226969978 0.566750291398498]     :log-Z-theta-best 833.0858274604582     :mean-theta-best 824.6971612982984     :std-dev-theta-best 39.47473764851365     :i-best 0
;;; :theta-next [-2.3216113772192113 0.34668212542158633]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 746.1346009914962
;;; :log-Z-i-best 833.0858274604582
;;; :theta-mean-best ([-2.398537226969978 0.566750291398498] 824.6971612982984)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.12846037022280773]
;;; :theta-best [-2.398537226969978 0.566750291398498]     :log-Z-theta-best 833.0858274604582     :mean-theta-best 826.5795790501979     :std-dev-theta-best 16.584934292542197     :i-best 1
;;; :theta-next [-2.217440843085523 0.939084122699758]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 911.4169977746108
;;; :log-Z-i-best 833.0858274604582
;;; :theta-mean-best ([-2.398537226969978 0.566750291398498] 826.5795790501979)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.08979017971860331]
;;; :theta-best [-2.217440843085523 0.939084122699758]     :log-Z-theta-best 911.4169977746108     :mean-theta-best 856.0665613110718     :std-dev-theta-best 131.73756363149258     :i-best 0
;;; :theta-next [-2.433789402265904 1.2645978847754078]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 1018.1389974817869
;;; :log-Z-i-best 911.4169977746108
;;; :theta-mean-best ([-2.217440843085523 0.939084122699758] 856.0665613110718)
;;; :BO-Iteration 4
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.1374950797139267]
;;; :theta-best [-2.433789402265904 1.2645978847754078]     :log-Z-theta-best 1018.1389974817869     :mean-theta-best 920.9251877697038     :std-dev-theta-best 172.9346628557226     :i-best 0
;;; :theta-next [-2.8863566850707585 1.1467454560297992]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 818.531429094612
;;; :log-Z-i-best 1018.1389974817869
;;; :theta-mean-best ([-2.433789402265904 1.2645978847754078] 920.9251877697038)
;;; :BO-Iteration 5
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.12544207873144766]
;;; :theta-best [-2.433789402265904 1.2645978847754078]     :log-Z-theta-best 1018.1389974817869     :mean-theta-best 1015.1380098995651     :std-dev-theta-best 9.750467008132057     :i-best 1
;;; :theta-next [-1.6672422336750374 1.263719010307455]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 868.4782995164372
;;; :log-Z-i-best 1018.1389974817869
;;; :theta-mean-best ([-2.433789402265904 1.2645978847754078] 1015.1380098995651)
;;; :BO-Iteration 6
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.02960899792986642]
;;; :theta-best [-2.433789402265904 1.2645978847754078]     :log-Z-theta-best 1018.1389974817869     :mean-theta-best 994.7998085235489     :std-dev-theta-best 51.03109121549199     :i-best 2
;;; :theta-next [-2.5269982419807553 1.4100774673516547]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 629.2429744894893
;;; :log-Z-i-best 1018.1389974817869
;;; :theta-mean-best ([-2.433789402265904 1.2645978847754078] 994.7998085235489)
;;; :BO-Iteration 7
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.290096494286079]
;;; :theta-best [-2.433789402265904 1.2645978847754078]     :log-Z-theta-best 1018.1389974817869     :mean-theta-best 978.9011851305912     :std-dev-theta-best 90.52386124240358     :i-best 3
;;; :theta-next [-2.357628089105889 1.1442202912372446]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 1123.3180096257458
;;; :log-Z-i-best 1018.1389974817869
;;; :theta-mean-best ([-2.433789402265904 1.2645978847754078] 978.9011851305912)
;;; :BO-Iteration 8
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.01893633975476174]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1112.509177976276     :std-dev-theta-best 19.51569716177706     :i-best 0
;;; :theta-next [-1.4558916723111364 2.9747724651493055]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 263.9115271139044
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1112.509177976276)
;;; :BO-Iteration 9
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.010424877834101665]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1084.1282190506572     :std-dev-theta-best 70.7754399886195     :i-best 1
;;; :theta-next [2.98689993866214 0.38952809500081365]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 637.450233357434
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1084.1282190506572)
;;; :BO-Iteration 10
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.0725919286693711]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1116.2776929681547     :std-dev-theta-best 14.772742523960128     :i-best 2
;;; :theta-next [-1.6068565622223032 1.0950823718162885]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 843.286434485451
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1116.2776929681547)
;;; :BO-Iteration 11
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03022162874027848]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1073.8929930320883     :std-dev-theta-best 61.088641278931505     :i-best 3
;;; :theta-next [2.620213101203042 2.414852473775169]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 282.2745479953092
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1073.8929930320883)
;;; :BO-Iteration 12
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.028976822184186814]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1122.924211287897     :std-dev-theta-best 6.709907755257001     :i-best 4
;;; :theta-next [-2.196635338146864 1.1789520432272247]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 1195.5812326205198
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1122.924211287897)
;;; :BO-Iteration 13
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.012612513556520719]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1195.0882085551802     :std-dev-theta-best 1.3396064831328198     :i-best 0
;;; :theta-next [1.4427750058947098 0.5080585518667765]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -706.4866946583766
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1195.0882085551802)
;;; :BO-Iteration 14
;;; :n-gps-in-acq-function 13
;;; :acq-opt [0.0024596472792707855]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1175.9545809775427     :std-dev-theta-best 32.93095562933427     :i-best 1
;;; :theta-next [-2.9912071253054666 0.36982290236578846]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 854.999416874624
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1175.9545809775427)
;;; :BO-Iteration 15
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07744124264513899]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1086.4036970344182     :std-dev-theta-best 194.32872650750065     :i-best 2
;;; :theta-next [-2.0961166073064303 1.226345574537462]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 984.1835623772664
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1086.4036970344182)
;;; :BO-Iteration 16
;;; :n-gps-in-acq-function 19
;;; :acq-opt [0.13385433441929662]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1150.030170933178     :std-dev-theta-best 43.42003948858712     :i-best 3
;;; :theta-next [-0.6847047930989043 1.7836241948112688]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 700.2832906310888
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1150.030170933178)
;;; :BO-Iteration 17
;;; :n-gps-in-acq-function 14
;;; :acq-opt [0.04960420826468309]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1069.4114864933226     :std-dev-theta-best 183.63364807452572     :i-best 4
;;; :theta-next [-2.542830916263952 2.092251460454982]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 666.6836192653711
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1069.4114864933226)
;;; :BO-Iteration 18
;;; :n-gps-in-acq-function 19
;;; :acq-opt [0.014228358452337431]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1131.843670891314     :std-dev-theta-best 28.571832592507956     :i-best 10
;;; :theta-next [-2.308643302267706 1.2530894742497791]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 936.2056004899861
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1131.843670891314)
;;; :BO-Iteration 19
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.11978126186012125]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1113.5262810474744     :std-dev-theta-best 166.79786652619902     :i-best 6
;;; :theta-next [-1.9642685387885324 1.07508192038059]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 750.3650335909489
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1113.5262810474744)
;;; :BO-Iteration 20
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.10653081803285282]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1157.2629410298377     :std-dev-theta-best 39.49505213396679     :i-best 7
;;; :theta-next [-0.9299855340465037 1.0792650301795976]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 722.9912531042947
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1157.2629410298377)
;;; :BO-Iteration 21
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.026083908768159563]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1089.713999248927     :std-dev-theta-best 77.50466058865067     :i-best 8
;;; :theta-next [2.542292012548896 1.6363475356071662]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 284.0848681242112
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1089.713999248927)
;;; :BO-Iteration 22
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.012677834197953588]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1161.7692682809982     :std-dev-theta-best 44.82975122815257     :i-best 9
;;; :theta-next [1.3144061070897308 2.9353328451569913]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 206.43777434108003
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1161.7692682809982)
;;; :BO-Iteration 23
;;; :n-gps-in-acq-function 18
;;; :acq-opt [0.06034594332295029]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1167.42051479134     :std-dev-theta-best 46.49925198262183     :i-best 10
;;; :theta-next [-2.582849036316664 1.0184652470763236]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 932.9856147421202
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1167.42051479134)
;;; :BO-Iteration 24
;;; :n-gps-in-acq-function 18
;;; :acq-opt [0.07495257026020467]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1090.4150084032858     :std-dev-theta-best 66.69614469223555     :i-best 11
;;; :theta-next [-2.7736348676191502 0.8474065529918291]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 842.0482016663791
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1090.4150084032858)
;;; :BO-Iteration 25
;;; :n-gps-in-acq-function 16
;;; :acq-opt [0.031039498901780958]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1194.9900131474706     :std-dev-theta-best 1.1758129409948939     :i-best 12
;;; :theta-next [-2.769160459835213 2.6624542742243573]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 96.5271898447732
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1194.9900131474706)
;;; :BO-Iteration 26
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.05007211459816568]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1193.7172658425723     :std-dev-theta-best 3.522297091798287     :i-best 13
;;; :theta-next [-1.2468616287984018 1.5817733490000527]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 728.9474852157834
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1193.7172658425723)
;;; :BO-Iteration 27
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.011800892893411944]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1173.4529781180836     :std-dev-theta-best 36.9802292351302     :i-best 14
;;; :theta-next [-0.13778105533336538 0.23300310505820895]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 883.0378493035488
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1173.4529781180836)
;;; :BO-Iteration 28
;;; :n-gps-in-acq-function 18
;;; :acq-opt [0.02523655903706328]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1100.2339470967368     :std-dev-theta-best 96.17322494806999     :i-best 15
;;; :theta-next [-0.20054069899222693 1.4741664298875323]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 776.7678421768736
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1100.2339470967368)
;;; :BO-Iteration 29
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.024115744054361954]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1193.1216670099175     :std-dev-theta-best 5.957697521945133     :i-best 16
;;; :theta-next [0.0559418527089286 2.9054890646841898]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 334.0043755011666
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1193.1216670099175)
;;; :BO-Iteration 30
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.0037784880372047506]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1127.1474711196302     :std-dev-theta-best 15.959245506783219     :i-best 22
;;; :theta-next [-1.9938366468828834 2.2036358367457147]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 634.9958991398499
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1127.1474711196302)
;;; :BO-Iteration 31
;;; :n-gps-in-acq-function 15
;;; :acq-opt [0.006332420933603692]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1173.8715322951666     :std-dev-theta-best 22.409810006955418     :i-best 18
;;; :theta-next [-0.036214537159720006 0.787365585672201]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 879.3798190869186
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1173.8715322951666)
;;; :BO-Iteration 32
;;; :n-gps-in-acq-function 14
;;; :acq-opt [0.02000238665191273]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1136.1718982084485     :std-dev-theta-best 66.38193623571435     :i-best 19
;;; :theta-next [2.912947354602972 0.05519422684167984]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 742.02996606258
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1136.1718982084485)
;;; :BO-Iteration 33
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.011497178625088696]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1123.0944697873465     :std-dev-theta-best 76.49456055234258     :i-best 20
;;; :theta-next [0.44248075915532636 1.1852185672035889]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 789.7517873453676
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1123.0944697873465)
;;; :BO-Iteration 34
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.005389948252437053]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1149.9671615500852     :std-dev-theta-best 52.60423606423803     :i-best 21
;;; :theta-next [-1.4062783176114744 0.7264148137352551]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 773.2249779472743
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1149.9671615500852)
;;; :BO-Iteration 35
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.020606908549764873]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1154.1648031129862     :std-dev-theta-best 63.79181863635301     :i-best 22
;;; :theta-next [2.9303713089554826 1.0420981570076178]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 566.7392529520166
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1154.1648031129862)
;;; :BO-Iteration 36
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.008890358141483912]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1161.592346052346     :std-dev-theta-best 37.61247029591826     :i-best 23
;;; :theta-next [-0.3191688376851949 0.5611138337778115]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 942.641580252025
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1161.592346052346)
;;; :BO-Iteration 37
;;; :n-gps-in-acq-function 13
;;; :acq-opt [0.006105078740548625]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1194.9199051743522     :std-dev-theta-best 1.4973220140298806     :i-best 24
;;; :theta-next [-2.8357718079502465 0.013748168801482037]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 803.6783417777983
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1194.9199051743522)
;;; :BO-Iteration 38
;;; :n-gps-in-acq-function 15
;;; :acq-opt [0.002884706232549153]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1161.911435236798     :std-dev-theta-best 57.333955836160044     :i-best 25
;;; :theta-next [0.27050346914483425 0.04737892942015393]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 730.6924825192917
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1161.911435236798)
;;; :BO-Iteration 39
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07243221702756787]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1149.5970583372273     :std-dev-theta-best 103.41320970361592     :i-best 26
;;; :theta-next [-2.256242297029935 1.1264072777364986]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 760.8536800215413
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1149.5970583372273)
;;; :BO-Iteration 40
;;; :n-gps-in-acq-function 13
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1187.970129764361     :std-dev-theta-best 13.183660092790385     :i-best 27
;;; :theta-next [1.1336860666894317 1.309423082694279]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -99.58586185595273
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1187.970129764361)
;;; :BO-Iteration 41
;;; :n-gps-in-acq-function 13
;;; :acq-opt [0.6431415311117799]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1142.0970205911688     :std-dev-theta-best 93.03118635404026     :i-best 28
;;; :theta-next [-0.9154604199167142 2.9467350451267036]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 163.88390902022093
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1142.0970205911688)
;;; :BO-Iteration 42
;;; :n-gps-in-acq-function 20
;;; :acq-opt [51136.64604174632]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1142.466734879957     :std-dev-theta-best 92.20581329062037     :i-best 29
;;; :theta-next [2.3554686690697952 1.250777476431731]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 435.1609844341028
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1142.466734879957)
;;; :BO-Iteration 43
;;; :n-gps-in-acq-function 17
;;; :acq-opt [6.502184718512822E31]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1085.0028013081928     :std-dev-theta-best 113.01224270722197     :i-best 30
;;; :theta-next [2.4669342139474075 0.7674616726082877]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 317.1361015489083
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1085.0028013081928)
;;; :BO-Iteration 44
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.24061282276179433]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1156.226877090919     :std-dev-theta-best 81.572521929459     :i-best 31
;;; :theta-next [-1.3352772621495657 0.917881327831281]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 587.390960263731
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1156.226877090919)
;;; :BO-Iteration 45
;;; :n-gps-in-acq-function 12
;;; :acq-opt [38181.61846182756]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1139.5608873798772     :std-dev-theta-best 97.4422148739712     :i-best 32
;;; :theta-next [1.479686860406117 2.8082556624439077]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 103.35712317165925
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1139.5608873798772)
;;; :BO-Iteration 46
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03031944929576789]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1114.8120929790693     :std-dev-theta-best 92.4518049007638     :i-best 33
;;; :theta-next [-2.2074360614902 1.931630081969296]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 538.7848874534095
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1114.8120929790693)
;;; :BO-Iteration 47
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.049785730113665085]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 961.0869715552815     :std-dev-theta-best 109.87780742820685     :i-best 39
;;; :theta-next [-2.510269896070101 1.1967340694124387]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 1033.8776790149532
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 961.0869715552815)
;;; :BO-Iteration 48
;;; :n-gps-in-acq-function 20
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1158.9425885258377     :std-dev-theta-best 43.29789788248903     :i-best 35
;;; :theta-next [2.724065970196797 0.8256362105628905]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 487.54800493727697
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1158.9425885258377)
;;; :BO-Iteration 49
;;; :n-gps-in-acq-function 20
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1091.997988690233     :std-dev-theta-best 143.25647531286344     :i-best 36
;;; :theta-next [1.9095755555104814 2.4016225218978864]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 125.20815103100479
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1091.997988690233)
;;; :BO-Iteration 50
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.09467841872019707]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1166.9155410487715     :std-dev-theta-best 45.539334424070724     :i-best 37
;;; :theta-next [-1.2638558702007108 0.5026605419845955]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 673.916026122377
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1166.9155410487715)
;;; :BO-Iteration 51
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.170648441859777]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1125.6951431401426     :std-dev-theta-best 68.31788003000854     :i-best 38
;;; :theta-next [-1.9909279156120592 0.7275380003804355]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 768.4875199213647
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1125.6951431401426)
;;; :BO-Iteration 52
;;; :n-gps-in-acq-function 19
;;; :acq-opt [9.34775638923886E151]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1146.597789145314     :std-dev-theta-best 75.95397373060854     :i-best 39
;;; :theta-next [-2.1308940883594065 0.5980718770960732]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 742.4205656090243
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1146.597789145314)
;;; :BO-Iteration 53
;;; :n-gps-in-acq-function 20
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1150.4674052107828     :std-dev-theta-best 13647.03193170898     :i-best 40
;;; :theta-next [-2.1883508484221714 0.4281279221399957]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 735.2451511814779
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1150.4674052107828)
;;; :BO-Iteration 54
;;; :n-gps-in-acq-function 14
;;; :acq-opt [0.0432176865097679]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1162.081012828271     :std-dev-theta-best 58.12321500471758     :i-best 41
;;; :theta-next [-0.3510805285047365 1.2387688439779934]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 958.4529314402226
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1162.081012828271)
;;; :BO-Iteration 55
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03887865070437953]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1082.5193139184516     :std-dev-theta-best 126.02216938196669     :i-best 42
;;; :theta-next [-0.9634907931057661 2.2649661838353214]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 770.1365438069153
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1082.5193139184516)
;;; :BO-Iteration 56
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.749602110819624]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1113.621978677878     :std-dev-theta-best 80.33304732072362     :i-best 43
;;; :theta-next [-2.861180891995353 1.5680332737101892]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 596.2949592865903
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1113.621978677878)
;;; :BO-Iteration 57
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.08772862363496121]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1111.3714766148596     :std-dev-theta-best 153.78971883677556     :i-best 44
;;; :theta-next [-2.4187649848667085 1.0103332559702383]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 960.4758457367618
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1111.3714766148596)
;;; :BO-Iteration 58
;;; :n-gps-in-acq-function 12
;;; :acq-opt [0.45073390026773535]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1195.2310550873572     :std-dev-theta-best 2.262267958279849     :i-best 45
;;; :theta-next [1.6536675250517119 1.9890406151629003]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 205.4508501632224
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1195.2310550873572)
;;; :BO-Iteration 59
;;; :n-gps-in-acq-function 14
;;; :acq-opt [0.18440698522968965]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1189.3583465073161     :std-dev-theta-best 8.613775182884744     :i-best 46
;;; :theta-next [-2.3926287182448274 0.8011036867185442]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 875.7131065708606
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1189.3583465073161)
;;; :BO-Iteration 60
;;; :n-gps-in-acq-function 20
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1141.887297647789     :std-dev-theta-best 84.94681765540749     :i-best 47
;;; :theta-next [0.5224405862275439 1.1161760298593488]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 764.6580344057994
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1141.887297647789)
;;; :BO-Iteration 61
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.28446607594766554]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1060.77861509157     :std-dev-theta-best 87.62833777541086     :i-best 48
;;; :theta-next [-1.1424482926908899 2.010906044349958]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 707.6898853121032
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1060.77861509157)
;;; :BO-Iteration 62
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04166097732781462]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1177.2470719985472     :std-dev-theta-best 27.892657066579634     :i-best 49
;;; :theta-next [-0.7424973908928245 0.20416661730965283]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 318.92358894035715
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1177.2470719985472)
;;; :BO-Iteration 63
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.022882498221166326]
;;; :theta-best [-2.510269896070101 1.1967340694124387]     :log-Z-theta-best 1033.8776790149532     :mean-theta-best 1044.7011861437654     :std-dev-theta-best 13.290090790512785     :i-best 15
;;; :theta-next [0.0019391602363896254 1.034546659333309]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 886.7932528078478
;;; :log-Z-i-best 1033.8776790149532
;;; :theta-mean-best ([-2.510269896070101 1.1967340694124387] 1044.7011861437654)
;;; :BO-Iteration 64
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.7015474525248919]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1115.344707888852     :std-dev-theta-best 63.003934702059134     :i-best 51
;;; :theta-next [-2.106180106193408 1.06163914685329]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 1163.2134804403074
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1115.344707888852)
;;; :BO-Iteration 65
;;; :n-gps-in-acq-function 17
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1147.8987973519854     :std-dev-theta-best 82.49447241415233     :i-best 52
;;; :theta-next [-0.2925819987803586 0.44943633289414303]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 930.7264229846563
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1147.8987973519854)
;;; :BO-Iteration 66
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.14103128148948432]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1080.1269486614974     :std-dev-theta-best 63.73621065663623     :i-best 53
;;; :theta-next [-2.1067016858677956 1.0191155271175374]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 943.6477771173332
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1080.1269486614974)
;;; :BO-Iteration 67
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.013556334127375187]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1071.198599165934     :std-dev-theta-best 61.3934044172399     :i-best 54
;;; :theta-next [-0.3681921770893002 1.0483603688867826]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 963.358793290074
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1071.198599165934)
;;; :BO-Iteration 68
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.029127017683271125]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1031.9454680237463     :std-dev-theta-best 77.47203588989238     :i-best 55
;;; :theta-next [-0.36078081368427073 0.7844409232927622]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 956.6898652908649
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1031.9454680237463)
;;; :BO-Iteration 69
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03249280376583735]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1155.9367609880373     :std-dev-theta-best 64.48724026243234     :i-best 56
;;; :theta-next [-0.7683759284589549 1.3875425429865273]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 872.4523632621862
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1155.9367609880373)
;;; :BO-Iteration 70
;;; :n-gps-in-acq-function 18
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1195.4772614796702     :std-dev-theta-best 5410.1637330983285     :i-best 57
;;; :theta-next [-1.5575847335129365 2.5045957684102547]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 470.47783139837213
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1195.4772614796702)
;;; :BO-Iteration 71
;;; :n-gps-in-acq-function 19
;;; :acq-opt [0.0991689665080797]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1107.669579055972     :std-dev-theta-best 82.17441171384351     :i-best 58
;;; :theta-next [2.829286743811662 1.965321115979254]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 95.14602762573561
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1107.669579055972)
;;; :BO-Iteration 72
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.68657440230881]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1108.8517696865679     :std-dev-theta-best 77.16810442602133     :i-best 59
;;; :theta-next [0.2372314494418042 1.3745582127263118]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 839.0075304990365
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1108.8517696865679)
;;; :BO-Iteration 73
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.06548959104554283]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1082.024347850384     :std-dev-theta-best 66.89200024837807     :i-best 60
;;; :theta-next [-2.1217568802641207 1.143420271609355]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 784.9308698858368
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1082.024347850384)
;;; :BO-Iteration 74
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.5470205396065311]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1146.3152360630456     :std-dev-theta-best 95.14158737101555     :i-best 61
;;; :theta-next [0.2054779052896003 1.5678727049689918]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 656.7506372480588
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1146.3152360630456)
;;; :BO-Iteration 75
;;; :n-gps-in-acq-function 19
;;; :acq-opt [0.11723332326462318]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1138.2004083911534     :std-dev-theta-best 61.89816090233163     :i-best 62
;;; :theta-next [-2.457546529448684 0.11959907266374975]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 621.2642758984189
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1138.2004083911534)
;;; :BO-Iteration 76
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.9269401046612509]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1052.2298280034497     :std-dev-theta-best 108.27715057183254     :i-best 63
;;; :theta-next [2.2813491952110967 1.018792547583525]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 40.676860602121906
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1052.2298280034497)
;;; :BO-Iteration 77
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.006690048807645925]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1025.5550596673174     :std-dev-theta-best 14.400402911940734     :i-best 69
;;; :theta-next [-0.022478675800917003 0.45507765878138606]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 868.7301762258827
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1025.5550596673174)
;;; :BO-Iteration 78
;;; :n-gps-in-acq-function 20
;;; :acq-opt [7.4833948620333261E18]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1143.967585645226     :std-dev-theta-best 84.33374336618787     :i-best 65
;;; :theta-next [-0.33937209155167425 0.8314738759773872]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 957.304747106531
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1143.967585645226)
;;; :BO-Iteration 79
;;; :n-gps-in-acq-function 9
;;; :acq-opt [0.1416433360409727]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1084.4898799113823     :std-dev-theta-best 101.62926993080579     :i-best 66
;;; :theta-next [-1.15372497835562 1.2407041765775104]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 682.9937208552161
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1084.4898799113823)
;;; :BO-Iteration 80
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.10962908281118067]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1008.602689461943     :std-dev-theta-best 127.65810247497893     :i-best 67
;;; :theta-next [2.3626914454444865 0.2047522373450763]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 421.78101240021607
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1008.602689461943)
;;; :BO-Iteration 81
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.0171416318196516E8]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1037.2212043106083     :std-dev-theta-best 43.001278858114404     :i-best 73
;;; :theta-next [2.045420761709508 1.809698830233693]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 309.93144669091544
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1037.2212043106083)
;;; :BO-Iteration 82
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.5628553928419069]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1099.4686633707556     :std-dev-theta-best 101.66880686932494     :i-best 69
;;; :theta-next [-2.23005336472265 1.6530421640014963]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 597.7554641004839
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1099.4686633707556)
;;; :BO-Iteration 83
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.0058620196706511375]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1021.879384364581     :std-dev-theta-best 19.20525791431786     :i-best 75
;;; :theta-next [-0.2978144281604016 2.074036921350249]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 743.4890279106992
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1021.879384364581)
;;; :BO-Iteration 84
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.12969361492682147]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1083.3570800467287     :std-dev-theta-best 89.05161465427386     :i-best 71
;;; :theta-next [2.774021979467335 2.3303954649008523]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 340.775062345005
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1083.3570800467287)
;;; :BO-Iteration 85
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.3217992146905419]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1056.0313256727677     :std-dev-theta-best 357.94177946342137     :i-best 72
;;; :theta-next [-1.9264976077523868 1.6793325780240653]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 854.5157246510294
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1056.0313256727677)
;;; :BO-Iteration 86
;;; :n-gps-in-acq-function 17
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1115.9284862034574     :std-dev-theta-best 672.9567555743554     :i-best 73
;;; :theta-next [-0.0603250962415296 0.5006416723631503]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 879.0060843529577
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1115.9284862034574)
;;; :BO-Iteration 87
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.29290501170283345]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1132.589809456501     :std-dev-theta-best 111.295561041495     :i-best 74
;;; :theta-next [0.44880214060436296 2.8167747492480943]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 576.5175569130088
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1132.589809456501)
;;; :BO-Iteration 88
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.05249883185396363]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 961.981600145527     :std-dev-theta-best 76.98306624293214     :i-best 80
;;; :theta-next [-2.470198756279432 1.1428259022695613]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 942.4812356791647
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 961.981600145527)
;;; :BO-Iteration 89
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.11093377104315881]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1096.9784385121939     :std-dev-theta-best 102.49552320309975     :i-best 76
;;; :theta-next [0.48170935257485503 0.23241406329809983]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 591.6208453291862
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1096.9784385121939)
;;; :BO-Iteration 90
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.4743357477650997E59]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1089.8803805323644     :std-dev-theta-best 109.08493282312818     :i-best 77
;;; :theta-next [1.9135675861786767 1.2270826919721904]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 406.30004709667816
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1089.8803805323644)
;;; :BO-Iteration 91
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.01828161182094796]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1047.87496927762     :std-dev-theta-best 81.35308856174247     :i-best 78
;;; :theta-next [-0.058599401679452434 2.4509222892019107]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 784.3078528527576
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1047.87496927762)
;;; :BO-Iteration 92
;;; :n-gps-in-acq-function 20
;;; :acq-opt [Infinity]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1062.9097202338082     :std-dev-theta-best 112.51895125062133     :i-best 79
;;; :theta-next [1.3496864262708046 0.41647010034472987]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -843.8388927105816
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1062.9097202338082)
;;; :BO-Iteration 93
;;; :n-gps-in-acq-function 10
;;; :acq-opt [4.155577262757159E8]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1139.7788474903587     :std-dev-theta-best 93.22476016669808     :i-best 80
;;; :theta-next [0.25941146693999695 1.6761564415599686]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 598.5629234971093
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1139.7788474903587)
;;; :BO-Iteration 94
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.6626621362325625E9]
;;; :theta-best [-2.357628089105889 1.1442202912372446]     :log-Z-theta-best 1123.3180096257458     :mean-theta-best 1015.5762204417507     :std-dev-theta-best 26.536546366423437     :i-best 86
;;; :theta-next [0.32722464512084404 1.3028826825765125]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 828.465117539396
;;; :log-Z-i-best 1123.3180096257458
;;; :theta-mean-best ([-2.357628089105889 1.1442202912372446] 1015.5762204417507)
;;; :BO-Iteration 95
;;; :n-gps-in-acq-function 19
;;; :acq-opt [0.15123049899089788]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1051.2818227961109     :std-dev-theta-best 85.19417481020798     :i-best 82
;;; :theta-next [2.048258941695126 0.009284736011677593]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 227.29338534785091
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1051.2818227961109)
;;; :BO-Iteration 96
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07787304528811385]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1036.3212780256956     :std-dev-theta-best 113.67636543864513     :i-best 83
;;; :theta-next [0.9863484598060301 0.811305575446579]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -815.6164330136986
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1036.3212780256956)
;;; :BO-Iteration 97
;;; :n-gps-in-acq-function 16
;;; :acq-opt [0.11458176447780159]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1107.3957116729785     :std-dev-theta-best 95.28553254698133     :i-best 84
;;; :theta-next [-0.5103372732505078 2.2186570611312364]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 760.7545219558673
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1107.3957116729785)
;;; :BO-Iteration 98
;;; :n-gps-in-acq-function 14
;;; :acq-opt [0.18005267270018294]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1193.2931684982539     :std-dev-theta-best 7.569379116933334     :i-best 85
;;; :theta-next [-2.849032878931107 2.40190675338458]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 511.61886167722946
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1193.2931684982539)
;;; :BO-Iteration 99
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.3293272107627172]
;;; :theta-best [-2.196635338146864 1.1789520432272247]     :log-Z-theta-best 1195.5812326205198     :mean-theta-best 1079.5756424160159     :std-dev-theta-best 78.02586344017547     :i-best 86
;;; :theta-next [-2.0526535645201616 1.3838981399424035]
;;; Calling original query with theta next  
;;; :log-Z-theta-next 575.76921703097
;;; :log-Z-i-best 1195.5812326205198
;;; :theta-mean-best ([-2.196635338146864 1.1789520432272247] 1079.5756424160159)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.chaos/samples</span>","value":"#'worksheets.chaos/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.8531507208526143</span>","value":"-1.8531507208526143"},{"type":"html","content":"<span class='clj-double'>1.763583387232576</span>","value":"1.763583387232576"}],"value":"[-1.8531507208526143 1.763583387232576]"},{"type":"html","content":"<span class='clj-double'>715.622733900344</span>","value":"715.622733900344"}],"value":"([-1.8531507208526143 1.763583387232576] 715.622733900344)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.398537226969978</span>","value":"-2.398537226969978"},{"type":"html","content":"<span class='clj-double'>0.566750291398498</span>","value":"0.566750291398498"}],"value":"[-2.398537226969978 0.566750291398498]"},{"type":"html","content":"<span class='clj-double'>824.6971612982984</span>","value":"824.6971612982984"}],"value":"([-2.398537226969978 0.566750291398498] 824.6971612982984)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.398537226969978</span>","value":"-2.398537226969978"},{"type":"html","content":"<span class='clj-double'>0.566750291398498</span>","value":"0.566750291398498"}],"value":"[-2.398537226969978 0.566750291398498]"},{"type":"html","content":"<span class='clj-double'>826.5795790501979</span>","value":"826.5795790501979"}],"value":"([-2.398537226969978 0.566750291398498] 826.5795790501979)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.217440843085523</span>","value":"-2.217440843085523"},{"type":"html","content":"<span class='clj-double'>0.939084122699758</span>","value":"0.939084122699758"}],"value":"[-2.217440843085523 0.939084122699758]"},{"type":"html","content":"<span class='clj-double'>856.0665613110718</span>","value":"856.0665613110718"}],"value":"([-2.217440843085523 0.939084122699758] 856.0665613110718)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.433789402265904</span>","value":"-2.433789402265904"},{"type":"html","content":"<span class='clj-double'>1.2645978847754078</span>","value":"1.2645978847754078"}],"value":"[-2.433789402265904 1.2645978847754078]"},{"type":"html","content":"<span class='clj-double'>920.9251877697038</span>","value":"920.9251877697038"}],"value":"([-2.433789402265904 1.2645978847754078] 920.9251877697038)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.433789402265904</span>","value":"-2.433789402265904"},{"type":"html","content":"<span class='clj-double'>1.2645978847754078</span>","value":"1.2645978847754078"}],"value":"[-2.433789402265904 1.2645978847754078]"},{"type":"html","content":"<span class='clj-double'>1015.1380098995651</span>","value":"1015.1380098995651"}],"value":"([-2.433789402265904 1.2645978847754078] 1015.1380098995651)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.433789402265904</span>","value":"-2.433789402265904"},{"type":"html","content":"<span class='clj-double'>1.2645978847754078</span>","value":"1.2645978847754078"}],"value":"[-2.433789402265904 1.2645978847754078]"},{"type":"html","content":"<span class='clj-double'>994.7998085235489</span>","value":"994.7998085235489"}],"value":"([-2.433789402265904 1.2645978847754078] 994.7998085235489)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.433789402265904</span>","value":"-2.433789402265904"},{"type":"html","content":"<span class='clj-double'>1.2645978847754078</span>","value":"1.2645978847754078"}],"value":"[-2.433789402265904 1.2645978847754078]"},{"type":"html","content":"<span class='clj-double'>978.9011851305912</span>","value":"978.9011851305912"}],"value":"([-2.433789402265904 1.2645978847754078] 978.9011851305912)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1112.509177976276</span>","value":"1112.509177976276"}],"value":"([-2.357628089105889 1.1442202912372446] 1112.509177976276)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1084.1282190506572</span>","value":"1084.1282190506572"}],"value":"([-2.357628089105889 1.1442202912372446] 1084.1282190506572)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1116.2776929681547</span>","value":"1116.2776929681547"}],"value":"([-2.357628089105889 1.1442202912372446] 1116.2776929681547)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1073.8929930320883</span>","value":"1073.8929930320883"}],"value":"([-2.357628089105889 1.1442202912372446] 1073.8929930320883)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1122.924211287897</span>","value":"1122.924211287897"}],"value":"([-2.357628089105889 1.1442202912372446] 1122.924211287897)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1195.0882085551802</span>","value":"1195.0882085551802"}],"value":"([-2.196635338146864 1.1789520432272247] 1195.0882085551802)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1175.9545809775427</span>","value":"1175.9545809775427"}],"value":"([-2.196635338146864 1.1789520432272247] 1175.9545809775427)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1086.4036970344182</span>","value":"1086.4036970344182"}],"value":"([-2.196635338146864 1.1789520432272247] 1086.4036970344182)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1150.030170933178</span>","value":"1150.030170933178"}],"value":"([-2.196635338146864 1.1789520432272247] 1150.030170933178)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1069.4114864933226</span>","value":"1069.4114864933226"}],"value":"([-2.196635338146864 1.1789520432272247] 1069.4114864933226)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1131.843670891314</span>","value":"1131.843670891314"}],"value":"([-2.357628089105889 1.1442202912372446] 1131.843670891314)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1113.5262810474744</span>","value":"1113.5262810474744"}],"value":"([-2.196635338146864 1.1789520432272247] 1113.5262810474744)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1157.2629410298377</span>","value":"1157.2629410298377"}],"value":"([-2.196635338146864 1.1789520432272247] 1157.2629410298377)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1089.713999248927</span>","value":"1089.713999248927"}],"value":"([-2.196635338146864 1.1789520432272247] 1089.713999248927)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1161.7692682809982</span>","value":"1161.7692682809982"}],"value":"([-2.196635338146864 1.1789520432272247] 1161.7692682809982)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1167.42051479134</span>","value":"1167.42051479134"}],"value":"([-2.196635338146864 1.1789520432272247] 1167.42051479134)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1090.4150084032858</span>","value":"1090.4150084032858"}],"value":"([-2.196635338146864 1.1789520432272247] 1090.4150084032858)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1194.9900131474706</span>","value":"1194.9900131474706"}],"value":"([-2.196635338146864 1.1789520432272247] 1194.9900131474706)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1193.7172658425723</span>","value":"1193.7172658425723"}],"value":"([-2.196635338146864 1.1789520432272247] 1193.7172658425723)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1173.4529781180836</span>","value":"1173.4529781180836"}],"value":"([-2.196635338146864 1.1789520432272247] 1173.4529781180836)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1100.2339470967368</span>","value":"1100.2339470967368"}],"value":"([-2.196635338146864 1.1789520432272247] 1100.2339470967368)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1193.1216670099175</span>","value":"1193.1216670099175"}],"value":"([-2.196635338146864 1.1789520432272247] 1193.1216670099175)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1127.1474711196302</span>","value":"1127.1474711196302"}],"value":"([-2.357628089105889 1.1442202912372446] 1127.1474711196302)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1173.8715322951666</span>","value":"1173.8715322951666"}],"value":"([-2.196635338146864 1.1789520432272247] 1173.8715322951666)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1136.1718982084485</span>","value":"1136.1718982084485"}],"value":"([-2.196635338146864 1.1789520432272247] 1136.1718982084485)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1123.0944697873465</span>","value":"1123.0944697873465"}],"value":"([-2.196635338146864 1.1789520432272247] 1123.0944697873465)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1149.9671615500852</span>","value":"1149.9671615500852"}],"value":"([-2.196635338146864 1.1789520432272247] 1149.9671615500852)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1154.1648031129862</span>","value":"1154.1648031129862"}],"value":"([-2.196635338146864 1.1789520432272247] 1154.1648031129862)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1161.592346052346</span>","value":"1161.592346052346"}],"value":"([-2.196635338146864 1.1789520432272247] 1161.592346052346)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1194.9199051743522</span>","value":"1194.9199051743522"}],"value":"([-2.196635338146864 1.1789520432272247] 1194.9199051743522)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1161.911435236798</span>","value":"1161.911435236798"}],"value":"([-2.196635338146864 1.1789520432272247] 1161.911435236798)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1149.5970583372273</span>","value":"1149.5970583372273"}],"value":"([-2.196635338146864 1.1789520432272247] 1149.5970583372273)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1187.970129764361</span>","value":"1187.970129764361"}],"value":"([-2.196635338146864 1.1789520432272247] 1187.970129764361)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1142.0970205911688</span>","value":"1142.0970205911688"}],"value":"([-2.196635338146864 1.1789520432272247] 1142.0970205911688)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1142.466734879957</span>","value":"1142.466734879957"}],"value":"([-2.196635338146864 1.1789520432272247] 1142.466734879957)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1085.0028013081928</span>","value":"1085.0028013081928"}],"value":"([-2.196635338146864 1.1789520432272247] 1085.0028013081928)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1156.226877090919</span>","value":"1156.226877090919"}],"value":"([-2.196635338146864 1.1789520432272247] 1156.226877090919)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1139.5608873798772</span>","value":"1139.5608873798772"}],"value":"([-2.196635338146864 1.1789520432272247] 1139.5608873798772)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1114.8120929790693</span>","value":"1114.8120929790693"}],"value":"([-2.196635338146864 1.1789520432272247] 1114.8120929790693)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>961.0869715552815</span>","value":"961.0869715552815"}],"value":"([-2.357628089105889 1.1442202912372446] 961.0869715552815)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1158.9425885258377</span>","value":"1158.9425885258377"}],"value":"([-2.196635338146864 1.1789520432272247] 1158.9425885258377)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1091.997988690233</span>","value":"1091.997988690233"}],"value":"([-2.196635338146864 1.1789520432272247] 1091.997988690233)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1166.9155410487715</span>","value":"1166.9155410487715"}],"value":"([-2.196635338146864 1.1789520432272247] 1166.9155410487715)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1125.6951431401426</span>","value":"1125.6951431401426"}],"value":"([-2.196635338146864 1.1789520432272247] 1125.6951431401426)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1146.597789145314</span>","value":"1146.597789145314"}],"value":"([-2.196635338146864 1.1789520432272247] 1146.597789145314)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1150.4674052107828</span>","value":"1150.4674052107828"}],"value":"([-2.196635338146864 1.1789520432272247] 1150.4674052107828)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1162.081012828271</span>","value":"1162.081012828271"}],"value":"([-2.196635338146864 1.1789520432272247] 1162.081012828271)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1082.5193139184516</span>","value":"1082.5193139184516"}],"value":"([-2.196635338146864 1.1789520432272247] 1082.5193139184516)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1113.621978677878</span>","value":"1113.621978677878"}],"value":"([-2.196635338146864 1.1789520432272247] 1113.621978677878)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1111.3714766148596</span>","value":"1111.3714766148596"}],"value":"([-2.196635338146864 1.1789520432272247] 1111.3714766148596)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1195.2310550873572</span>","value":"1195.2310550873572"}],"value":"([-2.196635338146864 1.1789520432272247] 1195.2310550873572)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1189.3583465073161</span>","value":"1189.3583465073161"}],"value":"([-2.196635338146864 1.1789520432272247] 1189.3583465073161)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1141.887297647789</span>","value":"1141.887297647789"}],"value":"([-2.196635338146864 1.1789520432272247] 1141.887297647789)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1060.77861509157</span>","value":"1060.77861509157"}],"value":"([-2.196635338146864 1.1789520432272247] 1060.77861509157)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1177.2470719985472</span>","value":"1177.2470719985472"}],"value":"([-2.196635338146864 1.1789520432272247] 1177.2470719985472)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.510269896070101</span>","value":"-2.510269896070101"},{"type":"html","content":"<span class='clj-double'>1.1967340694124387</span>","value":"1.1967340694124387"}],"value":"[-2.510269896070101 1.1967340694124387]"},{"type":"html","content":"<span class='clj-double'>1044.7011861437654</span>","value":"1044.7011861437654"}],"value":"([-2.510269896070101 1.1967340694124387] 1044.7011861437654)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1115.344707888852</span>","value":"1115.344707888852"}],"value":"([-2.196635338146864 1.1789520432272247] 1115.344707888852)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1147.8987973519854</span>","value":"1147.8987973519854"}],"value":"([-2.196635338146864 1.1789520432272247] 1147.8987973519854)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1080.1269486614974</span>","value":"1080.1269486614974"}],"value":"([-2.196635338146864 1.1789520432272247] 1080.1269486614974)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1071.198599165934</span>","value":"1071.198599165934"}],"value":"([-2.196635338146864 1.1789520432272247] 1071.198599165934)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1031.9454680237463</span>","value":"1031.9454680237463"}],"value":"([-2.196635338146864 1.1789520432272247] 1031.9454680237463)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1155.9367609880373</span>","value":"1155.9367609880373"}],"value":"([-2.196635338146864 1.1789520432272247] 1155.9367609880373)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1195.4772614796702</span>","value":"1195.4772614796702"}],"value":"([-2.196635338146864 1.1789520432272247] 1195.4772614796702)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1107.669579055972</span>","value":"1107.669579055972"}],"value":"([-2.196635338146864 1.1789520432272247] 1107.669579055972)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1108.8517696865679</span>","value":"1108.8517696865679"}],"value":"([-2.196635338146864 1.1789520432272247] 1108.8517696865679)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1082.024347850384</span>","value":"1082.024347850384"}],"value":"([-2.196635338146864 1.1789520432272247] 1082.024347850384)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1146.3152360630456</span>","value":"1146.3152360630456"}],"value":"([-2.196635338146864 1.1789520432272247] 1146.3152360630456)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1138.2004083911534</span>","value":"1138.2004083911534"}],"value":"([-2.196635338146864 1.1789520432272247] 1138.2004083911534)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1052.2298280034497</span>","value":"1052.2298280034497"}],"value":"([-2.196635338146864 1.1789520432272247] 1052.2298280034497)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1025.5550596673174</span>","value":"1025.5550596673174"}],"value":"([-2.357628089105889 1.1442202912372446] 1025.5550596673174)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1143.967585645226</span>","value":"1143.967585645226"}],"value":"([-2.196635338146864 1.1789520432272247] 1143.967585645226)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1084.4898799113823</span>","value":"1084.4898799113823"}],"value":"([-2.196635338146864 1.1789520432272247] 1084.4898799113823)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1008.602689461943</span>","value":"1008.602689461943"}],"value":"([-2.196635338146864 1.1789520432272247] 1008.602689461943)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1037.2212043106083</span>","value":"1037.2212043106083"}],"value":"([-2.357628089105889 1.1442202912372446] 1037.2212043106083)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1099.4686633707556</span>","value":"1099.4686633707556"}],"value":"([-2.196635338146864 1.1789520432272247] 1099.4686633707556)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1021.879384364581</span>","value":"1021.879384364581"}],"value":"([-2.357628089105889 1.1442202912372446] 1021.879384364581)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1083.3570800467287</span>","value":"1083.3570800467287"}],"value":"([-2.196635338146864 1.1789520432272247] 1083.3570800467287)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1056.0313256727677</span>","value":"1056.0313256727677"}],"value":"([-2.196635338146864 1.1789520432272247] 1056.0313256727677)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1115.9284862034574</span>","value":"1115.9284862034574"}],"value":"([-2.196635338146864 1.1789520432272247] 1115.9284862034574)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1132.589809456501</span>","value":"1132.589809456501"}],"value":"([-2.196635338146864 1.1789520432272247] 1132.589809456501)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>961.981600145527</span>","value":"961.981600145527"}],"value":"([-2.357628089105889 1.1442202912372446] 961.981600145527)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1096.9784385121939</span>","value":"1096.9784385121939"}],"value":"([-2.196635338146864 1.1789520432272247] 1096.9784385121939)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1089.8803805323644</span>","value":"1089.8803805323644"}],"value":"([-2.196635338146864 1.1789520432272247] 1089.8803805323644)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1047.87496927762</span>","value":"1047.87496927762"}],"value":"([-2.196635338146864 1.1789520432272247] 1047.87496927762)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1062.9097202338082</span>","value":"1062.9097202338082"}],"value":"([-2.196635338146864 1.1789520432272247] 1062.9097202338082)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1139.7788474903587</span>","value":"1139.7788474903587"}],"value":"([-2.196635338146864 1.1789520432272247] 1139.7788474903587)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.357628089105889</span>","value":"-2.357628089105889"},{"type":"html","content":"<span class='clj-double'>1.1442202912372446</span>","value":"1.1442202912372446"}],"value":"[-2.357628089105889 1.1442202912372446]"},{"type":"html","content":"<span class='clj-double'>1015.5762204417507</span>","value":"1015.5762204417507"}],"value":"([-2.357628089105889 1.1442202912372446] 1015.5762204417507)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1051.2818227961109</span>","value":"1051.2818227961109"}],"value":"([-2.196635338146864 1.1789520432272247] 1051.2818227961109)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1036.3212780256956</span>","value":"1036.3212780256956"}],"value":"([-2.196635338146864 1.1789520432272247] 1036.3212780256956)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1107.3957116729785</span>","value":"1107.3957116729785"}],"value":"([-2.196635338146864 1.1789520432272247] 1107.3957116729785)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1193.2931684982539</span>","value":"1193.2931684982539"}],"value":"([-2.196635338146864 1.1789520432272247] 1193.2931684982539)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-2.196635338146864</span>","value":"-2.196635338146864"},{"type":"html","content":"<span class='clj-double'>1.1789520432272247</span>","value":"1.1789520432272247"}],"value":"[-2.196635338146864 1.1789520432272247]"},{"type":"html","content":"<span class='clj-double'>1079.5756424160159</span>","value":"1079.5756424160159"}],"value":"([-2.196635338146864 1.1789520432272247] 1079.5756424160159)"}],"value":"[([-1.8531507208526143 1.763583387232576] 715.622733900344) ([-2.398537226969978 0.566750291398498] 824.6971612982984) ([-2.398537226969978 0.566750291398498] 826.5795790501979) ([-2.217440843085523 0.939084122699758] 856.0665613110718) ([-2.433789402265904 1.2645978847754078] 920.9251877697038) ([-2.433789402265904 1.2645978847754078] 1015.1380098995651) ([-2.433789402265904 1.2645978847754078] 994.7998085235489) ([-2.433789402265904 1.2645978847754078] 978.9011851305912) ([-2.357628089105889 1.1442202912372446] 1112.509177976276) ([-2.357628089105889 1.1442202912372446] 1084.1282190506572) ([-2.357628089105889 1.1442202912372446] 1116.2776929681547) ([-2.357628089105889 1.1442202912372446] 1073.8929930320883) ([-2.357628089105889 1.1442202912372446] 1122.924211287897) ([-2.196635338146864 1.1789520432272247] 1195.0882085551802) ([-2.196635338146864 1.1789520432272247] 1175.9545809775427) ([-2.196635338146864 1.1789520432272247] 1086.4036970344182) ([-2.196635338146864 1.1789520432272247] 1150.030170933178) ([-2.196635338146864 1.1789520432272247] 1069.4114864933226) ([-2.357628089105889 1.1442202912372446] 1131.843670891314) ([-2.196635338146864 1.1789520432272247] 1113.5262810474744) ([-2.196635338146864 1.1789520432272247] 1157.2629410298377) ([-2.196635338146864 1.1789520432272247] 1089.713999248927) ([-2.196635338146864 1.1789520432272247] 1161.7692682809982) ([-2.196635338146864 1.1789520432272247] 1167.42051479134) ([-2.196635338146864 1.1789520432272247] 1090.4150084032858) ([-2.196635338146864 1.1789520432272247] 1194.9900131474706) ([-2.196635338146864 1.1789520432272247] 1193.7172658425723) ([-2.196635338146864 1.1789520432272247] 1173.4529781180836) ([-2.196635338146864 1.1789520432272247] 1100.2339470967368) ([-2.196635338146864 1.1789520432272247] 1193.1216670099175) ([-2.357628089105889 1.1442202912372446] 1127.1474711196302) ([-2.196635338146864 1.1789520432272247] 1173.8715322951666) ([-2.196635338146864 1.1789520432272247] 1136.1718982084485) ([-2.196635338146864 1.1789520432272247] 1123.0944697873465) ([-2.196635338146864 1.1789520432272247] 1149.9671615500852) ([-2.196635338146864 1.1789520432272247] 1154.1648031129862) ([-2.196635338146864 1.1789520432272247] 1161.592346052346) ([-2.196635338146864 1.1789520432272247] 1194.9199051743522) ([-2.196635338146864 1.1789520432272247] 1161.911435236798) ([-2.196635338146864 1.1789520432272247] 1149.5970583372273) ([-2.196635338146864 1.1789520432272247] 1187.970129764361) ([-2.196635338146864 1.1789520432272247] 1142.0970205911688) ([-2.196635338146864 1.1789520432272247] 1142.466734879957) ([-2.196635338146864 1.1789520432272247] 1085.0028013081928) ([-2.196635338146864 1.1789520432272247] 1156.226877090919) ([-2.196635338146864 1.1789520432272247] 1139.5608873798772) ([-2.196635338146864 1.1789520432272247] 1114.8120929790693) ([-2.357628089105889 1.1442202912372446] 961.0869715552815) ([-2.196635338146864 1.1789520432272247] 1158.9425885258377) ([-2.196635338146864 1.1789520432272247] 1091.997988690233) ([-2.196635338146864 1.1789520432272247] 1166.9155410487715) ([-2.196635338146864 1.1789520432272247] 1125.6951431401426) ([-2.196635338146864 1.1789520432272247] 1146.597789145314) ([-2.196635338146864 1.1789520432272247] 1150.4674052107828) ([-2.196635338146864 1.1789520432272247] 1162.081012828271) ([-2.196635338146864 1.1789520432272247] 1082.5193139184516) ([-2.196635338146864 1.1789520432272247] 1113.621978677878) ([-2.196635338146864 1.1789520432272247] 1111.3714766148596) ([-2.196635338146864 1.1789520432272247] 1195.2310550873572) ([-2.196635338146864 1.1789520432272247] 1189.3583465073161) ([-2.196635338146864 1.1789520432272247] 1141.887297647789) ([-2.196635338146864 1.1789520432272247] 1060.77861509157) ([-2.196635338146864 1.1789520432272247] 1177.2470719985472) ([-2.510269896070101 1.1967340694124387] 1044.7011861437654) ([-2.196635338146864 1.1789520432272247] 1115.344707888852) ([-2.196635338146864 1.1789520432272247] 1147.8987973519854) ([-2.196635338146864 1.1789520432272247] 1080.1269486614974) ([-2.196635338146864 1.1789520432272247] 1071.198599165934) ([-2.196635338146864 1.1789520432272247] 1031.9454680237463) ([-2.196635338146864 1.1789520432272247] 1155.9367609880373) ([-2.196635338146864 1.1789520432272247] 1195.4772614796702) ([-2.196635338146864 1.1789520432272247] 1107.669579055972) ([-2.196635338146864 1.1789520432272247] 1108.8517696865679) ([-2.196635338146864 1.1789520432272247] 1082.024347850384) ([-2.196635338146864 1.1789520432272247] 1146.3152360630456) ([-2.196635338146864 1.1789520432272247] 1138.2004083911534) ([-2.196635338146864 1.1789520432272247] 1052.2298280034497) ([-2.357628089105889 1.1442202912372446] 1025.5550596673174) ([-2.196635338146864 1.1789520432272247] 1143.967585645226) ([-2.196635338146864 1.1789520432272247] 1084.4898799113823) ([-2.196635338146864 1.1789520432272247] 1008.602689461943) ([-2.357628089105889 1.1442202912372446] 1037.2212043106083) ([-2.196635338146864 1.1789520432272247] 1099.4686633707556) ([-2.357628089105889 1.1442202912372446] 1021.879384364581) ([-2.196635338146864 1.1789520432272247] 1083.3570800467287) ([-2.196635338146864 1.1789520432272247] 1056.0313256727677) ([-2.196635338146864 1.1789520432272247] 1115.9284862034574) ([-2.196635338146864 1.1789520432272247] 1132.589809456501) ([-2.357628089105889 1.1442202912372446] 961.981600145527) ([-2.196635338146864 1.1789520432272247] 1096.9784385121939) ([-2.196635338146864 1.1789520432272247] 1089.8803805323644) ([-2.196635338146864 1.1789520432272247] 1047.87496927762) ([-2.196635338146864 1.1789520432272247] 1062.9097202338082) ([-2.196635338146864 1.1789520432272247] 1139.7788474903587) ([-2.357628089105889 1.1442202912372446] 1015.5762204417507) ([-2.196635338146864 1.1789520432272247] 1051.2818227961109) ([-2.196635338146864 1.1789520432272247] 1036.3212780256956) ([-2.196635338146864 1.1789520432272247] 1107.3957116729785) ([-2.196635338146864 1.1789520432272247] 1193.2931684982539) ([-2.196635338146864 1.1789520432272247] 1079.5756424160159)]"}
;; <=

;; @@

;; @@
