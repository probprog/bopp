(ns examples.chaos
    (:require [anglican.core :refer [doquery]]
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


(defn square
  [x]
  (m/mul x x))


(defn chaos-data [n-data n-points]
  (->> (str "data/chaos/y" n-data ".csv")
       io/resource
       io/reader
       slurp
       csv/read-csv
       (into [])
       (mapv #(mapv read-string %))
       (take n-points)
       (into [])))

(defn transition-dist-general-clj
  [b c beta nu x y z]
  "Chaotic transition according to the pickover attractor"
  [(- (sin (* beta y)) (* z (cos (* b x))))
   (- (* z (sin (* c x))) (cos (* nu y)))
   (sin x)])

(def b 5/2)
(def c -3/2)

(def transition-dist-set-clj (partial transition-dist-general-clj b c))

(def obs-matrix
  "Observation matrix"
  (transpose (matrix
   [[0.0243087960714491	0.0168113871672383	0.0691772445397626	1.57387136968515e-05	0.303732278352682	0.0618563110870659	4.95623561031004e-06	1.89039711478951e-08	0.00113995322996380	0.00515810614718469	0.00341658550834822	7.85881265878172e-07	1.30076659987469e-16	0.291018005156189	3.67275441629293e-10	3.20483237445295e-05	1.63410289930715e-18	1.65003111575406e-05	0.211775125663058	0.0115361583403371]
    [2.03156858097514e-14	1.61789389343840e-08	0.0518864591774790	3.77174213741851e-10	0.000580831681939922	7.97405604743620e-30	0.0660259803416299	6.29137955059571e-12	0.000500238717434247	2.78157203235329e-05	3.67744070950491e-05	0.706675337480128	0.00104581188589214	7.94060469213274e-06	0.0634211844981511	7.77022364283257e-19	0.0558005175443096	0.0539415131008550	3.43420457005632e-05	1.52362319447776e-05]
    [2.02934951680463e-05	0.0819514426812765	0.00587746530272417	1.12942910915840e-08	7.07682025375137e-07	0.000219700599512478	0.229511658185978	0.209524004136316	2.35546751056699e-19	0.286199014662333	0.00126960750520102	0.0892967894625532	0.00508817870175021	2.91561618504951e-24	0.0895862772514099	6.47343828059709e-06	5.67384242548235e-09	0.000224094066002085	0.00122301803539802	1.25782593868869e-06]])))

(def mu0 [0 0 0]);[-0.2149 -0.0177 0.7630]
(def sig-0 1); 0.01
(def sig-q 0.01)
(def sig-y 0.2)

(defdist initial-dist
  []
  [dist (normal 0 sig-0)]
  (sample* [this] (into [] (mapv #(+ % (sample* dist)) mu0)))
  (observe* [this value]
           (reduce + (map #(observe* dist (- %1 %2))
                          mu0
                          value))))


(defdist obs-dist
  [x] ; distribution parameters
  [dist (normal 0 sig-y)]        ; auxiliary bindings
  (sample* [this] (mapv #(+ % (sample* dist)) (mmul obs-matrix x)))
  (observe* [this value]
           (reduce + (map #(observe* dist (- %1 %2))
                          (mmul obs-matrix x)
                          value))))

(defdist trans-dist
  [beta nu x]
  [dist (normal 0 sig-q)]
  (sample* [this] (mapv #(+ % (sample* dist))
                       (transition-dist-set-clj beta nu (first x) (second x) (nth x 2))))
  (observe* [this value]
           (reduce + (map #(observe* dist (- %1 %2))
                          (transition-dist-set-clj beta nu (first x) (second x) (nth x 2))
                          value))))

(defdist nested-query
  [alg query options n-samples print-folder]
  []
  (sample* [this] nil)
  (observe* [this value] (let [samples (->> (apply doquery alg query value options)
                                           (take n-samples)
                                           doall)
                              target-folder (str "chaos-results/" print-folder "/pmmh" )
                              log-Z (log-marginal samples)
                              now (System/currentTimeMillis)]
                          (.mkdir (java.io.File. (str "chaos-results/" print-folder)))
                          (.mkdir (java.io.File. target-folder))
                          (spit (str target-folder "/logZ-" now ".csv")
                                (write-csv (map (comp vector str) [[log-Z]])))
                          (print :log-Z-pmmh log-Z)
                          log-Z)))

(with-primitive-procedures [matrix initial-dist obs-dist trans-dist nested-query]
  (defquery kalman-chaos-query-fixed
    [observations beta nu]
    (let [trans-sample (fn [x] (sample (trans-dist beta nu x)))]
      (predict :beta beta)
      (predict :nu nu)
      (predict :states
               (matrix
                (reduce (fn [states obs]
                          (let [;; sample next state
                                prev-state (peek states)
                                state (if prev-state
                                        (trans-sample prev-state)
                                        (sample (initial-dist)))]
                            ;; observe next data point (when available)
                            (observe (count states) (obs-dist state) obs)
                            ;; append state to sequence and continue with next obs
                            (conj states state)))
                        ;; start with empty sequence
                        []
                        ;; loop over data
                        observations)))))

  (defquery kalman-chaos-query-sample
    [observations]
    (let [beta (sample (uniform-continuous -3 3))
          nu (sample (uniform-continuous 0 3))
          trans-sample (fn [x] (sample (trans-dist beta nu x)))]
      (predict :beta beta)
      (predict :nu nu)
      (predict :states
               (matrix
                (reduce (fn [states obs]
                          (let [;; sample next state
                                prev-state (peek states)
                                state (if prev-state
                                        (trans-sample prev-state)
                                        (sample (initial-dist)))]
                            ;; observe next data point (when available)
                            (observe (count states) (obs-dist state) obs)
                            ;; append state to sequence and continue with next obs
                            (conj states state)))
                        ;; start with empty sequence
                        []
                        ;; loop over data
                        observations)))))

  (defquery kalman-chaos-pmmh
    "Does PMMH for kalman chaos by using a nested query"
    [observations target-alg target-options n-samples-target print-folder]
    (let [beta (sample (uniform-continuous -3 3))
          nu (sample (uniform-continuous 0 3))]
      (predict :beta beta)
      (predict :nu nu)
      (observe (nested-query target-alg kalman-chaos-query-fixed target-options n-samples-target print-folder)
               [observations beta nu])))

  (defopt kalman-chaos-opt
    [observations]
    [beta nu]
    (let [beta (sample (uniform-continuous -3 3))
          nu (sample (uniform-continuous 0 3))
          trans-sample (fn [x] (sample (trans-dist beta nu x)))]
      (predict :beta beta)
      (predict :nu nu)
      (predict :states
               (matrix
                (reduce (fn [states obs]
                          (let [;; sample next state
                                prev-state (peek states)
                                state (if prev-state
                                        (trans-sample prev-state)
                                        (sample (initial-dist)))]
                            ;; observe next data point (when available)
                            (observe (count states) (obs-dist state) obs)
                            ;; append state to sequence and continue with next obs
                            (conj states state)))
                        ;; start with empty sequence
                        []
                        ;; loop over data
                        observations))))))

(def number-of-particles 500)
(def number-of-samples 75000)
(def T 500)

(defn do-inf-get-best [folder-name n-data inf-alg extract-keys]
  (let [;; Runs SMC inference, extracts highest weight and rerurns

        folder-name (str folder-name "-" n-data)
        _ (.mkdir (java.io.File. (str "chaos-results/" folder-name)))

        observations (into [] (chaos-data n-data T))
        samples-inference
          (->> (doquery inf-alg kalman-chaos-query-sample
                        [observations]
                        :number-of-particles number-of-particles)
               (take number-of-samples)
               doall
             time)
        _ (println :out)
        _ (println :log-marginal (log-marginal samples-inference))
        log-weights-all (map :log-weight samples-inference)
        log-weights-sweeps (take-nth number-of-particles log-weights-all)
        _ (println :log-weights-sweeps log-weights-sweeps)

        weights
            (let [max-log-weights (apply max log-weights-all)
                  log-weights-c (m/sub log-weights-all max-log-weights)
                  weights (m/exp log-weights-c)
                  weights (m/div weights (reduce + weights))]
              weights)
        best-ind (argmax weights)
        ;(println :id-best best-ind)
        best-sample (nth samples-inference best-ind)

        _ (println :log-weight-best (:log-weight best-sample))
        ;(println :best-sample best-sample)
        preds-best (get-predicts best-sample)

        extracts (mapv #(get preds-best %) extract-keys)
        _ (println :extracts-best extracts)


        samples-set
            (->> (doquery :smc kalman-chaos-query-fixed
                          (into [] (concat [observations]
                                           extracts))
                        :number-of-particles number-of-particles)
                 (take (max number-of-particles (int (/ number-of-samples 10))))
                 doall
                 time)

        preds-samples-set (get-predicts (last samples-set))

        _ (println :log-marginal-set-sample (log-marginal samples-set))

        extract-preds-inference (map (fn [x] (mapv #(get (get-predicts x) %) extract-keys)) (take-nth number-of-particles samples-inference))

        output [{:log-marginal-query (log-marginal samples-inference)}
                {:mean-log-marginal-query (mean (mapv #(:log-weight %) samples-inference))}
                {:best-log-weight (:log-weight best-sample)}
                {:extracts-best extracts}
                {:log-marginal-resampled (log-marginal samples-set)}
                {:mean-log-resampled (mean (mapv #(:log-weight %) samples-set))}
                {:log-weights-sweeps log-weights-sweeps}
                {:extract-preds-inference extract-preds-inference}
                {:predicts preds-samples-set}]
        _ (println :extract-preds-inference extract-preds-inference)
        _ (println :mean-log-marginal-query (mean (mapv #(:log-weight %) samples-inference)))]
    (spit (str "chaos-results/" folder-name "/inference-results-" n-data ".json") (json/write-str (str output)))
  nil))


(defn pmmh-chaos [folder-name n-data outer-inf-alg extract-keys]
  (let [;; Runs PMMH inference

        folder-name (str folder-name "-" n-data)

        _ (.mkdir (java.io.File. (str "chaos-results/" folder-name)))

        observations (into [] (chaos-data n-data T))

        samples-pmmh
          (->> (doquery outer-inf-alg
                        kalman-chaos-pmmh
                        [observations :smc
                         [:number-of-particles number-of-particles]
                         number-of-particles folder-name])
               (take (int (/ number-of-samples number-of-particles)))
               doall
             time)

        extract-preds-pmmh (map (fn [x] (mapv #(get (get-predicts x) %) extract-keys))
                                samples-pmmh)

        output [{:extract-preds-pmmh extract-preds-pmmh}]
        _ (println :output (str output))]
    (spit (str "chaos-results/" folder-name "/pmmh-results-" (subs (str outer-inf-alg) 1) "-" n-data ".json")
          (json/write-str (str output)))
  nil))

(defn get-samples-opt [folder-name n-data]
  (let [; Runs BOPP

        folder-name (str folder-name "-" n-data)

        _ (.mkdir (java.io.File. (str "chaos-results/" folder-name)))

        observations (into [] (chaos-data n-data T))

        samples (->> (doopt :smc kalman-chaos-opt
                            [observations]
                            number-of-particles
                            :bo-options {:verbose 1
                                         :bo-debug-folder (str folder-name n-data)})
                     (take (int (/ number-of-samples number-of-particles)))
                     doall)
        points (mapv #(vec [(first %) (second %)]) samples)

        output [{:theta-log-Z-pairs points}]
        ]
  (spit (str "chaos-results/" folder-name "/bopp-results-" n-data ".json") (json/write-str (str output)))
  nil))

(def folder-base "light-test")

(defn -main [folder-name n-run]
  (let [folder-name (str folder-name "-N-" number-of-particles "-T-" T "-" (System/currentTimeMillis) "-" n-run)]
      (get-samples-opt folder-name 1)
      (do-inf-get-best folder-name 1 :smc [:beta :nu])
      (pmmh-chaos folder-name 1 :rmh [:beta :nu])
      (pmmh-chaos folder-name 1 :lmh [:beta :nu])))

;(-main "quick-test" "1")
