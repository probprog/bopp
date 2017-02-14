(ns examples.hmm
  (:require [anglican.core :refer [doquery]]
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
          [state :only [get-predicts get-log-weight]]
          [inference :only [log-marginal rand-roulette]]]))


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

(defopt hmm-simple-opt
  [observations]
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
  (predict
    :states
    (reduce
      (fn [states obs]
        (let [state (sample (trans-dist-mem (peek states)))]
          (observe (normal (nth mus state) sig) obs)
          (conj states state)))
      [(sample init-dist)]
      observations))
  (predict
   :n-states n-states)
  (predict
   :mus mus)
  (predict
   :transition-dist {0 (trans-dist-mem 0)
                     1 (trans-dist-mem 1)
                     2 (trans-dist-mem 2)
                     3 (trans-dist-mem 3)
                     4 (trans-dist-mem 4)})))

; Commented out as currently unused
(defquery hmm-simple-query
  [observations]
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
  (predict
    :states
    (reduce
      (fn [states obs]
        (let [state (sample (trans-dist-mem (peek states)))]
          (observe (normal (nth mus state) sig) obs)
          (conj states state)))
      [(sample init-dist)]
      observations))
  (predict
   :n-states n-states)
  (predict
   :mus mus)
  (predict
   :transition-dist {0 (trans-dist-mem 0)
                     1 (trans-dist-mem 1)
                     2 (trans-dist-mem 2)
                     3 (trans-dist-mem 3)
                     4 (trans-dist-mem 4)})))

(defquery hmm-set-query
  [observations n-states mus]
  (let [init-dist (discrete (apply conj [1] (repeat (dec  n-states) 0)))
        trans-dist-mem (mem (fn [n] (discrete (sample (dirichlet (into [] (repeat n-states 1))))))) ]
  (predict
    :states
    (reduce
      (fn [states obs]
        (let [state (sample (trans-dist-mem (peek states)))]
          (observe (normal (nth mus state) sig) obs)
          (conj states state)))
      [(sample init-dist)]
      observations))
    (predict
   :n-states n-states)
  (predict
   :mus mus)
  (predict
   :transition-dist {0 (trans-dist-mem 0)
                     1 (trans-dist-mem 1)
                     2 (trans-dist-mem 2)
                     3 (trans-dist-mem 3)
                     4 (trans-dist-mem 4)})))

(defdist nested-query
  [alg query options n-samples print-folder]
  []
  (sample* [this] nil)
  (observe* [this value] (let [samples (->> (apply doquery alg query value options)
                                           (take n-samples)
                                           doall)
                              target-folder (str "hmm-results/" print-folder "/pmmh" )
                              log-Z (log-marginal samples)
                              now (System/currentTimeMillis)]
                          (.mkdir (java.io.File. (str "hmm-results/" print-folder)))
                          (.mkdir (java.io.File. target-folder))
                          (spit (str target-folder "/logZ-" now ".csv")
                                (write-csv (map (comp vector str) [[log-Z]])))
                          (print :log-Z-pmmh log-Z)
                          log-Z)))

(with-primitive-procedures [nested-query]
  (defquery hmm-pmmh
    "Does PMMH for hmm by using a nested query"
    [observations target-alg target-options n-samples-target print-folder]
    (let [opt-min (apply min observations)
        opt-max (apply max observations)
        n-states (sample (uniform-discrete 1 6))

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
      (predict :n-states n-states)
      (predict :mus mus)
      (observe (nested-query target-alg hmm-set-query target-options n-samples-target print-folder)
               [observations n-states mus]))))

(def number-of-particles 500)
(def number-of-samples 75000)

(defn pmmh-hmm [folder-name n-data outer-inf-alg extract-keys]
  (let [;; Runs PMMH inference

        folder-name (str folder-name "-" n-data)

        _ (.mkdir (java.io.File. (str "hmm-results/" folder-name)))

        observations (second (hmm-data n-data))

        samples-pmmh
          (->> (doquery outer-inf-alg
                        hmm-pmmh
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
    (spit (str "hmm-results/" folder-name "/pmmh-results-" (subs (str outer-inf-alg) 1) "-" n-data ".json")
          (json/write-str (str output)))
  nil))

(defn do-inf-get-best [folder-name n-data target inf-alg extract-keys]
  (let [folder-name (str folder-name "-" n-data)
        _ (.mkdir (java.io.File. (str "hmm-results/" folder-name)))

        observations (second (hmm-data n-data))

        samples-query
          (->> (doquery inf-alg target
                        [observations]
                        :number-of-particles number-of-particles)
               (take number-of-samples)
               doall
             time)
        _ (println :log-marginal (log-marginal samples-query))
        log-weights (map :log-weight samples-query)
        log-weights-sweeps (take-nth number-of-particles log-weights)
        weights
            (let [max-log-weights (apply max log-weights)
                  log-weights-c (m/sub log-weights max-log-weights)
                  weights (m/exp log-weights-c)
                  weights (m/div weights (reduce + weights))]
              weights)
        best-ind (argmax weights)
        ;(println :id-best best-ind)
        best-sample (nth samples-query best-ind)
        _ (println :log-weight-best (:log-weight best-sample))
        ;(println :best-sample best-sample)
        preds-best (get-predicts best-sample)

        extracts (mapv #(get preds-best %) extract-keys)
        _ (println :extracts-best extracts)

        n-states (:n-states preds-best)
        mus (:mus preds-best)
        samples-set
            (->> (doquery :smc hmm-set-query
                          (into [] (concat [observations]
                                           extracts))
                        :number-of-particles number-of-particles)
                 (take number-of-samples)
                 doall
                 time)

        preds-samples-set (get-predicts (last samples-set))

        _ (println :log-marginal-set-sample (log-marginal samples-set))

        extract-preds-inference (map (fn [x] (mapv #(get (get-predicts x) %) extract-keys)) (take-nth number-of-particles samples-query))

        output [{:log-marginal-query (log-marginal samples-query)}
                {:mean-log-marginal-query (mean (mapv #(:log-weight %) samples-query))}
                {:best-log-weight (:log-weight best-sample)}
                {:extracts-best extracts}
                {:log-marginal-resampled (log-marginal samples-set)}
                {:mean-log-resampled (mean (mapv #(:log-weight %) samples-set))}
                {:log-weights-sweeps log-weights-sweeps}
                {:extract-preds-inference extract-preds-inference}
                {:predicts preds-samples-set}]
        _ (println :output (str output))]
    (spit (str "hmm-results/" folder-name "/inference-results-" n-data ".json") (json/write-str (str output)))
  output))

(defn get-samples-opt [folder-name n-data]
  (let [folder-name (str folder-name "-" n-data)

        _ (.mkdir (java.io.File. (str "hmm-results/" folder-name)))
        samples (->> (doopt :smc hmm-simple-opt
                            [(second (hmm-data n-data))]
                            number-of-particles
                            :bo-options {:verbose 1
                                         :debug-folder (str folder-name n-data)})
                     (take (int (/ number-of-samples number-of-particles)))
                     doall)
        points (mapv #(vec [(first %) (second %)]) samples)

        output [{:theta-log-Z-pairs points}]
        ]
  (spit (str "hmm-results/" folder-name "/bopp-results-" n-data ".json") (json/write-str (str output)))
  nil))

(def target-folder "hmm-var-nstates")

(defn -main [folder-name n-run]
  (let [folder-name (str folder-name "-N-" number-of-particles "-" (System/currentTimeMillis) "-" n-run)]
      (.mkdir (java.io.File. "hmm-results"))
      (get-samples-opt folder-name 5)
      (do-inf-get-best folder-name 5 hmm-simple-query :smc [:n-states :mus])
      (pmmh-hmm folder-name 5 :rmh [:n-states :mus])
      (pmmh-hmm folder-name 5 :lmh [:n-states :mus])))

;(-main "dudtest" "1")
