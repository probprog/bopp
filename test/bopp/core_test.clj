(ns bopp.core-test
  (:require [clojure.test :refer :all]
            [bopp.core :refer :all]
            [anglican.runtime :refer :all]
            [anglican.smc]
            [anglican.importance]
            [anglican.inference :refer [infer]]))

(defopt q [y] [theta]
  (let [a (sample (normal 0 1))
        theta (sample (gamma 1 1))
        b (sample (normal a theta))]
    (observe (normal (+ a b) theta) y)
    [a b]))

(deftest doopt-test
  (let [mmap-states (doopt :smc q [2] 100 :bo-options {:verbose true} :opt-type :mmap)
        ml2-states (doopt :smc q [2] 100 :bo-options {:verbose true} :opt-type :ml2)]
    (is (first mmap-states)
        "testing whether mmap runs.")
    (is (first ml2-states)
        "testing whether ml2 runs.")))

;; ;; This should throw an Exception
;; (defopt q2 [y] [theta]
;;   (let [a (sample (normal 0 1))
;;         theta (sample (gamma 1 1))
;;         b (sample (normal a theta))
;;         theta (sample (discrete 1 1))]
;;     (observe (normal (+ a b) theta) y)
;;     [a b]))

;; ;; This shouldn't throw an Exception
;; (defopt q3 [y] [theta]
;;   (let [a (sample (normal 0 1))
;;         theta (sample (gamma 1 1))
;;         b (sample (normal a theta))
;;         theta (sample (normal 1 1))]
;;     (observe (normal (+ a b) theta) y)
;;     [a b]))

;; ;; This should detect multiple instances of declarations of theta during runtime
;; (defopt q4 [y] [theta]
;;   (let [a (sample (normal 0 1))
;;         theta (sample (gamma 1 1))
;;         b (sample (normal a theta))
;;         theta (sample (normal 1 1))]
;;     (observe (normal (+ a b) theta) y)
;;     [a b]))
;; (def mmap-states (doopt :smc q4 [2] 100 :bo-options {:verbose true} :opt-type :mmap))

;; ;; This should be slow because there is an expensive operation before sampling all optim-vars
;; (defopt q5 [y] [theta psi]
;;   (let [a (sample (normal 0 1))
;;         psi (sample (normal 1 2))
;;         x (doall (repeatedly 1e5 #(sample (mvn [1 0] [[1 0] [0 1]]))))
;;         theta (sample (gamma 1 1))
;;         b (sample (normal a theta))]
;;     (observe (normal (+ a b) theta) y)
;;     [a b]))

;; (time (first (infer :importance (:prior-query q5) [2])))

;; ;; This should be fast
;; (defopt q6 [y] [theta psi]
;;   (let [a (sample (normal 0 1))
;;         psi (sample (normal 1 2))
;;         theta (sample (gamma 1 1))
;;         x (doall (repeatedly 1e5 #(sample (mvn [1 0] [[1 0] [0 1]]))))
;;         b (sample (normal a theta))]
;;     (observe (normal (+ a b) theta) y)
;;     [a b]))

;; (time (first (infer :importance (:prior-query q6) [2])))

;; (run-tests 'bopp.core-test)
