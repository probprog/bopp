(ns bopp.error-handling-test
  (:require [clojure.test :refer :all])
  (:require [bopp.error-handling :refer :all]))

(deftest changing-types?-test
  (is (= true
         (apply changing-types? '([y] [theta]
                            (let [a (sample (normal 0 10))
                                  theta (sample (normal a 20))
                                  b (sample (normal a (* theta theta)))
                                  theta (sample (discrete [a 20]))]
                              (observe (normal (* a theta) (* b b)) y)
                              [a b]))))
      (str "normal and discrete are different types."))
  (is (= false
         (apply changing-types? '([y] [theta]
                            (let [a (sample (normal 0 10))
                                  theta (sample (normal a 20))
                                  b (sample (normal a (* theta theta)))
                                  theta (sample (wishart nil nil))]
                              (observe (normal (* a theta) (* b b)) y)
                              [a b]))))
      (str "normal and wishart are same type."))
  (is (= false
         (apply changing-types? '([y] [theta]
                            (let [a (sample (discrete [0 10]))
                                  theta (sample (discrete a 20))
                                  b (sample (normal a (* theta theta)))
                                  theta (sample (categorical nil nil))]
                              (observe (normal (* a theta) (* b b)) y)
                              [a b]))))
      (str "discrete and categorical is different types."))
  (is (= true
         (apply changing-types? '([y] [theta]
                            (let [a (sample (discrete [0 10]))
                                  b (sample (normal a (* theta theta)))]
                              (if condition
                                (let [theta (sample (discrete a 20))]
                                  (+ 1 theta))
                                (let [theta (sample (normal nil nil))]
                                  (+ 2 theta)))
                              (observe (normal (* a theta) (* b b)) y)
                              [a b]))))
      (str "discrete and normal in different branches are different types.")))

;; (run-tests 'bopp.error-handling-test)
