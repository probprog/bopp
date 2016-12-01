(ns bopp.program-transformations-test
  (:require [bopp.program-transformations :refer :all]
            [anglican.runtime :refer [normal flip]]))

;;; Prior transformation

;; (def a (prior-query [y] [theta]
;;                     (let [a (sample (normal 0 10))
;;                           theta (sample (normal a 20))
;;                           b (sample (normal a (* theta theta)))]
;;                       (observe (normal (* a theta) (* b b)) y)
;;                       [a b])))

;; (:source (meta a)) ;; Returns:
;; ;; (query [y]
;; ;;        (do
;; ;;          (let [a (sample (normal 0 10))
;; ;;                theta (let [value (sample (normal a 20))]
;; ;;                        (if (retrieve (symbol "OPTIM_ARGS")
;; ;;                                      (symbol "theta"))
;; ;;                          (throw-exception "WARNING: Multiple instances of declaration of optimization variable theta detected!"))

;; ;;                        (store (symbol "OPTIM_ARGS")
;; ;;                               (symbol "theta")
;; ;;                               value)

;; ;;                        (if (= (set (keys (retrieve (symbol "OPTIM_ARGS"))))
;; ;;                               (set [(symbol "theta")]))
;; ;;                          (do
;; ;;                            (println "exiting the program early")
;; ;;                            (return (map (fn [x] (retrieve (symbol "OPTIM_ARGS") x))
;; ;;                                         [(symbol "theta")]))))

;; ;;                        value)
;; ;;                b (sample (normal a (* theta theta)))]
;; ;;            nil
;; ;;            [a b])
;; ;;          (map (fn [x] (retrieve (symbol "OPTIM_ARGS") x)) [(symbol "theta")])))

;;; Acq transformation

;; (def c (acq-query [y] [theta]
;;         (let [a (sample (normal 0 10))
;;               theta (sample (normal a 20))
;;               b (sample (normal a (* theta theta)))]
;;           (observe (normal (* a theta) (* b b)) y)
;;           [a b])))

;; (:source (meta c)) ; Returns:
;; ;; (query [y ACQ_F]
;; ;;        (do (let [a (sample (normal 0 10))
;; ;;                  theta (let [value (sample (normal a 20))]
;; ;;                          (if (retrieve (symbol "OPTIM_ARGS")
;; ;;                                        (symbol "theta"))
;; ;;                            (throw-exception "WARNING: Multiple instances of declaration of optimization variable theta detected!"))

;; ;;                          (store (symbol "OPTIM_ARGS")
;; ;;                                 (symbol "theta")
;; ;;                                 value)

;; ;;                          (if (= (set (keys (retrieve (symbol "OPTIM_ARGS"))))
;; ;;                                 (set [(symbol "theta")]))
;; ;;                            (do (println "exiting the program early")
;; ;;                              (return
;; ;;                               (do
;; ;;                                 (observe (factor) (ACQ_F (map (fn [x] (retrieve (symbol "OPTIM_ARGS") x)) [(symbol "theta")])))
;; ;;                                 (map (fn [x] (retrieve (symbol "OPTIM_ARGS") x))
;; ;;                                      [(symbol "theta")])))))

;; ;;                          value)
;; ;;                  b (sample (normal a (* theta theta)))]
;; ;;              nil
;; ;;              [a b])
;; ;;          (observe (factor) (ACQ_F (map (fn [x] (retrieve (symbol "OPTIM_ARGS") x)) [(symbol "theta")])))
;; ;;          (map (fn [x] (retrieve (symbol "OPTIM_ARGS") x)) [(symbol "theta")])))

;;; MMAP transformation

;; (def d (mmap-query [y] [theta]
;;                    (let [a (sample (normal 0 10))
;;                          c (sample (normal a 20))
;;                          b (sample (normal a (* c c)))
;;                          d (let [theta (sample (normal a (* b c)))]
;;                              theta)]
;;                      (observe (normal (* a d) (* b b)) y)
;;                      [a b])))

;; (:source (meta d))

;;; ML2 transformation

;; (def e (ml2-query [y] [theta]
;;                    (if (sample (flip 0.3))
;;                      (let [a (sample (normal 0 1))
;;                            theta (sample (normal a 1))]
;;                        (observe (normal a (* theta theta)) 1)
;;                        a)
;;                      (let [a (sample (normal 0 12))
;;                            theta (sample (normal a 21))]
;;                        (observe (normal a (* theta theta)) 31)
;;                        theta))))

;; (meta e)


