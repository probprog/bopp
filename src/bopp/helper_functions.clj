(ns bopp.helper-functions
  "Helper functions"
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.linear :as ml]))

;; matrix library uses vectorz for protocol implementations
(m/set-current-implementation :vectorz)

(defn argmax
  "Index of maximum of a collection"
  [coll]
  (first
   (apply max-key
          second
          (map vector
               (range (count coll))
               coll))))

(defn indexed-max
  "Returns an indexed maximum. Accepts a function f and a collection
  coll. Returns a pair [y-max i-max] in which y-max is the largest
  value (f x-max) and i-max is the index such that (nth coll i-max)
  returns x-max."
  [f coll]
  (loop [best [(/ 1 (- 0.0)) 0]
         i 0
         xs coll]
    (if-let [x (first xs)]
      (let [[y-max i-max] best
            y (f x)]
        (recur (if (> y y-max)
                 [y i]
                 best)
               (inc i)
               (rest xs)))
      best)))

(defn cartesian [colls]
  (if (empty? colls)
    '(())
    (for [x (first colls)
          more (cartesian (rest colls))]
      (cons x more))))

;; Functions ported from Anglican

;; erf
(defn erf
  "error function"
  [x]
  (org.apache.commons.math3.special.Erf/erf x))

;; mean
(defn sum
  "sums array slices along specified dimension"
  ([a dimension]
   (reduce
     m/add
     (m/slices a dimension)))
  ([a]
   (sum a 0)))

(defn mean
  "mean of array slices along specified dimension"
  ([a dimension]
   (m/div (sum a dimension)
          (get (m/shape a) dimension)))
  ([a]
   (mean a 0)))

;; distributions
(defprotocol distribution
  "random distribution"
  (sample* [this]
    "draws a sample from the distribution")
  (observe* [this value]
    "return the probability [density] of the value"))

(def RNG
  "random number generator;
  used by Apache Commons Math distribution objects"
  (org.apache.commons.math3.random.SynchronizedRandomGenerator.
    (org.apache.commons.math3.random.Well19937c.)))

(defn ^:private qualify
  "accepts a symbol, returns the qualified symbol;
  intended to be called from a macro"
  [s]
  (symbol (format "%s/%s" *ns* s)))

(defmacro defdist
  "defines distribution"
  [name & args]
  (let [[docstring parameters & args]
        (if (string? (first args))
          args
          `(~(format "%s distribution" name) ~@args))
        [bindings & methods]
        (if (vector? (first args))
          args
          `[[] ~@args])
        record-name (symbol (format "%s-distribution" name))
        variables (take-nth 2 bindings)]
    `(do
       (declare ~name)
       (defrecord ~record-name [~@parameters ~@variables]
         Object
         (toString [~'this]
           (str (list '~(qualify name) ~@parameters)))
         distribution
         ~@methods)
       (defn ~name ~docstring ~parameters
         (let ~bindings
           (~(symbol (format "->%s" record-name))
                     ~@parameters ~@variables)))
       (defmethod print-method ~record-name
         [~'o ~'m]
         (print-simple (str ~'o) ~'m)))))

(defmacro ^:private from-apache
  "wraps Apache Commons Math distribution"
  [name args type [apache-name & apache-args]]
  (let [dist (gensym "dist")]
    `(defdist ~(symbol name)
       ~(format "%s distribution (imported from apache)" name)
       ~args
       [~dist (~(symbol (format "org.apache.commons.math3.distribution.%sDistribution." apache-name))
                        RNG ~@apache-args)]
       (~'sample* [~'this] (.sample ~dist))
       (~'observe* [~'this ~'value]
         ~(case type
            :discrete `(~'.logProbability ~dist ~'value)
            :continuous `(~'.logDensity ~dist ~'value))))))

(defprotocol multivariate-distribution
  "additional methods for multivariate distributions"
  (transform-sample [this samples]
    "accepts a vector of random values and generates
    a sample from the multivariate distribution"))

(from-apache normal [mean sd] :continuous
  (Normal (double mean) (double sd)))

(defdist mvn
  "multivariate normal"
  [mean cov] [k (m/ecount mean)     ; number of dimensions
              Lcov (:L (ml/cholesky (m/matrix cov)))
              unit-normal (normal 0 1)
              Z (delay (let [|Lcov| (reduce * (m/diagonal Lcov))]
                         (+ (* 0.5 k (Math/log (* 2 Math/PI)))
                            (Math/log |Lcov|))))
              iLcov (delay (m/inverse Lcov))
              transform-sample (fn [samples]
                                 (m/add mean (m/mmul Lcov samples)))]
  (sample* [this] (transform-sample
                   (repeatedly k #(sample* unit-normal))))
  (observe* [this value]
           (let [dx (m/mmul @iLcov (m/sub value mean))]
             (- (* -0.5 (m/dot dx dx)) @Z)))
  multivariate-distribution
  (transform-sample [this samples] (transform-sample samples)))
