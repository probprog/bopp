(ns bopp.core
  "User interface for BOPP"
  (require [deodorant.core :as bo]
           [bopp.program-transformations :as pt]
           [bopp.helper-functions :refer [argmax]]
           [bopp.error-handling :refer [changing-types?]]
           [anglican.inference :refer [log-marginal infer]]
           [anglican.state :refer [get-log-weight get-result]]
           [anglican.trap :refer [primitive-procedure-cps]]
           [anglican.ais]
           [anglican.importance]))

(defn get-result-vector [sample]
  (into [] (get-result sample)))

(defmacro opt [& args]
  (if (apply changing-types? args)
    (throw (Exception. "Query contains distribution assignments to optimization variables that have different type (continuous, discrete)."))
    `(with-meta
       {:prior-query (pt/prior-query ~@args)
        :acq-query (pt/acq-query ~@args)
        :mmap-query (pt/mmap-query ~@args)
        :ml2-query (pt/ml2-query ~@args)}
       {:source '(~'opt ~@args)})))

(defmacro defopt
  "Binds variable to a BOPP query. Synatx is the
   same as Anglican's defquery, except that an
   additional input is specified to dictate the
   variables to be optimized.

  Example:

  (defopt opt-query
     <docstring>
     [y] [theta]
     (let [a (sample (normal 0 10))
           theta (sample (normal a 20))
           b (sample (normal a (* theta theta)))]
       (observe (normal (* a theta) (* b b)) y)
       [a b]))


  For more details please see the README and the
  paper.  https://github.com/twgr/bopp/blob/master/README.md
  http://papers.nips.cc/paper/6421-bayesian-optimization-for-probabilistic-programs"
  [name & args]
  (let [[docstring source]
        (if (string? (first args))
          [(first args) (rest args)]
          [(format "BOPP program '%s'" name) args])]
    `(def ~(with-meta name {:doc docstring})
       (opt ~@source))))

(defn bo-target
  "Calculates the value of a the BO target given a marginal query and a
   point to evaluate.  Target will typically be the log-marginal for
   MMAP problems and (- (log-marginal samples)) for risk minimization.
   However, also allows overide capabilities from different transformations
   of the target (e.g. WASABI), to using completely different opt-sample-summarizers for
   the sample populations (e.g. variance on the weights).  Also provides
   flexibility in the additional return values taken (e.g. get-result-vector
   or get-predicts).

  Inputs: marg-q = transformed marginal query
          query-args = fixed arguments to the defopt (i.e. its first argument)
          algorithm = inner inference algorithm to use
          num-samples = number of samples to take in inner inference algorithm
          num-particles = number of particles to use in inner inferenec algorithm
          output-extractor = extracts additional outputs from the inner query,
                 e.g get-result-vector or get-predicts
          opt-bo-target-transformation = transformation function applied to the output of
                 the opt-sample-summarizer calculation.
          opt-sample-summarizer = function used to calculate the desired untransformed target
                 from the samples.
          f-theta-inf = Value of f-theta set to if evaluates as +/- infinity (sign is maintained)
  Outputs: f-theta = Evaluated value as used for the bayes opt.  Note that the bayes
               opt always carries out maximization.  Minimization is thus done
               by maximiziation a function where f-theta = -target.
           inner-outputs = additional outputs from the inner query."
  [theta
   marg-q
   query-args
   algorithm
   num-samples
   num-particles

   opt-bo-target-transformation
   opt-sample-summarizer

   f-theta-inf
   output-extractor]
  (let [samples (take num-samples
                      (infer algorithm
                             marg-q
                             (concat query-args theta)
                             :number-of-particles
                             num-particles))

        f-theta (opt-bo-target-transformation (opt-sample-summarizer samples)) ;; Could be +-Infinity
        f-theta (if (= f-theta (/ -1. 0))
                  (do (println "BOPP WARNING: Target evaluation has underflown!")
                    (- f-theta-inf))
                  (if (= f-theta (/ 1. 0))
                    (do (println "BOPP WARNING: Target evaluation has overflown!")
                      f-theta-inf)
                    f-theta))
        outputs (map output-extractor samples)]
    [f-theta outputs]))

(defn- distributed-aq-optimizer
  "Performs parallelized calls of infer on an optimization query and takes
   the best sample from any of the runs."
  [acq-fn opt-query-args acq-q acq-opt-num-starts ais-num-steps ais-start-exponent ais-end-exponent]
  (let [log-acq-fn (fn [theta & args]
                     (Math/log (apply acq-fn theta args)))
        cpsd-acq-fn (fn [cont state & args]
                             (fn []
                               (cont (apply log-acq-fn args) state)))
        outputs (pmap (fn [unused] (first (infer :ais
                                                 acq-q
                                                 (conj opt-query-args cpsd-acq-fn)
                                                 :number-of-steps ais-num-steps
                                                 :start-exponent ais-start-exponent
                                                 :end-exponent ais-end-exponent)))
                      (range acq-opt-num-starts))
        log-weights (map get-log-weight outputs)
        i-best (argmax log-weights)
        best-output (get-result-vector (nth outputs i-best))]
    best-output))

(defn doopt
  "Calls BOPP on a target query as defined by defopt, returning a lazy
   infinite sequence of increasingly optimal estimates.

  For more details please see the README and the
  paper.  https://github.com/twgr/bopp/blob/master/README.md
  http://papers.nips.cc/paper/6421-bayesian-optimization-for-probabilistic-programs

  Input:
    algorithm ... Inference algorithm
      [:smc, :pcascade, :importance]
    opt-query ... Query to be optimized, defined by defopt.
    opt-query-args ... Arguments of opt-query.
    num-samples ... Number of samples for constructing the estimator to be
      optimized

  Optional keyword arguments:
    Speed option:
      speed-option ... Set the default options to a set of pre-tuned values to
        be able to switch between the different modes of running depending on
        speed. Can be overriden by user post-hoc.
        [:careful, :fast] (Default: :fast)

    Bayesian optimization (BO) options:
      bo-options ... look at documentation of Deodorant, this should be a map
        with those options, i.e.:

          {:initial-points <initial-points>
           :num-scaling-thetas <num-scaling-thetas>
           :num-initial-points <num-initial-points>
           :cov-fn-form <cov-fn-form>
           ... etc.}

        You only need to provide options whose defaults that you want to
        override.

        (Default:
          if :speed-option is :careful
            {:hmc-num-leapfrog-steps 5
             :hmc-num-mcmc-steps 50
             :hmc-num-opt-steps 15
             :hmc-num-chains 8
             :hmc-max-gps 50}

          if :speed-option is :fast
            {:hmc-num-leapfrog-steps 2
             :hmc-num-mcmc-steps 20
             :hmc-num-opt-steps 10
             :hmc-num-chains 4
             :hmc-max-gps 20})

    Optimization options:
      opt-type ... Optimization type. Use :custom to customize the
        optimization; must provide opt-program-transformation,
        opt-sample-summarizer, opt-bo-target-transformation in this case.
        [:mmap, :ml2, :risk-minimization :custom] (Default: :mmap)
      opt-program-transformation ... Custom program transformation for
        optimization.
        [nil, :mmap, :ml2] (Default: nil)
      opt-sample-summarizer ... Custom sample summarizer for optimizer. Takes
        in a sequence of Anglican samples, outputs target scalar for BO.
        (Default: nil; Suggested: anglican.inference/log-marginal)
      opt-bo-target-transformation ... Custom BO target transformation for
        optimization. Takes in BO target, outputs alternative target.
        (Default: nil; Suggested: identity, - for minimization)

    Inference options:
      num-particles ... Number of particles during inference.
        (Default: num-samples)

    Acquisition optimization/AIS options:
      acq-opt-num-starts ... Number of parallel optimizations to prevent
        multimodality. (Default: 4 if :speed-option is :fast 16 if :careful)
      ais-num-steps ... Number of MCMC steps performed during acquisition
        function optimization.
        (Default: (max 100 (int (* 2 (/ num-samples acq-opt-num-starts)))))
      ais-start-exponent ... Start exponent for annealed importance sampling
        (Default: 0.001)
      ais-end-exponent ... End exponent for annealed importance sampling
        (Default: 10)

    Other options
      output-extractor ... Takes in an Anglican sample, returns its output.
        (Default: #(into [] (anglican.state/get-result sample));
         Alternatives: anglican.state/get-result, anglican.state/get-predicts)
      f-theta-inf ... Replace -Infinity during BO by this. (Default: -1e5)

  Output:
    Lazy list of increasingly optimal (in the MMAP sense) triples
    (opt-vars, log-marginal, program-results):
      opt-vars ... optimization variables
      log-marginal ... the log marginal likelihood p(observes, opt-vars)
      program-results ... program results associated with opt-vars

  See paper for details."
  [algorithm opt-query opt-query-args num-samples &
   {:keys [;; Speed option
           speed-option

           ;; BO options
           bo-options

           ;; Optimization options
           opt-type opt-program-transformation opt-sample-summarizer opt-bo-target-transformation

           ;; Inference options
           num-particles ;; FIXME: might want to include all inference options

           ;; Acquisition optimization/AIS options
           acq-opt-num-starts ais-num-steps ais-start-exponent ais-end-exponent

           ;; Other options
           output-extractor f-theta-inf]
    :or {speed-option :fast

         ;; Default optimization options
         opt-type :mmap

         ;; Default acquisition optimization/AIS options
         ais-start-exponent 0.001
         ais-end-exponent 10

         ;; Default other options
         output-extractor get-result-vector
         f-theta-inf -1e5 ; FIXME: Don't replace -Infinity with f-theta-inf; resample a random point instead.
         }}]
  (let [;; Set default parameters based on speed-option
        [acq-opt-num-starts-default
         bo-options-default]
        (case speed-option
          :careful [16 {:hmc-num-leapfrog-steps 5
                       :hmc-num-mcmc-steps 50
                       :hmc-num-opt-steps 15
                       :hmc-num-chains 8
                       :hmc-max-gps 50}]
          :fast [4 {:hmc-num-leapfrog-steps 2
                    :hmc-num-mcmc-steps 20
                    :hmc-num-opt-steps 10
                    :hmc-num-chains 4
                    :hmc-max-gps 20}]
          (throw (Exception. (str "speed-option must be one of [:careful, :fast]. Supplied speed-option: " speed-option))))

        ;; Override these parameters by ones provided by the user
        [acq-opt-num-starts bo-options]
        [(or acq-opt-num-starts acq-opt-num-starts-default)
         (merge bo-options-default bo-options)]

        ;; Optimization options error checking
        _ (if (and (not (= opt-type :custom))
                   (not (some nil? [opt-program-transformation opt-sample-summarizer opt-bo-target-transformation])))
            (throw (Exception. ("If opt-type is not :custom, opt-program-transformation, opt-sample-summarizer- opt-bo-target-transformation must not be supplied."))))

        ;; Override Optimization options
        [opt-program-transformation opt-sample-summarizer opt-bo-target-transformation]
        (case opt-type
          :mmap [:mmap log-marginal identity]
          :ml2 [:ml2 log-marginal identity]
          :risk-minimization [:ml2 log-marginal -]
          :custom [opt-program-transformation opt-sample-summarizer opt-bo-target-transformation]
          (throw (Exception. (str "Optimization type " opt-type " not supported. Supported: :mmap, :ml2, :risk-minimization, :custom."))))

        ;; Optimization options error checking
        _ (if (some nil? [opt-program-transformation opt-sample-summarizer opt-bo-target-transformation])
            (throw (Exception. ("If opt-type is :custom, opt-program-transformation, opt-sample-summarizer- opt-bo-target-transformation must be supplied."))))

        ;; Compiled transformed queries
        prior-q (:prior-query opt-query)
        acq-q (:acq-query opt-query)
        marg-q (case opt-program-transformation
                 :mmap (:mmap-query opt-query)
                 :ml2 (:ml2-query opt-query)
                 (throw (Exception. (str "opt-program-transformation " opt-program-transformation " not supported. Supported: :mmap, :ml2."))))

        ;; Inference option
        num-particles (or num-particles num-samples)

        ;; Setup the target function for BO
        f (fn [theta]
            (bo-target
             theta marg-q opt-query-args algorithm num-samples num-particles opt-bo-target-transformation opt-sample-summarizer f-theta-inf output-extractor))

        ais-num-steps (or ais-num-steps
                          (max 100 (int (* 2 (/ num-samples acq-opt-num-starts)))))

        ;; Setup the acqusition function optimizer
        aq-optimizer (fn [acq-fn]
                       (distributed-aq-optimizer acq-fn opt-query-args acq-q acq-opt-num-starts ais-num-steps ais-start-exponent ais-end-exponent))

        ;; Setup the sampler for cheaply sampling valid theta's
        theta-sampler (fn [n-samples]
                        (map get-result-vector
                             (take n-samples
                                   (infer :importance prior-q opt-query-args))))

        ;; Setup BO options to pass in an apply
        bo-options (flatten (into [] (filter second bo-options)))]
    (apply bo/deodorant f aq-optimizer theta-sampler bo-options)))
