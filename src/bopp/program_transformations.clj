(ns bopp.program-transformations
  "Program transformations for BOPP"
  (:require [anglican.emit :refer [query with-primitive-procedures]]
            [anglican.runtime :refer [defdist]]))

(def optim-args-key "OPTIM_ARGS")
(def acq-f-name "ACQ_F")

(defn- update-values
  "Updates each value of a map m by applying f.

  Adapted from
  http://blog.jayfields.com/2011/08/clojure-apply-function-to-each-value-of.html"
  [m f & args]
  (reduce (fn [r [k v]] (assoc r k (apply f v args))) {} m))

(defn- print-ast [ast level]
  (if (list? ast)
    (do
      (prn (str (apply str (repeat level "-")) (first ast)))
      (loop [arguments (rest ast)]
        (if (not (empty? arguments))
          (do
            (print-ast (first arguments) (inc level))
            (recur (rest arguments))))))
    (prn (str (apply str (repeat level "-")) ast))))

(defn- append-list [lst item]
  (apply list (conj (vec lst) item)))

(defn- retrieve-optim-args-code [optim-args]
  `(~'map (~'fn [~'x]
                (~'retrieve (~'symbol ~optim-args-key)
                            ~'x))
          ~(mapv #(list 'symbol (str %)) optim-args)))

(defn- append-optim-predict [asts optim-args]
  "Appends a predict statement with optimization arguments to the query
  represented by an abstract syntax tree.
  Used in prior-transformation."
  (append-list asts
               (retrieve-optim-args-code optim-args)))

(defn- acq-prologue-code [optim-args acq-f]
  `(do
     (~'observe (~'factor) (~acq-f ~(retrieve-optim-args-code optim-args)))
     ~(retrieve-optim-args-code optim-args)))

(defn- prior-prologue-code [optim-args]
  (retrieve-optim-args-code optim-args))

(defn- acq-prior-helper-transformation
  "Helper transformation used in prior-transformation and acq-transformation.

  Takes in
    expr         a list representing the abstract syntax tree of the program
    optim-args   list of symbols of optimization variables
    mode         either \"acq\" or \"prior\"; indicates which transformation
                 called this function. if mode is \"acq\", a keyword argument
                 acq-f must be provided.

  Outputs an abstract syntax tree in which
    1) all observe statements are removed,
    2) error checking during optim-arg binding is done,
    3) sampled optim-arg optim-arg are stored in a store-retrieve hashmap,
    4) early stopping via return statement is triggered when all optim-args
       are sampled."
  [expr optim-args mode & {:keys [acq-f] :or {acq-f nil}}]
  (cond
   (vector? expr) (mapv #(acq-prior-helper-transformation % optim-args mode :acq-f acq-f) expr)
   (map? expr) (update-values expr #(acq-prior-helper-transformation % optim-args mode :acq-f acq-f))
   (set? expr) (set (map #(acq-prior-helper-transformation % optim-args mode :acq-f acq-f) expr))
   (seq? expr)
   (let [[kwd & args] expr]
     (case kwd
       observe nil
       let (loop [bindings (first args)
                  transformed-bindings []]
             (if (empty? bindings)
               (apply list kwd (vec transformed-bindings) (map #(acq-prior-helper-transformation % optim-args mode :acq-f acq-f) (rest args)))
               (let [[name value] (take 2 bindings)]
                 (recur (drop 2 bindings)
                        (concat transformed-bindings [name (if (some #(= name %) optim-args)
                                                             `(~'let [~'value ~(acq-prior-helper-transformation value optim-args mode :acq-f acq-f)]

                                                                     ;; Error checking for detecting multiple instances of declaration of optim-vars during runtime.
                                                                     (~'if (~'retrieve (~'symbol ~optim-args-key)
                                                                                       (~'symbol ~(str name)))
                                                                           (~'throw-exception ~(str "BOPP ERROR: Multiple instances of declaration of optimization variable " (str name) " detected!")))

                                                                     ;; Store the sampled value in store-retrieve hashmap
                                                                     (~'store (~'symbol ~optim-args-key)
                                                                              (~'symbol ~(str name))
                                                                              ~'value)

                                                                     (~'if (~'= (~'set (~'keys (~'retrieve (~'symbol ~optim-args-key))))
                                                                                (~'set ~(mapv #(list 'symbol (str %)) optim-args)))
                                                                           (~'return ~(cond (= mode "acq")
                                                                                            (acq-prologue-code optim-args acq-f)

                                                                                            (= mode "prior")
                                                                                            (prior-prologue-code optim-args))))

                                                                     ~'value)
                                                             (acq-prior-helper-transformation value optim-args mode :acq-f acq-f))])))))
       loop (loop [bindings (first args)
                   transformed-bindings []]
              (if (empty? bindings)
                (apply list kwd (vec transformed-bindings) (map #(acq-prior-helper-transformation % optim-args mode :acq-f acq-f) (rest args)))
                (let [[name value & bindings] bindings]
                  (recur bindings
                         (concat transformed-bindings
                                 [name (acq-prior-helper-transformation value optim-args mode :acq-f acq-f)])))))

       ;; else
       (map #(acq-prior-helper-transformation % optim-args mode :acq-f acq-f) expr)))
   :else expr))

(defn- prior-transformation
  [asts optim-args]
  "Takes in a list of lists representing the abstract syntax trees of the
  program.
  Outputs an abstract syntax tree in which all observe statements are removed
  and a predict statement with optimization variables is appended."
  (append-optim-predict
   (acq-prior-helper-transformation (apply list (conj asts 'do)) optim-args "prior")
   optim-args))

(defn throw-exception [msg] (throw (Exception. msg)))

(defmacro prior-query
  "Returns CPS'd query with the prior transformation applied."
  [& args]
  (let [value (first args)
        optim-args (second args)
        source (drop 2 args)
        transformed-source (prior-transformation source optim-args)]
    `(with-primitive-procedures [throw-exception]
       (query ~value
              ~@(list transformed-source)))))

(defn- append-acq-factor [asts optim-args acq-f]
  "Appends an observe statement which factors an acquisition function output
  on optim-args input.
  Used in acq-transformation."
  (append-list asts `(~'observe (~'factor) (~acq-f ~(retrieve-optim-args-code optim-args)))))

(defn- acq-transformation
  [asts optim-args acq-f]
  "Takes in a list of lists representing the abstract syntax trees of the
  program.
  Outputs an abstract syntax tree in which all observe statements are removed
  and the following is appended:
  1) observe statement which factors an acquisition function output on
  optim-args input, and
  2) predict statement with optimization variables."
  (-> (acq-prior-helper-transformation (apply list 'do asts) optim-args "acq" :acq-f acq-f)
      (append-acq-factor optim-args acq-f)
      (append-optim-predict optim-args)))

(defdist factor
  "Factor the value to the log-likelihood terms.
  e.g. (observe (factor) 1) factors the term 1 in the log likelihood as if
  there was an observed value with the likelihood e."
  [] []
  (sample* [this] nil)
  (observe* [this value] value))

(defmacro acq-query
  "Returns CPS'd query with the acquisition transformation applied."
  [& args]
  (let [value (first args)
        acq-f (symbol acq-f-name)
        optim-args (second args)
        source (drop 2 args)
        transformed-source (acq-transformation source optim-args acq-f)]
    `(with-primitive-procedures [factor throw-exception]
       (query ~(conj value acq-f)
              ~@(list transformed-source)))))

(defn- marg-transformation
  "Marginal transformation used in mmap-query and ml2-query.

  Takes in
    expr         a list representing the abstract syntax tree of the program
    optim-args   list of symbols of optimization variables
    mode         either \"mmap\" or \"ml2\"; indicates which transformation
                 should be performed

  Outputs an abstract syntax tree in which

  1) mode = \"mmap\":
     sample statements being assigned to any of the optimization variables are
     replaced by a value-returning observe.
  2) mode = \"ml2\":
     sample statements being assigned to any of the optimization variables are
     replaced by a deterministic assignment."
  [expr optim-args mode]
  (cond
    (vector? expr) (mapv #(marg-transformation % optim-args mode) expr)
    (map? expr) (update-values expr #(marg-transformation % optim-args mode))
    (set? expr) (set (map #(marg-transformation % optim-args mode) expr))
    (seq? expr)
      (let [[kwd & args] expr]
        (case kwd
          let (loop [bindings (first args)
                     transformed-bindings []]
                (if (empty? bindings)
                  (apply list kwd (vec transformed-bindings) (map #(marg-transformation % optim-args mode) (rest args)))
                  (let [[name value & bindings] bindings]
                    (recur bindings
                           (concat transformed-bindings
                                   [name (if (some #(= name %) optim-args)
                                           (if (= (first value) 'sample)
                                             `(~'do
                                               ;; Error checking for detecting instances of declaration of optim-vars during runtime.
                                               (~'if (~'retrieve (~'symbol ~optim-args-key)
                                                                 (~'symbol ~(str name)))
                                                     (~'throw-exception ~(str "BOPP ERROR: Multiple instances of declaration of optimization variable " (str name) " detected!")))

                                               ;; Store the sampled value in store-retrieve hashmap
                                               (~'store (~'symbol ~optim-args-key)
                                                        (~'symbol ~(str name))
                                                        ~(symbol (str (str name) "-hat")))

                                               ~(cond (= mode "mmap")
                                                      `(~'observe ~(marg-transformation (second value) optim-args mode) ~(symbol (str (str name) "-hat")))

                                                      (= mode "ml2")
                                                      nil)

                                               ~(symbol (str (str name) "-hat")))
                                             (throw-exception "BOPP ERROR: Prior of optimization variables must be defined directly by a sample statement."))
                                           (marg-transformation value optim-args mode))])))))
          loop (loop [bindings (first args)
                      transformed-bindings []]
                 (if (empty? bindings)
                   (apply list kwd (vec transformed-bindings) (map #(marg-transformation % optim-args mode) (rest args)))
                   (let [[name value & bindings] bindings]
                     (recur bindings
                            (concat transformed-bindings
                                    [name (marg-transformation value optim-args mode)])))))

          ;; else
          (map #(marg-transformation % optim-args mode) expr)))
      :else expr))

(defmacro mmap-query
  "Returns CPS'd query with the MMAP transformation applied."
  [& args]
  (let [value (first args)
        optim-args (second args)
        source (drop 2 args)
        transformed-arguments (vec (concat value
                                           (mapv #(symbol (str (str %) "-hat")) optim-args)))
        transformed-source (marg-transformation (apply list (conj source 'do)) optim-args "mmap")]
    `(with-primitive-procedures [throw-exception]
       (query ~transformed-arguments
              ~@(list transformed-source)))))

(defmacro ml2-query
  "Returns CPS'd query with the maximum likelihood (type 2) transformation applied."
  [& args]
  (let [value (first args)
        optim-args (second args)
        source (drop 2 args)
        transformed-arguments (vec (concat value
                                           (mapv #(symbol (str (str %) "-hat")) optim-args)))
        transformed-source (marg-transformation (apply list (conj source 'do)) optim-args "ml2")]
    `(with-primitive-procedures [throw-exception]
       (query ~transformed-arguments
              ~@(list transformed-source)))))
