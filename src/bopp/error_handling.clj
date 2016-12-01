(ns bopp.error-handling
  "Error handling for BOPP")

(defn- dist-type [dist-symbol]
  "Given a distribution object symbol, returns either 'continuous or
  'discrete."
  (cond (some #(= dist-symbol %)
              ['bernoulli 'binomial 'categorical 'discrete 'flip 'poisson
               'uniform-discrete])
        'discrete

        (some #(= dist-symbol %)
              ['beta 'gamma 'dirichlet 'exponential 'chi-squared 'normal
               'uniform-continuous 'mvn 'multivariate-t 'wishart])
        'continuous

        :else
        (throw (Exception. (str "Distribution " dist-symbol " not supported.")))))

(defn- detect-types [ast optim-args types]
  "Given an abstract syntax tree of an optimization query, collect types of
  distribution objects being assigned to optimization variables in optim-args.

  Returns a map with optim-arg as key and {list of types of distributions being
  assigned to this optim-arg} as value. E.g.

  {'theta ['continuous 'continuous 'discrete]
   'psi ['continuous]}.

  Used as a helper function for changing-types?."
  (if (list? ast)
    (let [root (first ast)]
      (if (= root 'let)
        (loop [bindings (second ast)
               ;; Initialize new-types to contain types from statemets after
               ;; the binding vector
               new-types (reduce (partial merge-with concat) {}
                                 (map #(detect-types % optim-args {}) (drop 2 ast)))]
          (if (empty? bindings)
            (merge-with concat types new-types)
            (let [[name value] (take 2 bindings)]
              (recur (drop 2 bindings)
                     (if (some #(= name %) optim-args)
                       (update new-types name #(conj % (dist-type (first (second value)))))
                       (merge-with concat new-types (detect-types value optim-args {})))))))
        (reduce (partial merge-with concat)
                types
                (map #(detect-types % optim-args {}) (rest ast)))))
    types))

(defn changing-types?
  [& args]
  "Given an optimization query, detect whether there are assignments to the
  same optimization variables that change type from continuous to discrete or
  vice-versa.
  Returns true or false."
  (let [value (first args)
        optim-args (second args)
        source (drop 2 args)
        types (detect-types (apply list (conj source 'do)) optim-args {})]
    (not (every? (partial apply =)
                 (vals types)))))
