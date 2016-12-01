(ns examples.simple-bimodal
  (require [bopp.core :refer :all]
           [clojure.data.json :as json]
           [anglican.runtime :refer :all]
           [taoensso.tufte :as tufte :refer (defnp p profiled profile)]))


(defopt simple-bimodal [y] [theta]
  (let [theta (sample (normal 0 0.5))]
    (observe
     (normal (sqrt (* theta theta)) 0.5) y)))

(defopt simple-bimodal-noisy [y] [theta]
  (let [theta (sample (normal 0 5))
        sig-n (sample (normal 0 0.35))
        t (+ theta sig-n)]
    (observe
     (normal (sqrt (* t t)) 0.5) y)))

(defn -main [folder-name & opts]
  (let [[num-steps num-samples num-init plot-aq model] (map read-string (take 5 opts))
        num-steps (or num-steps 1)
        num-samples (or num-samples 1000)
        num-init (or num-init 2)
        plot-aq (or plot-aq false)
        model (or (eval model) simple-bimodal)
        _ (println :model-type (type model))]
    (->> (doopt
          :importance model [5]
          num-samples
          :bo-num-initial-thetas num-init
          :bo-debug-folder folder-name
          :bo-plot-aq plot-aq
          :bo-verbose true)
         (take num-steps)
         (json/write-str)
         (spit (str folder-name ".json")))))

(tufte/add-basic-println-handler! {})

(tufte/profile
 {}
   (-main "simple-bimodal-dump" "10" "2" "5" "false" "simple-bimodal"))
