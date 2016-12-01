(defproject bopp "0.1.0-SNAPSHOT"
  :description "Bayesian Optimization for Probabilistic Programs."
  :url "http://www.robots.ox.ac.uk/~twgr/publications/"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :plugins [[lein-codox "0.10.1"]
            [lein-gorilla "0.3.6"]]
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojars.tuananhle/anglican "1.1.0-BOPP"]
                 [deodorant "0.1.0-SNAPSHOT"]
                 [clatrix "0.5.0"]
                 [org.apache.commons/commons-math3 "3.6.1"]
                 [com.taoensso/tufte "1.0.0-RC2"]
                 [org.clojure/data.csv "0.1.3"]
                 [clojure-csv/clojure-csv "2.0.1"]])
