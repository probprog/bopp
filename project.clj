(defproject bopp "0.1.2"
  :description "Bayesian Optimization for Probabilistic Programs."
  :url "http://github.com/probprog/bopp"
  :license {:name "GNU General Public License Version 3"
            :url "http://www.gnu.org/licenses/gpl.html"}
  :plugins [[lein-codox "0.10.2"]
            [lein-gorilla "0.3.6"]]
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojars.tuananhle/anglican "1.1.1-BOPP"]
                 [deodorant "0.1.1"]
                 [clatrix "0.5.0"]
                 [org.apache.commons/commons-math3 "3.6.1"]
                 [com.taoensso/tufte "1.0.0-RC2"]
                 [org.clojure/data.csv "0.1.3"]
                 [clojure-csv/clojure-csv "2.0.1"]])
