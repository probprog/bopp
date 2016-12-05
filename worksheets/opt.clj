;; gorilla-repl.fileformat = 1

;; **
;;; # Optimization Benchmarks
;; **

;; @@
(ns worksheets.opt
  (:require [gorilla-plot.core :as plot]
            [clojure.core.matrix :as mat]
            [bopp.core :refer :all])
  (:use [anglican runtime emit]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; In this worksheet we consider simply using BOPP as an optimizer on some classic Bayesian optimization benchmark problems.  Results for some popular packages on the same problems are presented in the paper.
;; **

;; @@
(defdist factor
  []
  (sample* [this] 0)
  (observe* [this value] value))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>#multifn[print-method 0x546f29b7]</span>","value":"#multifn[print-method 0x546f29b7]"}
;; <=

;; **
;;; ## Branin
;; **

;; **
;;; @@f(x) = \left(x\_2-\frac{5.1}{4\pi^2}x\_1^2+\frac{5}{\pi}x\_1-6\right)^2+10(1-\frac{1}{8\pi}\cos(x\_1))+10@@
;;; 
;;; Bounds @@-5 \le x\_1 \le 10@@ and @@0 \le x\_2 \le 15@@.
;;; 
;;; There are three global minimum, @@f(x^\*) = 0.397887@@ at @@x^\* = (-\pi,12,1275), (\pi,2.275) ~\rm{and}~ (9.42478,2.475)@@.
;; **

;; @@
(with-primitive-procedures [factor]
  (defopt branin-opt [] [x1 x2]
   (let [x1 (sample (uniform-continuous -5 10))
         x2 (sample (uniform-continuous 0 15))
         t1 (- (+ (- x2 (* (pow x1 2) (/ 5.1 (* 4 (pow Math/PI 2)))))
                  (* x1 (/ 5 Math/PI)))
               6)
         t2 (* 10
               (- 1 (/ 1 (* 8 Math/PI)))
               (cos x1))
         branin (- 0 (+ (pow t1 2) t2 10))
         prior-correction (+ (log 15) (log 15))]
     (observe (factor) (+ branin prior-correction)))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.opt/branin-opt</span>","value":"#'worksheets.opt/branin-opt"}
;; <=

;; @@
(def branin-samples (->> (doopt :importance branin-opt [] 1 :bo-verbose true)
                		 (take 50) ;; Number of optimization iterations to do
                		 doall
                		 (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[-0.5321583108529393 10.266428590549186] [7.2796862694825215 10.659164321881876] [7.540676463014467 0.029820907770213845] [5.689912136636714 4.618525995636671] [2.086297125152694 8.790352945538359]]
;;; :initial-log-Zs [-29.71820748217702 -103.56054242680013 -14.687282524874554 -30.15490917692295 -36.05234647979107]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.09768502061409035]
;;; :theta-best [7.540676463014467 0.029820907770213845]     :log-Z-theta-best -14.687282524874554     :mean-theta-best -14.778085847774804     :std-dev-theta-best 1.3435394025226302     :i-best 2
;;; :theta-next [6.549185062439843 1.0402443035631679]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -19.270393085340537
;;; :log-Z-i-best -14.687282524874554
;;; :theta-mean-best ([7.540676463014467 0.029820907770213845] -14.778085847774804)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.08917113503234315]
;;; :theta-best [7.540676463014467 0.029820907770213845]     :log-Z-theta-best -14.687282524874554     :mean-theta-best -14.692043917535216     :std-dev-theta-best 0.3716158465258613     :i-best 3
;;; :theta-next [8.250640857418134 0.8357346825930164]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -6.974173400192438
;;; :log-Z-i-best -14.687282524874554
;;; :theta-mean-best ([7.540676463014467 0.029820907770213845] -14.692043917535216)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.045581965406129986]
;;; :theta-best [8.250640857418134 0.8357346825930164]     :log-Z-theta-best -6.974173400192438     :mean-theta-best -6.987795988469244     :std-dev-theta-best 0.5294939061467452     :i-best 0
;;; :theta-next [8.892372550768464 3.0639989797937783]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.729884102513581
;;; :log-Z-i-best -6.974173400192438
;;; :theta-mean-best ([8.250640857418134 0.8357346825930164] -6.987795988469244)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04959324056291642]
;;; :theta-best [8.892372550768464 3.0639989797937783]     :log-Z-theta-best -2.729884102513581     :mean-theta-best -2.788161445262503     :std-dev-theta-best 1.0467183013532368     :i-best 0
;;; :theta-next [8.747376417020762 2.5735301562552193]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.890893224973918
;;; :log-Z-i-best -2.729884102513581
;;; :theta-mean-best ([8.892372550768464 3.0639989797937783] -2.788161445262503)
;;; :BO-Iteration 4
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.07074167353272263]
;;; :theta-best [8.892372550768464 3.0639989797937783]     :log-Z-theta-best -2.729884102513581     :mean-theta-best -2.7408772880451693     :std-dev-theta-best 0.27803434955296075     :i-best 1
;;; :theta-next [7.947791963413026 3.014562518316745]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -11.36140767069929
;;; :log-Z-i-best -2.729884102513581
;;; :theta-mean-best ([8.892372550768464 3.0639989797937783] -2.7408772880451693)
;;; :BO-Iteration 5
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.061567791688655144]
;;; :theta-best [8.892372550768464 3.0639989797937783]     :log-Z-theta-best -2.729884102513581     :mean-theta-best -2.7642485242290604     :std-dev-theta-best 0.3278627502964007     :i-best 2
;;; :theta-next [9.336827101176265 2.4311175523609476]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -0.4358601745809043
;;; :log-Z-i-best -2.729884102513581
;;; :theta-mean-best ([8.892372550768464 3.0639989797937783] -2.7642485242290604)
;;; :BO-Iteration 6
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.028385872511322296]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.5442380061303851     :std-dev-theta-best 1.4040878090828166     :i-best 0
;;; :theta-next [9.654576698288723 3.566782418176338]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -1.4444007895914464
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.5442380061303851)
;;; :BO-Iteration 7
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.024267952314446447]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.6297459106889534     :std-dev-theta-best 2.5745810028727987     :i-best 1
;;; :theta-next [9.603529601432438 2.657880426158515]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -0.5516658213340797
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.6297459106889534)
;;; :BO-Iteration 8
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.037067184139001556]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4416888965115362     :std-dev-theta-best 0.29814786566754203     :i-best 2
;;; :theta-next [2.801248553364741 0.046629185101557474]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -7.2425189294416805
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4416888965115362)
;;; :BO-Iteration 9
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.0482464999377279]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4605459972559629     :std-dev-theta-best 0.5916850440553203     :i-best 3
;;; :theta-next [3.7736774899790646 0.6617528979431472]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -3.6264435868744647
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4605459972559629)
;;; :BO-Iteration 10
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.019475387159281588]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.5000874271024003     :std-dev-theta-best 2.214215892894577     :i-best 2
;;; :theta-next [0.4689369352518131 1.3770959344312603]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.814392154688036
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.5000874271024003)
;;; :BO-Iteration 11
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.020067188700679216]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.47985530191910186     :std-dev-theta-best 0.7547760418299188     :i-best 5
;;; :theta-next [-2.8671468002468377 7.06503692829799]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -21.55206495033756
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.47985530191910186)
;;; :BO-Iteration 12
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04811614224232135]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.49297948385407153     :std-dev-theta-best 2.1762109157322853     :i-best 4
;;; :theta-next [9.87952908494879 1.7199104931510107]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.73189641452646
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.49297948385407153)
;;; :BO-Iteration 13
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03344168481088605]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.49997705400923564     :std-dev-theta-best 3.3475839282291666     :i-best 7
;;; :theta-next [3.413307501319391 2.156949413061427]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -0.7572773152605858
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.49997705400923564)
;;; :BO-Iteration 14
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.054241934178099964]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.43652489608261646     :std-dev-theta-best 0.23423883600909887     :i-best 8
;;; :theta-next [3.7627374857235 1.3760814547006857]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.4070653721433954
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.43652489608261646)
;;; :BO-Iteration 15
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.036120123919314566]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.39762800753061356     :std-dev-theta-best 1.4342947331798621     :i-best 7
;;; :theta-next [3.5555175620358686 4.738935045754122]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -8.851831082422095
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.39762800753061356)
;;; :BO-Iteration 16
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.027673839470871233]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.43604211584039376     :std-dev-theta-best 0.069224908625896     :i-best 10
;;; :theta-next [3.9566378836193916 2.6727013882479804]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -4.312295875910351
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.43604211584039376)
;;; :BO-Iteration 17
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.029709938800952226]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4794234074435053     :std-dev-theta-best 0.8006500605595425     :i-best 11
;;; :theta-next [8.961650580312412 0.8680301252604974]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.9569635797951808
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4794234074435053)
;;; :BO-Iteration 18
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.030454950403040252]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.7134297129247784     :std-dev-theta-best 2.9614411934950096     :i-best 12
;;; :theta-next [2.269615652372547 3.6644195047635186]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -4.196360052298528
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.7134297129247784)
;;; :BO-Iteration 19
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.0497684867711327]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.47371574198696464     :std-dev-theta-best 1.2348892420810516     :i-best 13
;;; :theta-next [0.736205816190056 5.033252315536716]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -17.133582412307298
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.47371574198696464)
;;; :BO-Iteration 20
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.017523687537137907]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4480813383535889     :std-dev-theta-best 0.3297625357369885     :i-best 14
;;; :theta-next [-4.87176543155768 9.684245875923098]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -62.43906171552868
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4480813383535889)
;;; :BO-Iteration 21
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04811310518836771]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.45415768643727006     :std-dev-theta-best 0.373918383775264     :i-best 15
;;; :theta-next [2.768293050763975 1.5985678769733556]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.030513198343817
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.45415768643727006)
;;; :BO-Iteration 22
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.012412424497497783]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.44475699410646996     :std-dev-theta-best 0.18679466903711214     :i-best 16
;;; :theta-next [-1.9834139265532094 14.280821785686726]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.456165027076064
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.44475699410646996)
;;; :BO-Iteration 23
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.030524251868782278]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.3521505765415469     :std-dev-theta-best 1.0379167355347558     :i-best 17
;;; :theta-next [2.9624449603820953 3.073178123486435]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -0.9796977246313592
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.3521505765415469)
;;; :BO-Iteration 24
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.034328333441257075]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -1.233597046114454     :std-dev-theta-best 5.622769901897606     :i-best 18
;;; :theta-next [8.882316366933786 1.8404190904825712]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -1.822594460503435
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -1.233597046114454)
;;; :BO-Iteration 25
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.025542249752886118]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.4503661750881065     :std-dev-theta-best 0.9619712497968259     :i-best 17
;;; :theta-next [2.5021612786951346 2.351581550436604]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.5204613986505358
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.4503661750881065)
;;; :BO-Iteration 26
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.012665953156203988]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.2016517350432565     :std-dev-theta-best 2.318731389770724     :i-best 18
;;; :theta-next [-1.2574270180721547 6.004905122286396]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -17.802679680806854
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.2016517350432565)
;;; :BO-Iteration 27
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.057589271776702596]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -6.868181524864667     :std-dev-theta-best 15.25236934808428     :i-best 21
;;; :theta-next [8.425781916495474 5.361926325681319]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -17.76871061766759
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -6.868181524864667)
;;; :BO-Iteration 28
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.018621638563558283]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4464944029087974     :std-dev-theta-best 0.33313001325227676     :i-best 22
;;; :theta-next [4.984842812110273 2.491834254886361]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -14.061074138039164
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4464944029087974)
;;; :BO-Iteration 29
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.01694731456060301]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4762564847550834     :std-dev-theta-best 0.7494280305512042     :i-best 23
;;; :theta-next [1.6426274310230644 14.70662404056983]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -129.70394022667256
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4762564847550834)
;;; :BO-Iteration 30
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.02539475237113733]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.44237183269669345     :std-dev-theta-best 0.14874046917599026     :i-best 24
;;; :theta-next [-0.7760832375258637 8.087869981886946]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -17.4531185294357
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.44237183269669345)
;;; :BO-Iteration 31
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.009556472245138026]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4484272737249171     :std-dev-theta-best 0.6946326275204951     :i-best 25
;;; :theta-next [2.4028969652731087 5.9122208534238005]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -11.84468985861266
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4484272737249171)
;;; :BO-Iteration 32
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.027651324099657912]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.48679922794413244     :std-dev-theta-best 0.8234582018533046     :i-best 24
;;; :theta-next [-4.522977436082874 13.750582778568942]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -12.56325183033696
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.48679922794413244)
;;; :BO-Iteration 33
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.020271597908412875]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4396144064978955     :std-dev-theta-best 0.4555754948300857     :i-best 27
;;; :theta-next [-3.920086158381068 14.801819164845323]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -3.4972076213955443
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4396144064978955)
;;; :BO-Iteration 34
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.015842893171715777]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.9018201708738047     :std-dev-theta-best 3.4971512232487054     :i-best 28
;;; :theta-next [-4.072156077215631 3.7740919393224654]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -121.96775219494351
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.9018201708738047)
;;; :BO-Iteration 35
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.019637185188151402]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.49494083381905796     :std-dev-theta-best 2.663007463266893     :i-best 29
;;; :theta-next [-3.0257547126556985 13.058317930063927]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -1.585776098543974
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.49494083381905796)
;;; :BO-Iteration 36
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03434253780362841]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.699722940680374     :std-dev-theta-best 3.5945906163672303     :i-best 30
;;; :theta-next [-2.4143178241043732 10.245209726658057]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.950040668819528
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.699722940680374)
;;; :BO-Iteration 37
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.08680784872111941]
;;; :theta-best [9.603529601432438 2.657880426158515]     :log-Z-theta-best -0.5516658213340797     :mean-theta-best -0.4462983581398561     :std-dev-theta-best 1.5431900125056977     :i-best 29
;;; :theta-next [-2.4581911289852316 11.634004907972994]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -3.4398121526168612
;;; :log-Z-i-best -0.5516658213340797
;;; :theta-mean-best ([9.603529601432438 2.657880426158515] -0.4462983581398561)
;;; :BO-Iteration 38
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.06297959434470314]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.44167584602081433     :std-dev-theta-best 0.8881616304881389     :i-best 32
;;; :theta-next [-3.640932804762329 11.87290827238454]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -4.241374105756842
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.44167584602081433)
;;; :BO-Iteration 39
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.016733317041225697]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.44683479855450514     :std-dev-theta-best 0.279317742558385     :i-best 33
;;; :theta-next [9.902103280185482 6.376258574309147]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -13.506433355717334
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.44683479855450514)
;;; :BO-Iteration 40
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.024393946580231966]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4736328790237394     :std-dev-theta-best 0.491016832522226     :i-best 34
;;; :theta-next [-2.66119661541875 10.77888287729697]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -1.6226870278224013
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4736328790237394)
;;; :BO-Iteration 41
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.046211699634938916]
;;; :theta-best [9.336827101176265 2.4311175523609476]     :log-Z-theta-best -0.4358601745809043     :mean-theta-best -0.4766235747647585     :std-dev-theta-best 1.6736345026841222     :i-best 35
;;; :theta-next [-3.193568981916109 12.533703620289828]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -0.42866172529558355
;;; :log-Z-i-best -0.4358601745809043
;;; :theta-mean-best ([9.336827101176265 2.4311175523609476] -0.4766235747647585)
;;; :BO-Iteration 42
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.024076806864073857]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.4258298860613934     :std-dev-theta-best 0.15760550273434204     :i-best 0
;;; :theta-next [3.0692144890207214 6.928170116951229]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -21.546689848438604
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.4258298860613934)
;;; :BO-Iteration 43
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.008453387914366742]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.3399367053864921     :std-dev-theta-best 1.1114228274894846     :i-best 1
;;; :theta-next [3.2137170433590505 3.6717284055696924]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -2.532036061488803
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.3399367053864921)
;;; :BO-Iteration 44
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.031176793199524207]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.42746750344637974     :std-dev-theta-best 0.11974896879692885     :i-best 2
;;; :theta-next [9.643315964636704 4.819911967401188]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -5.267712126772877
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.42746750344637974)
;;; :BO-Iteration 45
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.024756515473092218]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.42713846154124724     :std-dev-theta-best 0.836745977599468     :i-best 3
;;; :theta-next [9.683411759570932 0.4718545379259422]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -5.689927265698573
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.42713846154124724)
;;; :BO-Iteration 46
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.02803442499428095]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.41956075449485297     :std-dev-theta-best 0.6882057033936357     :i-best 4
;;; :theta-next [-3.1719186876317695 14.039383959531417]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -3.263083673585344
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.41956075449485297)
;;; :BO-Iteration 47
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.01774573837809687]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.4199436987798748     :std-dev-theta-best 0.15436036766642414     :i-best 5
;;; :theta-next [-2.0337005634452283 8.02468707011138]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -8.761918161269058
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.4199436987798748)
;;; :BO-Iteration 48
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.007481554467927541]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.3705122790284179     :std-dev-theta-best 0.6106777615331909     :i-best 6
;;; :theta-next [-3.135108804611513 11.697953418216882]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -0.7133374973138746
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.3705122790284179)
;;; :BO-Iteration 49
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.01515544117908262]
;;; :theta-best [-3.193568981916109 12.533703620289828]     :log-Z-theta-best -0.42866172529558355     :mean-theta-best -0.38397791092030786     :std-dev-theta-best 1.7861016771897575     :i-best 7
;;; :theta-next [1.485373037684476 2.2704987354964707]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -13.543333394010041
;;; :log-Z-i-best -0.42866172529558355
;;; :theta-mean-best ([-3.193568981916109 12.533703620289828] -0.38397791092030786)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.opt/branin-samples</span>","value":"#'worksheets.opt/branin-samples"}
;; <=

;; @@
branin-samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>7.540676463014467</span>","value":"7.540676463014467"},{"type":"html","content":"<span class='clj-double'>0.029820907770213845</span>","value":"0.029820907770213845"}],"value":"[7.540676463014467 0.029820907770213845]"},{"type":"html","content":"<span class='clj-double'>-14.778085847774804</span>","value":"-14.778085847774804"}],"value":"([7.540676463014467 0.029820907770213845] -14.778085847774804)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>7.540676463014467</span>","value":"7.540676463014467"},{"type":"html","content":"<span class='clj-double'>0.029820907770213845</span>","value":"0.029820907770213845"}],"value":"[7.540676463014467 0.029820907770213845]"},{"type":"html","content":"<span class='clj-double'>-14.692043917535216</span>","value":"-14.692043917535216"}],"value":"([7.540676463014467 0.029820907770213845] -14.692043917535216)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>8.250640857418134</span>","value":"8.250640857418134"},{"type":"html","content":"<span class='clj-double'>0.8357346825930164</span>","value":"0.8357346825930164"}],"value":"[8.250640857418134 0.8357346825930164]"},{"type":"html","content":"<span class='clj-double'>-6.987795988469244</span>","value":"-6.987795988469244"}],"value":"([8.250640857418134 0.8357346825930164] -6.987795988469244)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>8.892372550768464</span>","value":"8.892372550768464"},{"type":"html","content":"<span class='clj-double'>3.0639989797937783</span>","value":"3.0639989797937783"}],"value":"[8.892372550768464 3.0639989797937783]"},{"type":"html","content":"<span class='clj-double'>-2.788161445262503</span>","value":"-2.788161445262503"}],"value":"([8.892372550768464 3.0639989797937783] -2.788161445262503)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>8.892372550768464</span>","value":"8.892372550768464"},{"type":"html","content":"<span class='clj-double'>3.0639989797937783</span>","value":"3.0639989797937783"}],"value":"[8.892372550768464 3.0639989797937783]"},{"type":"html","content":"<span class='clj-double'>-2.7408772880451693</span>","value":"-2.7408772880451693"}],"value":"([8.892372550768464 3.0639989797937783] -2.7408772880451693)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>8.892372550768464</span>","value":"8.892372550768464"},{"type":"html","content":"<span class='clj-double'>3.0639989797937783</span>","value":"3.0639989797937783"}],"value":"[8.892372550768464 3.0639989797937783]"},{"type":"html","content":"<span class='clj-double'>-2.7642485242290604</span>","value":"-2.7642485242290604"}],"value":"([8.892372550768464 3.0639989797937783] -2.7642485242290604)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.5442380061303851</span>","value":"-0.5442380061303851"}],"value":"([9.336827101176265 2.4311175523609476] -0.5442380061303851)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.6297459106889534</span>","value":"-0.6297459106889534"}],"value":"([9.336827101176265 2.4311175523609476] -0.6297459106889534)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4416888965115362</span>","value":"-0.4416888965115362"}],"value":"([9.336827101176265 2.4311175523609476] -0.4416888965115362)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4605459972559629</span>","value":"-0.4605459972559629"}],"value":"([9.336827101176265 2.4311175523609476] -0.4605459972559629)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.5000874271024003</span>","value":"-0.5000874271024003"}],"value":"([9.603529601432438 2.657880426158515] -0.5000874271024003)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.47985530191910186</span>","value":"-0.47985530191910186"}],"value":"([9.336827101176265 2.4311175523609476] -0.47985530191910186)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.49297948385407153</span>","value":"-0.49297948385407153"}],"value":"([9.603529601432438 2.657880426158515] -0.49297948385407153)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.49997705400923564</span>","value":"-0.49997705400923564"}],"value":"([9.336827101176265 2.4311175523609476] -0.49997705400923564)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.43652489608261646</span>","value":"-0.43652489608261646"}],"value":"([9.336827101176265 2.4311175523609476] -0.43652489608261646)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.39762800753061356</span>","value":"-0.39762800753061356"}],"value":"([9.603529601432438 2.657880426158515] -0.39762800753061356)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.43604211584039376</span>","value":"-0.43604211584039376"}],"value":"([9.336827101176265 2.4311175523609476] -0.43604211584039376)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4794234074435053</span>","value":"-0.4794234074435053"}],"value":"([9.336827101176265 2.4311175523609476] -0.4794234074435053)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.7134297129247784</span>","value":"-0.7134297129247784"}],"value":"([9.336827101176265 2.4311175523609476] -0.7134297129247784)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.47371574198696464</span>","value":"-0.47371574198696464"}],"value":"([9.336827101176265 2.4311175523609476] -0.47371574198696464)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4480813383535889</span>","value":"-0.4480813383535889"}],"value":"([9.336827101176265 2.4311175523609476] -0.4480813383535889)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.45415768643727006</span>","value":"-0.45415768643727006"}],"value":"([9.336827101176265 2.4311175523609476] -0.45415768643727006)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.44475699410646996</span>","value":"-0.44475699410646996"}],"value":"([9.336827101176265 2.4311175523609476] -0.44475699410646996)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.3521505765415469</span>","value":"-0.3521505765415469"}],"value":"([9.336827101176265 2.4311175523609476] -0.3521505765415469)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-1.233597046114454</span>","value":"-1.233597046114454"}],"value":"([9.336827101176265 2.4311175523609476] -1.233597046114454)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.4503661750881065</span>","value":"-0.4503661750881065"}],"value":"([9.603529601432438 2.657880426158515] -0.4503661750881065)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.2016517350432565</span>","value":"-0.2016517350432565"}],"value":"([9.603529601432438 2.657880426158515] -0.2016517350432565)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-6.868181524864667</span>","value":"-6.868181524864667"}],"value":"([9.336827101176265 2.4311175523609476] -6.868181524864667)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4464944029087974</span>","value":"-0.4464944029087974"}],"value":"([9.336827101176265 2.4311175523609476] -0.4464944029087974)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4762564847550834</span>","value":"-0.4762564847550834"}],"value":"([9.336827101176265 2.4311175523609476] -0.4762564847550834)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.44237183269669345</span>","value":"-0.44237183269669345"}],"value":"([9.336827101176265 2.4311175523609476] -0.44237183269669345)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4484272737249171</span>","value":"-0.4484272737249171"}],"value":"([9.336827101176265 2.4311175523609476] -0.4484272737249171)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.48679922794413244</span>","value":"-0.48679922794413244"}],"value":"([9.603529601432438 2.657880426158515] -0.48679922794413244)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4396144064978955</span>","value":"-0.4396144064978955"}],"value":"([9.336827101176265 2.4311175523609476] -0.4396144064978955)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.9018201708738047</span>","value":"-0.9018201708738047"}],"value":"([9.336827101176265 2.4311175523609476] -0.9018201708738047)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.49494083381905796</span>","value":"-0.49494083381905796"}],"value":"([9.336827101176265 2.4311175523609476] -0.49494083381905796)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.699722940680374</span>","value":"-0.699722940680374"}],"value":"([9.336827101176265 2.4311175523609476] -0.699722940680374)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.603529601432438</span>","value":"9.603529601432438"},{"type":"html","content":"<span class='clj-double'>2.657880426158515</span>","value":"2.657880426158515"}],"value":"[9.603529601432438 2.657880426158515]"},{"type":"html","content":"<span class='clj-double'>-0.4462983581398561</span>","value":"-0.4462983581398561"}],"value":"([9.603529601432438 2.657880426158515] -0.4462983581398561)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.44167584602081433</span>","value":"-0.44167584602081433"}],"value":"([9.336827101176265 2.4311175523609476] -0.44167584602081433)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.44683479855450514</span>","value":"-0.44683479855450514"}],"value":"([9.336827101176265 2.4311175523609476] -0.44683479855450514)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4736328790237394</span>","value":"-0.4736328790237394"}],"value":"([9.336827101176265 2.4311175523609476] -0.4736328790237394)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>9.336827101176265</span>","value":"9.336827101176265"},{"type":"html","content":"<span class='clj-double'>2.4311175523609476</span>","value":"2.4311175523609476"}],"value":"[9.336827101176265 2.4311175523609476]"},{"type":"html","content":"<span class='clj-double'>-0.4766235747647585</span>","value":"-0.4766235747647585"}],"value":"([9.336827101176265 2.4311175523609476] -0.4766235747647585)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.4258298860613934</span>","value":"-0.4258298860613934"}],"value":"([-3.193568981916109 12.533703620289828] -0.4258298860613934)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.3399367053864921</span>","value":"-0.3399367053864921"}],"value":"([-3.193568981916109 12.533703620289828] -0.3399367053864921)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.42746750344637974</span>","value":"-0.42746750344637974"}],"value":"([-3.193568981916109 12.533703620289828] -0.42746750344637974)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.42713846154124724</span>","value":"-0.42713846154124724"}],"value":"([-3.193568981916109 12.533703620289828] -0.42713846154124724)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.41956075449485297</span>","value":"-0.41956075449485297"}],"value":"([-3.193568981916109 12.533703620289828] -0.41956075449485297)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.4199436987798748</span>","value":"-0.4199436987798748"}],"value":"([-3.193568981916109 12.533703620289828] -0.4199436987798748)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.3705122790284179</span>","value":"-0.3705122790284179"}],"value":"([-3.193568981916109 12.533703620289828] -0.3705122790284179)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-3.193568981916109</span>","value":"-3.193568981916109"},{"type":"html","content":"<span class='clj-double'>12.533703620289828</span>","value":"12.533703620289828"}],"value":"[-3.193568981916109 12.533703620289828]"},{"type":"html","content":"<span class='clj-double'>-0.38397791092030786</span>","value":"-0.38397791092030786"}],"value":"([-3.193568981916109 12.533703620289828] -0.38397791092030786)"}],"value":"[([7.540676463014467 0.029820907770213845] -14.778085847774804) ([7.540676463014467 0.029820907770213845] -14.692043917535216) ([8.250640857418134 0.8357346825930164] -6.987795988469244) ([8.892372550768464 3.0639989797937783] -2.788161445262503) ([8.892372550768464 3.0639989797937783] -2.7408772880451693) ([8.892372550768464 3.0639989797937783] -2.7642485242290604) ([9.336827101176265 2.4311175523609476] -0.5442380061303851) ([9.336827101176265 2.4311175523609476] -0.6297459106889534) ([9.336827101176265 2.4311175523609476] -0.4416888965115362) ([9.336827101176265 2.4311175523609476] -0.4605459972559629) ([9.603529601432438 2.657880426158515] -0.5000874271024003) ([9.336827101176265 2.4311175523609476] -0.47985530191910186) ([9.603529601432438 2.657880426158515] -0.49297948385407153) ([9.336827101176265 2.4311175523609476] -0.49997705400923564) ([9.336827101176265 2.4311175523609476] -0.43652489608261646) ([9.603529601432438 2.657880426158515] -0.39762800753061356) ([9.336827101176265 2.4311175523609476] -0.43604211584039376) ([9.336827101176265 2.4311175523609476] -0.4794234074435053) ([9.336827101176265 2.4311175523609476] -0.7134297129247784) ([9.336827101176265 2.4311175523609476] -0.47371574198696464) ([9.336827101176265 2.4311175523609476] -0.4480813383535889) ([9.336827101176265 2.4311175523609476] -0.45415768643727006) ([9.336827101176265 2.4311175523609476] -0.44475699410646996) ([9.336827101176265 2.4311175523609476] -0.3521505765415469) ([9.336827101176265 2.4311175523609476] -1.233597046114454) ([9.603529601432438 2.657880426158515] -0.4503661750881065) ([9.603529601432438 2.657880426158515] -0.2016517350432565) ([9.336827101176265 2.4311175523609476] -6.868181524864667) ([9.336827101176265 2.4311175523609476] -0.4464944029087974) ([9.336827101176265 2.4311175523609476] -0.4762564847550834) ([9.336827101176265 2.4311175523609476] -0.44237183269669345) ([9.336827101176265 2.4311175523609476] -0.4484272737249171) ([9.603529601432438 2.657880426158515] -0.48679922794413244) ([9.336827101176265 2.4311175523609476] -0.4396144064978955) ([9.336827101176265 2.4311175523609476] -0.9018201708738047) ([9.336827101176265 2.4311175523609476] -0.49494083381905796) ([9.336827101176265 2.4311175523609476] -0.699722940680374) ([9.603529601432438 2.657880426158515] -0.4462983581398561) ([9.336827101176265 2.4311175523609476] -0.44167584602081433) ([9.336827101176265 2.4311175523609476] -0.44683479855450514) ([9.336827101176265 2.4311175523609476] -0.4736328790237394) ([9.336827101176265 2.4311175523609476] -0.4766235747647585) ([-3.193568981916109 12.533703620289828] -0.4258298860613934) ([-3.193568981916109 12.533703620289828] -0.3399367053864921) ([-3.193568981916109 12.533703620289828] -0.42746750344637974) ([-3.193568981916109 12.533703620289828] -0.42713846154124724) ([-3.193568981916109 12.533703620289828] -0.41956075449485297) ([-3.193568981916109 12.533703620289828] -0.4199436987798748) ([-3.193568981916109 12.533703620289828] -0.3705122790284179) ([-3.193568981916109 12.533703620289828] -0.38397791092030786)]"}
;; <=

;; **
;;; ## Hartmann 6D
;; **

;; **
;;; @@\begin{align}
;;; f(x) &= -\sum\_{i=1}^{4} \alpha\_i \exp \left(-\sum\_{j=1}^{6} A\_{i,j}(x\_j-P\_{i,j})^2\right), \quad \rm{where} \\\\
;;; \alpha &= (1,1.2,3,3.2)^T \\\\
;;; A &= \left(\begin{matrix}
;;; 10 & 3 & 17 & 3.5 & 1.6 & 8 \\\\
;;; 0.05 & 10 & 17 & 0.1 & 8 & 14 \\\\
;;; 3 & 3.5 & 1.7 & 10 & 17 & 8 \\\\
;;; 17 & 8 & 0.05 & 10 & 0.1 & 14
;;; \end{matrix}\right) \\\\
;;; P &= 10^{-4} \left(
;;; \begin{matrix}
;;; 1312 & 1696 & 5569 & 124 & 8283 & 5886\\\\
;;; 2329 & 4135 & 8307 & 3736 & 1004 & 9991 \\\\
;;; 2358 & 1451 & 3522 & 2883 & 3047 & 6650 \\\\
;;; 4047 & 8828 & 8732 & 5743 & 1091 & 381
;;; \end{matrix}
;;; \right)
;;; \end{align}@@
;;; 
;;; Bounds @@0 \le x\_i \le 1, \; \forall i \in \\{1,\dots,6\\}@@
;;; 
;;; The function has 6 local minima and one global minima at
;;; @@\begin{align}
;;; x^\* &= (0.20169,0.150011,0.476874,0.275332,0.311652,0.6573),
;;; \end{align}@@
;;; with @@f(x^\*) = -3.32237@@.  
;; **

;; @@
(def A
  [[10 3 17 3.5 1.7 8]
   [0.05 10 17 0.1 8 14]
   [3 3.5 1.7 10 17 8]
   [17 8 0.05 10 0.1 14]])

(def P
  [[0.1312 0.1696 0.5569 0.0124 0.8283 0.5886]
   [0.2329 0.4135 0.8307 0.3736 0.1004 0.9991]
   [0.2348 0.1451 0.3522 0.2883 0.3047 0.6650]
   [0.4047 0.8828 0.8732 0.5743 0.1091 0.0381]])

(def alpha
  [1 1.2 3 3.2])

(defn hartman-6d [x1 x2 x3 x4 x5 x6]
  (let [x [x1 x2 x3 x4 x5 x6]
        dxP (mat/pow
             (mat/sub P x)
             2)
        terms (mat/transpose
               (mat/mul A dxP))
        terms-reduced (reduce mat/add terms)
        exp-terms (mat/exp
                   (mat/sub 0 terms-reduced))
        f (- 0
             (reduce + (mat/mul alpha exp-terms)))]
    (- f)))

(with-primitive-procedures [factor hartman-6d]
 (defopt h6 [] [x1 x2 x3 x4 x5 x6]
   (let [x1 (sample (uniform-continuous 0 1))
         x2 (sample (uniform-continuous 0 1))
         x3 (sample (uniform-continuous 0 1))
         x4 (sample (uniform-continuous 0 1))
         x5 (sample (uniform-continuous 0 1))
         x6 (sample (uniform-continuous 0 1))
         f (hartman-6d x1 x2 x3 x4 x5 x6)
         prior-correction 0]
     (observe (factor)
              (+ f prior-correction)))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;worksheets.opt/h6</span>","value":"#'worksheets.opt/h6"}
;; <=

;; @@

;; @@
