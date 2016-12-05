;; gorilla-repl.fileformat = 1

;; **
;;; # Simple Bimodal Modal
;; **

;; @@
(ns examples.simple-bimodal
  (require [bopp.core :refer :all]
           [anglican.runtime :refer :all]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; In this simple example, we demonstrate the ability of BOPP to carry out unbounded optimization on a 1D problem with a significant prior-posterior mismatch, namely @@p \left(\theta\right)={\rm Normal}(0, 0.5)@@ and @@p \left(Y|\theta\right)={\rm Normal}(5-\left|\theta\right|,0.5)@@ where we want to optimize @@p(\theta | Y)@@ with respec to @@\theta@@.  This requires BOPP to adapt to the target and effectively establish maxima in the presence of multiple modes.  It should be noted that this is a deterministic optimization problem.
;; **

;; @@
(defopt simple-bimodal [y] [theta]
  (let [theta (sample (normal 0 0.5))]
    (observe
     (normal (sqrt (* theta theta)) 0.5) y)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;examples.simple-bimodal/simple-bimodal</span>","value":"#'examples.simple-bimodal/simple-bimodal"}
;; <=

;; @@
(def samples (->> (doopt :importance 
                         simple-bimodal
                         [5]
                         1 ;; Number of particles
                         :bo-verbose true) ;; So that BOPP spits out some information to follow its progress
                (take 50) ;; Number of optimization iterations to do
                doall
                (mapv #(take 2 %))))
;; @@
;; ->
;;; :initial-thetas [[0.2505956957899326] [0.11086693241781502] [-1.4248822047652165] [1.351404150830059] [-0.4886839333683613]]
;;; :initial-log-Zs [-45.69086160048455 -48.283409963748106 -30.07509579981146 -30.72867240421112 -41.63315198485172]
;;; :BO-Iteration 0
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.031707112310298075]
;;; :theta-best [-1.4248822047652165]     :log-Z-theta-best -30.07509579981146     :mean-theta-best -30.083381179578573     :std-dev-theta-best 0.1571199437592305     :i-best 2
;;; :theta-next [-1.2418512774120993]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.783335537887723
;;; :log-Z-i-best -30.07509579981146
;;; :theta-mean-best ([-1.4248822047652165] -30.083381179578573)
;;; :BO-Iteration 1
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.060088814436476595]
;;; :theta-best [-1.4248822047652165]     :log-Z-theta-best -30.07509579981146     :mean-theta-best -30.12782126743909     :std-dev-theta-best 0.45481953479236203     :i-best 3
;;; :theta-next [-1.543133703807011]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.11395514044981
;;; :log-Z-i-best -30.07509579981146
;;; :theta-mean-best ([-1.4248822047652165] -30.12782126743909)
;;; :BO-Iteration 2
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.03066849454828423]
;;; :theta-best [-1.543133703807011]     :log-Z-theta-best -29.11395514044981     :mean-theta-best -29.241429564881507     :std-dev-theta-best 0.5118576922002365     :i-best 0
;;; :theta-next [1.1056085827767945]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.228892402992415
;;; :log-Z-i-best -29.11395514044981
;;; :theta-mean-best ([-1.543133703807011] -29.241429564881507)
;;; :BO-Iteration 3
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.017147069334507752]
;;; :theta-best [-1.543133703807011]     :log-Z-theta-best -29.11395514044981     :mean-theta-best -29.1420594612676     :std-dev-theta-best 0.2142032351988405     :i-best 1
;;; :theta-next [1.517643767574858]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.311677774827533
;;; :log-Z-i-best -29.11395514044981
;;; :theta-mean-best ([-1.543133703807011] -29.1420594612676)
;;; :BO-Iteration 4
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.013668959223999074]
;;; :theta-best [-1.543133703807011]     :log-Z-theta-best -29.11395514044981     :mean-theta-best -29.56356228288295     :std-dev-theta-best 1.15748224015815     :i-best 2
;;; :theta-next [1.441959188995377]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.929384136294736
;;; :log-Z-i-best -29.11395514044981
;;; :theta-mean-best ([-1.543133703807011] -29.56356228288295)
;;; :BO-Iteration 5
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.01349572760627439]
;;; :theta-best [-1.543133703807011]     :log-Z-theta-best -29.11395514044981     :mean-theta-best -29.678762847408148     :std-dev-theta-best 1.5956394118349095     :i-best 3
;;; :theta-next [1.3906072641457308]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.374591674754335
;;; :log-Z-i-best -29.11395514044981
;;; :theta-mean-best ([-1.543133703807011] -29.678762847408148)
;;; :BO-Iteration 6
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.05524002645797168]
;;; :theta-best [-1.543133703807011]     :log-Z-theta-best -29.11395514044981     :mean-theta-best -29.122961630797775     :std-dev-theta-best 0.1502936144692005     :i-best 4
;;; :theta-next [-1.6717662182040178]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.19546749452175
;;; :log-Z-i-best -29.11395514044981
;;; :theta-mean-best ([-1.543133703807011] -29.122961630797775)
;;; :BO-Iteration 7
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.04241269381684206]
;;; :theta-best [-1.6717662182040178]     :log-Z-theta-best -28.19546749452175     :mean-theta-best -28.480068501778426     :std-dev-theta-best 0.801583123057516     :i-best 0
;;; :theta-next [2.043442272659671]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -26.285362538866117
;;; :log-Z-i-best -28.19546749452175
;;; :theta-mean-best ([-1.6717662182040178] -28.480068501778426)
;;; :BO-Iteration 8
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.7106270834564292E-4]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.28769463307808     :std-dev-theta-best 0.10093246236347242     :i-best 0
;;; :theta-next [-0.9444332283485077]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -35.130734629554425
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.28769463307808)
;;; :BO-Iteration 9
;;; :n-gps-in-acq-function 20
;;; :acq-opt [5.741982594817627E-5]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.340927493309696     :std-dev-theta-best 0.4929742121183837     :i-best 1
;;; :theta-next [-1.6682722836069186]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.21866668215526
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.340927493309696)
;;; :BO-Iteration 10
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.016947485672542822]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.302408990460034     :std-dev-theta-best 0.24323657085154807     :i-best 2
;;; :theta-next [-2.0305175460807554]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -26.333237803441598
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.302408990460034)
;;; :BO-Iteration 11
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.001811122306550161]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.29135508269127     :std-dev-theta-best 0.18856444216650975     :i-best 3
;;; :theta-next [0.7439927936305826]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -37.78582794057476
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.29135508269127)
;;; :BO-Iteration 12
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.1657368761487604E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.353294691878464     :std-dev-theta-best 0.6491923278178224     :i-best 4
;;; :theta-next [-1.501337363982227]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.440890947601325
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.353294691878464)
;;; :BO-Iteration 13
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.4620023509001614E-6]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.2934672021875     :std-dev-theta-best 0.17145686330892923     :i-best 5
;;; :theta-next [-0.7043944308225996]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -38.34838014553304
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.2934672021875)
;;; :BO-Iteration 14
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.7543030285151716E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.332084339700987     :std-dev-theta-best 0.3452588207627258     :i-best 6
;;; :theta-next [1.258775328729079]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.614137443575878
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.332084339700987)
;;; :BO-Iteration 15
;;; :n-gps-in-acq-function 20
;;; :acq-opt [8.90738630847298E-5]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.572398822979498     :std-dev-theta-best 1.0657727126291616     :i-best 7
;;; :theta-next [0.534651010222146]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -40.90196931177278
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.572398822979498)
;;; :BO-Iteration 16
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.4800927472587289E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.292886175950795     :std-dev-theta-best 0.17560134836421554     :i-best 8
;;; :theta-next [-1.0491506600755434]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.87143793392638
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.292886175950795)
;;; :BO-Iteration 17
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.8977550348736982E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.287047734809185     :std-dev-theta-best 0.0749079401515374     :i-best 9
;;; :theta-next [-1.3882387508460883]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.395635005770522
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.287047734809185)
;;; :BO-Iteration 18
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.251899173607304E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.31591192840172     :std-dev-theta-best 0.3971722702271094     :i-best 10
;;; :theta-next [0.9113146154714733]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -35.54726770934766
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.31591192840172)
;;; :BO-Iteration 19
;;; :n-gps-in-acq-function 20
;;; :acq-opt [4.0846402008468347E-4]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.340092699623156     :std-dev-theta-best 0.5315096677830029     :i-best 11
;;; :theta-next [1.5956798398160625]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.722762513749863
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.340092699623156)
;;; :BO-Iteration 20
;;; :n-gps-in-acq-function 20
;;; :acq-opt [3.9154088744329226E-5]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.29383329113715     :std-dev-theta-best 0.12504496797844536     :i-best 12
;;; :theta-next [-1.7583962790686707]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.651487020886226
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.29383329113715)
;;; :BO-Iteration 21
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.5923435778614897E-4]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.299724433695104     :std-dev-theta-best 0.23256012520867758     :i-best 13
;;; :theta-next [1.6714247806178113]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.197730281986424
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.299724433695104)
;;; :BO-Iteration 22
;;; :n-gps-in-acq-function 20
;;; :acq-opt [6.945223635503506E-5]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.46301448642859     :std-dev-theta-best 0.8575755792372701     :i-best 14
;;; :theta-next [-1.5827977626974836]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.81662248174042
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.46301448642859)
;;; :BO-Iteration 23
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.016381284405726E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.304765199318076     :std-dev-theta-best 0.21305498813714988     :i-best 15
;;; :theta-next [1.4664593373827097]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.724407910423004
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.304765199318076)
;;; :BO-Iteration 24
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.001368021617077002]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.31620451698342     :std-dev-theta-best 0.38401996721073156     :i-best 16
;;; :theta-next [1.7249502414111055]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.85439121844427
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.31620451698342)
;;; :BO-Iteration 25
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.7743520040878977E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.29651147062012     :std-dev-theta-best 0.17507615359460538     :i-best 17
;;; :theta-next [1.3055263131411428]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.158652259682217
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.29651147062012)
;;; :BO-Iteration 26
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.4277196166603447E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.285774668867358     :std-dev-theta-best 0.03223042601952019     :i-best 18
;;; :theta-next [-1.0850016699506226]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.46046380145956
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.285774668867358)
;;; :BO-Iteration 27
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.0697662109432795E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.289518753348204     :std-dev-theta-best 0.11954087599823496     :i-best 19
;;; :theta-next [-1.5063031831436062]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.401316160611575
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.289518753348204)
;;; :BO-Iteration 28
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.7520921751626475E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.321399060018408     :std-dev-theta-best 0.3525748455439477     :i-best 20
;;; :theta-next [-1.283437387807642]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.371681062826628
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.321399060018408)
;;; :BO-Iteration 29
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.1029300367851674E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.291719012171818     :std-dev-theta-best 0.12004293036446144     :i-best 21
;;; :theta-next [1.5304696136889024]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.211539385211637
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.291719012171818)
;;; :BO-Iteration 30
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.1285144763415508E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.291642014843802     :std-dev-theta-best 0.1372225369048798     :i-best 22
;;; :theta-next [0.8608765240025393]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -36.198485783553444
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.291642014843802)
;;; :BO-Iteration 31
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.6838602840197635E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.3870495943686     :std-dev-theta-best 0.5063945748569144     :i-best 23
;;; :theta-next [-1.220129860194312]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -32.003853004354376
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.3870495943686)
;;; :BO-Iteration 32
;;; :n-gps-in-acq-function 20
;;; :acq-opt [9.09729490311332E-4]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.840389579270344     :std-dev-theta-best 1.7217975958453342     :i-best 24
;;; :theta-next [1.416993317729938]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.143196600655884
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.840389579270344)
;;; :BO-Iteration 33
;;; :n-gps-in-acq-function 20
;;; :acq-opt [3.943281391948341E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.373344049215163     :std-dev-theta-best 0.5008042084823106     :i-best 25
;;; :theta-next [1.6312426890883591]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -28.470539766339158
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.373344049215163)
;;; :BO-Iteration 34
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.645259167847283E-7]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.545867889168637     :std-dev-theta-best 1.073343038265487     :i-best 26
;;; :theta-next [1.4007871389482425]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.28465836089582
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.545867889168637)
;;; :BO-Iteration 35
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.52882115031281E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.285533043921927     :std-dev-theta-best 0.016041403648592963     :i-best 27
;;; :theta-next [1.801518553173016]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.403088031535525
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.285533043921927)
;;; :BO-Iteration 36
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.005047157546226585]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.29185864191977     :std-dev-theta-best 0.10879085784431673     :i-best 28
;;; :theta-next [-1.8508404422256413]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -27.137215231088657
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.29185864191977)
;;; :BO-Iteration 37
;;; :n-gps-in-acq-function 15
;;; :acq-opt [1.1616119586788921E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.28634495850172     :std-dev-theta-best 0.04427928331489927     :i-best 29
;;; :theta-next [-0.8825744456489]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -35.91584440076131
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.28634495850172)
;;; :BO-Iteration 38
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.6956533225346757E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.333293939614407     :std-dev-theta-best 0.2825569788540448     :i-best 30
;;; :theta-next [-1.25088791174527]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.692706741385823
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.333293939614407)
;;; :BO-Iteration 39
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.6495253369238775E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.290124288248105     :std-dev-theta-best 0.1176228163653029     :i-best 31
;;; :theta-next [-1.2351728193614593]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -31.850733892817615
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.290124288248105)
;;; :BO-Iteration 40
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.530916530978197E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.351358810316565     :std-dev-theta-best 0.48270622053985635     :i-best 32
;;; :theta-next [-1.1304035122439182]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -32.95476086238463
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.351358810316565)
;;; :BO-Iteration 41
;;; :n-gps-in-acq-function 20
;;; :acq-opt [0.011337178712543057]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.364149621608078     :std-dev-theta-best 0.3907213064700095     :i-best 33
;;; :theta-next [2.0106036481706226]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -26.409617862025065
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.364149621608078)
;;; :BO-Iteration 42
;;; :n-gps-in-acq-function 18
;;; :acq-opt [1.9516005238480018E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.29116825482903     :std-dev-theta-best 0.05341378153054276     :i-best 34
;;; :theta-next [-1.434951936500747]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.988892215543487
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.29116825482903)
;;; :BO-Iteration 43
;;; :n-gps-in-acq-function 20
;;; :acq-opt [5.67258936263066E-5]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.38294371653386     :std-dev-theta-best 0.46475659573933487     :i-best 35
;;; :theta-next [-0.23797369118443168]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -45.918634792384594
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.38294371653386)
;;; :BO-Iteration 44
;;; :n-gps-in-acq-function 17
;;; :acq-opt [1.675710679752486E-5]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.548672603788525     :std-dev-theta-best 0.9348027605558952     :i-best 36
;;; :theta-next [1.49081005784819]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.52544006265095
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.548672603788525)
;;; :BO-Iteration 45
;;; :n-gps-in-acq-function 20
;;; :acq-opt [2.0537268034976214E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.31471502276748     :std-dev-theta-best 0.24000037857226408     :i-best 37
;;; :theta-next [-1.4839550996345854]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -29.580971663523716
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.31471502276748)
;;; :BO-Iteration 46
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.4483874590829698E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.294373520297267     :std-dev-theta-best 0.12917007837412448     :i-best 38
;;; :theta-next [1.086803483289433]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -33.44008028466097
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.294373520297267)
;;; :BO-Iteration 47
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.828654738486437E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.301562571569303     :std-dev-theta-best 0.1594773109153447     :i-best 39
;;; :theta-next [1.3442630308822188]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.794494472431676
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.301562571569303)
;;; :BO-Iteration 48
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.9634493706602717E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.3233540892188     :std-dev-theta-best 0.22174457947446377     :i-best 40
;;; :theta-next [1.4314418787066017]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -30.018848539617764
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.3233540892188)
;;; :BO-Iteration 49
;;; :n-gps-in-acq-function 20
;;; :acq-opt [1.6285707122345373E-8]
;;; :theta-best [2.043442272659671]     :log-Z-theta-best -26.285362538866117     :mean-theta-best -26.321426937865628     :std-dev-theta-best 0.3015235055771824     :i-best 41
;;; :theta-next [1.192884795856358]
;;; Calling original query with theta next  
;;; :log-Z-theta-next -32.28578333290336
;;; :log-Z-i-best -26.285362538866117
;;; :theta-mean-best ([2.043442272659671] -26.321426937865628)
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;examples.simple-bimodal/samples</span>","value":"#'examples.simple-bimodal/samples"}
;; <=

;; @@
samples
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.4248822047652165</span>","value":"-1.4248822047652165"}],"value":"[-1.4248822047652165]"},{"type":"html","content":"<span class='clj-double'>-30.083381179578573</span>","value":"-30.083381179578573"}],"value":"([-1.4248822047652165] -30.083381179578573)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.4248822047652165</span>","value":"-1.4248822047652165"}],"value":"[-1.4248822047652165]"},{"type":"html","content":"<span class='clj-double'>-30.12782126743909</span>","value":"-30.12782126743909"}],"value":"([-1.4248822047652165] -30.12782126743909)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.543133703807011</span>","value":"-1.543133703807011"}],"value":"[-1.543133703807011]"},{"type":"html","content":"<span class='clj-double'>-29.241429564881507</span>","value":"-29.241429564881507"}],"value":"([-1.543133703807011] -29.241429564881507)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.543133703807011</span>","value":"-1.543133703807011"}],"value":"[-1.543133703807011]"},{"type":"html","content":"<span class='clj-double'>-29.1420594612676</span>","value":"-29.1420594612676"}],"value":"([-1.543133703807011] -29.1420594612676)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.543133703807011</span>","value":"-1.543133703807011"}],"value":"[-1.543133703807011]"},{"type":"html","content":"<span class='clj-double'>-29.56356228288295</span>","value":"-29.56356228288295"}],"value":"([-1.543133703807011] -29.56356228288295)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.543133703807011</span>","value":"-1.543133703807011"}],"value":"[-1.543133703807011]"},{"type":"html","content":"<span class='clj-double'>-29.678762847408148</span>","value":"-29.678762847408148"}],"value":"([-1.543133703807011] -29.678762847408148)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.543133703807011</span>","value":"-1.543133703807011"}],"value":"[-1.543133703807011]"},{"type":"html","content":"<span class='clj-double'>-29.122961630797775</span>","value":"-29.122961630797775"}],"value":"([-1.543133703807011] -29.122961630797775)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-1.6717662182040178</span>","value":"-1.6717662182040178"}],"value":"[-1.6717662182040178]"},{"type":"html","content":"<span class='clj-double'>-28.480068501778426</span>","value":"-28.480068501778426"}],"value":"([-1.6717662182040178] -28.480068501778426)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.28769463307808</span>","value":"-26.28769463307808"}],"value":"([2.043442272659671] -26.28769463307808)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.340927493309696</span>","value":"-26.340927493309696"}],"value":"([2.043442272659671] -26.340927493309696)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.302408990460034</span>","value":"-26.302408990460034"}],"value":"([2.043442272659671] -26.302408990460034)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.29135508269127</span>","value":"-26.29135508269127"}],"value":"([2.043442272659671] -26.29135508269127)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.353294691878464</span>","value":"-26.353294691878464"}],"value":"([2.043442272659671] -26.353294691878464)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.2934672021875</span>","value":"-26.2934672021875"}],"value":"([2.043442272659671] -26.2934672021875)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.332084339700987</span>","value":"-26.332084339700987"}],"value":"([2.043442272659671] -26.332084339700987)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.572398822979498</span>","value":"-26.572398822979498"}],"value":"([2.043442272659671] -26.572398822979498)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.292886175950795</span>","value":"-26.292886175950795"}],"value":"([2.043442272659671] -26.292886175950795)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.287047734809185</span>","value":"-26.287047734809185"}],"value":"([2.043442272659671] -26.287047734809185)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.31591192840172</span>","value":"-26.31591192840172"}],"value":"([2.043442272659671] -26.31591192840172)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.340092699623156</span>","value":"-26.340092699623156"}],"value":"([2.043442272659671] -26.340092699623156)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.29383329113715</span>","value":"-26.29383329113715"}],"value":"([2.043442272659671] -26.29383329113715)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.299724433695104</span>","value":"-26.299724433695104"}],"value":"([2.043442272659671] -26.299724433695104)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.46301448642859</span>","value":"-26.46301448642859"}],"value":"([2.043442272659671] -26.46301448642859)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.304765199318076</span>","value":"-26.304765199318076"}],"value":"([2.043442272659671] -26.304765199318076)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.31620451698342</span>","value":"-26.31620451698342"}],"value":"([2.043442272659671] -26.31620451698342)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.29651147062012</span>","value":"-26.29651147062012"}],"value":"([2.043442272659671] -26.29651147062012)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.285774668867358</span>","value":"-26.285774668867358"}],"value":"([2.043442272659671] -26.285774668867358)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.289518753348204</span>","value":"-26.289518753348204"}],"value":"([2.043442272659671] -26.289518753348204)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.321399060018408</span>","value":"-26.321399060018408"}],"value":"([2.043442272659671] -26.321399060018408)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.291719012171818</span>","value":"-26.291719012171818"}],"value":"([2.043442272659671] -26.291719012171818)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.291642014843802</span>","value":"-26.291642014843802"}],"value":"([2.043442272659671] -26.291642014843802)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.3870495943686</span>","value":"-26.3870495943686"}],"value":"([2.043442272659671] -26.3870495943686)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.840389579270344</span>","value":"-26.840389579270344"}],"value":"([2.043442272659671] -26.840389579270344)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.373344049215163</span>","value":"-26.373344049215163"}],"value":"([2.043442272659671] -26.373344049215163)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.545867889168637</span>","value":"-26.545867889168637"}],"value":"([2.043442272659671] -26.545867889168637)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.285533043921927</span>","value":"-26.285533043921927"}],"value":"([2.043442272659671] -26.285533043921927)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.29185864191977</span>","value":"-26.29185864191977"}],"value":"([2.043442272659671] -26.29185864191977)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.28634495850172</span>","value":"-26.28634495850172"}],"value":"([2.043442272659671] -26.28634495850172)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.333293939614407</span>","value":"-26.333293939614407"}],"value":"([2.043442272659671] -26.333293939614407)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.290124288248105</span>","value":"-26.290124288248105"}],"value":"([2.043442272659671] -26.290124288248105)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.351358810316565</span>","value":"-26.351358810316565"}],"value":"([2.043442272659671] -26.351358810316565)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.364149621608078</span>","value":"-26.364149621608078"}],"value":"([2.043442272659671] -26.364149621608078)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.29116825482903</span>","value":"-26.29116825482903"}],"value":"([2.043442272659671] -26.29116825482903)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.38294371653386</span>","value":"-26.38294371653386"}],"value":"([2.043442272659671] -26.38294371653386)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.548672603788525</span>","value":"-26.548672603788525"}],"value":"([2.043442272659671] -26.548672603788525)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.31471502276748</span>","value":"-26.31471502276748"}],"value":"([2.043442272659671] -26.31471502276748)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.294373520297267</span>","value":"-26.294373520297267"}],"value":"([2.043442272659671] -26.294373520297267)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.301562571569303</span>","value":"-26.301562571569303"}],"value":"([2.043442272659671] -26.301562571569303)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.3233540892188</span>","value":"-26.3233540892188"}],"value":"([2.043442272659671] -26.3233540892188)"},{"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>2.043442272659671</span>","value":"2.043442272659671"}],"value":"[2.043442272659671]"},{"type":"html","content":"<span class='clj-double'>-26.321426937865628</span>","value":"-26.321426937865628"}],"value":"([2.043442272659671] -26.321426937865628)"}],"value":"[([-1.4248822047652165] -30.083381179578573) ([-1.4248822047652165] -30.12782126743909) ([-1.543133703807011] -29.241429564881507) ([-1.543133703807011] -29.1420594612676) ([-1.543133703807011] -29.56356228288295) ([-1.543133703807011] -29.678762847408148) ([-1.543133703807011] -29.122961630797775) ([-1.6717662182040178] -28.480068501778426) ([2.043442272659671] -26.28769463307808) ([2.043442272659671] -26.340927493309696) ([2.043442272659671] -26.302408990460034) ([2.043442272659671] -26.29135508269127) ([2.043442272659671] -26.353294691878464) ([2.043442272659671] -26.2934672021875) ([2.043442272659671] -26.332084339700987) ([2.043442272659671] -26.572398822979498) ([2.043442272659671] -26.292886175950795) ([2.043442272659671] -26.287047734809185) ([2.043442272659671] -26.31591192840172) ([2.043442272659671] -26.340092699623156) ([2.043442272659671] -26.29383329113715) ([2.043442272659671] -26.299724433695104) ([2.043442272659671] -26.46301448642859) ([2.043442272659671] -26.304765199318076) ([2.043442272659671] -26.31620451698342) ([2.043442272659671] -26.29651147062012) ([2.043442272659671] -26.285774668867358) ([2.043442272659671] -26.289518753348204) ([2.043442272659671] -26.321399060018408) ([2.043442272659671] -26.291719012171818) ([2.043442272659671] -26.291642014843802) ([2.043442272659671] -26.3870495943686) ([2.043442272659671] -26.840389579270344) ([2.043442272659671] -26.373344049215163) ([2.043442272659671] -26.545867889168637) ([2.043442272659671] -26.285533043921927) ([2.043442272659671] -26.29185864191977) ([2.043442272659671] -26.28634495850172) ([2.043442272659671] -26.333293939614407) ([2.043442272659671] -26.290124288248105) ([2.043442272659671] -26.351358810316565) ([2.043442272659671] -26.364149621608078) ([2.043442272659671] -26.29116825482903) ([2.043442272659671] -26.38294371653386) ([2.043442272659671] -26.548672603788525) ([2.043442272659671] -26.31471502276748) ([2.043442272659671] -26.294373520297267) ([2.043442272659671] -26.301562571569303) ([2.043442272659671] -26.3233540892188) ([2.043442272659671] -26.321426937865628)]"}
;; <=

;; @@

;; @@
