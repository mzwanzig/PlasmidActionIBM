; Plasmid Action Model
;
; An individual-based model considering plasmid incompatibility and
; the regulation of bacterial conjugation by the switch of transfer competence
;
; Copyright 2017 Martin Zwanzig (neé Werisch)
; See Info tab for reference

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; EXTENSIONS ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

extensions [vid]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; PARAMETERS ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals [ ; Parameters that have an equal effect on all or a subgroup of individuals
;; Parameters defining characteristics of bacteria-plasmid-asscociations or the simulated environment ;;;;;
  ;here annotated with ';' as they are already determined in the general user interface (GUI):
  ;alpha_N                    ;burden to bear a non-transmissible plasmid (N)
  ;alpha_R                    ;burden to bear a repressed transmissible plasmid (R)
  ;alpha_D                    ;burden to bear a derepressed transmissible plasmid (D)
  ;kappa_N                    ;Probability for N to become excluded by segregation due to the co-occurrence of R or D in the same mother cell
  ;kappa_RD                   ;Probability for R and D to become excluded by segregation due to the co-occurrence of N in the same mother cell
  ;rho_R                      ;parameter of a function that determines the probability to become repressed
  ;rho_D                      ;parameter of a function that determines the probability to become derepressed
  ;gamma                      ;probability to perform a successful transfer attempt to a neighboring recipient
  ;omega                      ;probability for a bacterium and its plasmid(s) to become washed out
;; Setup parameters ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;simulation-period-in-generations  ; Duration of the simulation in maximum number of generations (constant)
  ;initial-N-proportion       ;proportion of bacteria bearing non-transmissible plasmids on all plasmid-bearers
  ;define-relative-plasmid-burden?   ;the plasmid-burdens alpha_N, alpha_R and alpha_D can be defined in relation to each other
                              ;if so, then alpha_D is the maximal plasmid burden and a reference for alpha_R, which is a reference for alpha_N
  ;alpha_N-alpha_R-ratio      ;alpha_N is defined in relation to alpha_R
  ;alpha_R-alpha_D-ratio      ;alpha_R is defined in relation to alpha_D
  ;MODEL ROBUSTNESS TESTS
  ;initial-demixing?          ;initialize the model in demixed state (clustered) instead of random positioning of F, R, and N?
  ;permanently-derepressed?   ;simulate permanently derepression of bacteria bearing transmissible plasmids?
                              ;every transmissible plasmid is and remains a derepressed one, starting from initialization
  ;incompatibility-exclusion? ;account for incompatibility exclusion? (= missegregation due to co-occurrence of either N and R or N and D)
  ;segregation?               ;account for missegregation of single plasmid types during bacterial fission
  ;global-dispersal?          ;simulate global dispersal of bacteria or not?
;; Procedure helper variables (store (local) outcomes of procedure substeps) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;this variables are used in the context of a single cell (and recycled for the next cell)
  target.                     ;for recolonization procedure
  success.                    ;for recolonization procedure
  receiver.                   ;for recolonization (incompatibility) procedure
  ;Defines the subgroup of plasmid-free cells (F):
  F                           ;colorcode of plasmid-free cells
  ;Defines the subgroup of cells bearing any plasmid (N, R or D):
  P                           ;colorcode of plasmid-bearing cells
  vid-counter                 ;used to make a frame for the video only every x ticks
;; Reporter variables ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Fc                          ;counts plasmid-free cells F
  rFc                         ;as a measure for the relative frequency of F to
  Rc                          ;counts cells bearing repressed transmissible plasmids R
  Dc                          ;counts cells bearing derepressed transmissible plasmids D
  Nc                          ;counts cells bearing non-transmissible plasmids N
  NTc                         ;counts cells inhibiting both N and R or D plasmids
  K                           ;a constant for the maximal attainable cell number (= patch count; constant)
  growth-rate                 ;a list of the growth rates per timestep (used to calculate generations)
  generations                 ;counts simulated generations
  fission                     ;counts successful cell division events per timestep
  transmission                ;counts successful conjugation events per timestep
  repression                  ;counts repression events per timestep
  derepression                ;counts derepression events per timestep
  R-exclusion                 ;counts R plasmids excluded for incompatibility reasons per timestep
  D-exclusion                 ;counts D plasmids excluded for incompatibility reasons per timestep
  N-exclusion                 ;counts N plasmids excluded for incompatibility reasons per timestep
]

breed [N Ns]
breed [R Rs]
breed [D Ds]

turtles-own [
  pb. ;plasmid-burden
]

to reset ;to set parameter values to default, press in interface before pressing "setup"
  set record-video? FALSE
  ; initial conditions and robustness tests
  set initial-N-proportion 0.5
  set initial-spatial-distribution "randomly-mixed"
  set permanently-derepressed? FALSE
  set incompatibility-exclusion? TRUE
  set consider-segregation? TRUE
  set suppressed_HGT_of_NT? FALSE
  set global-dispersal? FALSE
  ; stop conditions
  set simulation-period-in-generations 10000
  set stop_when_N=0_? TRUE
  set stop_when_R+D=0_? TRUE
  ; process rates
  set tau -4
  set alpha_N 0.0625
  set alpha_R 0.125
  set alpha_D 0.5
  set kappa_N 0
  set kappa_RD 0.5
  set rho_R -1
  set rho_D -3
  set gamma 0.3
  set omega 0.2
  set define-relative-plasmid-burden? FALSE
  set alpha_R-alpha_D-ratio 0.25
  set alpha_N-alpha_R-ratio 0.5
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; SETUP PROCEDURES ;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all

  ; define some global variables
  set F 3 ;define color of plasmid-free bacteria (F)
  set P 5 ;define color of plasmid-bearing bacteria (N, R or D)
  set success. FALSE ;helper variable used for recolonization procedure
  set K world-width * world-height ;defines the number of patches
  set growth-rate (list 0) ;creates a list, where the single growth rates of each timestep are stored
  if define-relative-plasmid-burden? ;alpha_R and alpha_N are either defined directly or in relation to alpha_D and each other
  [
    set alpha_R alpha_D * alpha_R-alpha_D-ratio ;alpha_R is a certain fraction of alpha_D
    set alpha_N alpha_R * alpha_N-alpha_R-ratio ;alpha_N is a certain fraction of alpha_R
  ]
  ;spread bacteria and plasmids
  ask patches [ ;any patch can become a bacterium by chance
    ifelse random-float 1 < 1 - omega  ;the probability depends on the setting for washout & mortality (this determines the total initial bacterial density)
    [
      set pcolor F ;become a plasmid-free bacterium by chance
      if initial-spatial-distribution = "randomly-mixed"
      [
        if random 2 = 1 ; become a plasmid-bearer by chance
        [
          become-a-plasmid-bearer ;procedure (see below)
        ]
      ]
    ]
    [ set pcolor white ] ;remain empty
  ]

  ;by default, the initial spatial distribution of bacteria and plasmid-bearers is random
  ;it can be changed to create a separation of single types
  if initial-spatial-distribution = "portioned" [ ;
    ifelse initial-N-proportion = 0
    [;in a system comprising F and R/D, but no N, a central circle-shaped cluster is created
      ask patch (max-pxcor / 2) (max-pxcor / 2) [
        become-a-plasmid-bearer
        ask patches in-radius (max-pxcor / 4) [
          if pcolor != white [
            set pcolor P
            ifelse permanently-derepressed?
            [sprout-D 1 [define-D] ]
            [sprout-R 1 [define-R] ]
          ]
        ]
      ]
    ]
    [;in a system comprising F, R/D and N three stripe-shaped clusters are created
      ask patches with [pxcor < (max-pxcor / 3)]
      [ if pcolor != white
        [
          set pcolor P
          ifelse permanently-derepressed?
          [sprout-D 1 [define-D] ]
          [sprout-R 1 [define-R] ]
        ]
      ]
      ask patches with [pxcor > (max-pxcor / 3) and pxcor < 2 * (max-pxcor / 3)]
      [ if pcolor != white
        [
          set pcolor P
          sprout-N 1 [define-N]
        ]
      ]
    ]
  ]
  reset-ticks

  if record-video? [ vid:start-recorder vid:record-interface ]
  make-video

end

to become-a-plasmid-bearer
  set pcolor P ;color of all plasmid-bearers
  ifelse random-float 1 <= initial-N-proportion ;get a non-transmissible plasmid by chance
  [ sprout-N 1 [define-N] ] ;get a non-transmissible plasmid 'N'
  [ ifelse permanently-derepressed? ;otherwise, get a transmissible plasmid 'R' or 'D'
    [sprout-D 1 [define-D] ] ;if transmissible plasmids are permanently-derepressed, get a 'D'
    [sprout-R 1 [define-R] ] ;if not, transmissible plasmids are initially in repressed state 'R'
  ]
end

to define-N ;characteristics of a non-transmissible plasmid N
  set shape "line" set heading 0 set color blue set pb. alpha_N
end

to define-R ;characteristics of a repressed transmissible plasmid R
  set shape "line" set heading 45 set color orange set pb. alpha_R
end

to define-D ;characteristics of a derepressed transmissible plasmid D
  set shape "line" set heading 90 set color lime set pb. alpha_D
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; RUNTIME PROCEDURES ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  measure-and-do-plots ;count the current frequency of single cell types and update plots
  ;stop the simulation, if any of the following conditions is fullfilled:
  if (generations >= simulation-period-in-generations) [stop] ;maximum simulation time is reached
  if stop_when_N=0_? [if(Nc = 0)  [stop]] ;non-transmissible plasmids N died out
  if stop_when_R+D=0_? [if (Rc + Dc = 0) [stop]] ;transmissible plasmids R and D died out
  ask patches [ ;asynchronous updating: focal points are selected randomly and updated one after the other, but each point only once in a single timestep
    ifelse pcolor != white
    [ ;each bacterium performs (by chance) single actions in this sequence:
      if permanently-derepressed? = FALSE [if any? R-here or any? D-here [sense-neighborhood-and-switch-transfer-competence]] ;a bacterium may switch from repressed to derepressed or vice versa
      if any? D-here [ifelse not suppressed_HGT_of_NT? [attempt-horizontal-gene-transfer] ;a transfer competent bacterium (carrying D) may perform horizontal gene transfer (default mode)
        [if not any? N-here [attempt-horizontal-gene-transfer]]] ;but HGT can be supressed by the co-occurrence of a non-transmissible plasmid in the same cell (a model robustness test)
      if random-float 1 <= omega [become-washed-out-or-die] ;a bacterium 'gets washed out or dies' by a probability omega
    ]
    [recolonization] ;an empty patch enables a neighboring bacterium to perform fission (by chance)
  ]
  if global-dispersal? [ ask patches [ ;global dispersal occurrs once in a timestep after performing the real asynchronous updating of the other processes
      get-dispersed ;bacteria move to a random empty patch and leave their original patch empty
    ]
  ]
  tick
  make-video
end

to sense-neighborhood-and-switch-transfer-competence ;SWITCH TRANSFER COMPETENCE
  ;is processed, if transmissible plasmids are not permanently derepressed and the focal bacterium is either repressed (carries 'R') or derepressed (carries 'D')
  let cB count neighbors with [pcolor != white] ;count all local bacteria
  if random (1 + cB) > 0 [;the more bacteria are present in the neighborhood, the stronger the signal and its sensing (high washout decreases the signal)
    let cRD count R-on neighbors + count D-on neighbors ;count all local bacteria bearing a transmissible plasmid R or D
    let FNprop ((cB - cRD) / cB) ;calculate the relative density of bacteria not bearing a transmissible plasmid R or D
    ifelse any? R-here ; either R or D are present in the focal cell, since this was the prerequisite checked in the 'go' procedure
    [
      if random-float 1 <= (10 ^ (rho_D * (1 - FNprop))) [ ;the probability to become derepressed increases with the local relative frequency of F + N
        ask R-here [ die ] sprout-D 1 [ define-D ] ;R is replaced by D
        set derepression derepression + 1 ;derepression event is counted
      ]
    ]
    [
      if random-float 1 <= (10 ^ (rho_R * FNprop)) [ ;the probability to become repressed decreases with the local relative frequency of F + N
        ask D-here [ die ] sprout-R 1 [ define-R ] ;D is replaced by R
        set repression repression + 1 ;repression event is counted
      ]
    ]
  ]
  if cB = 0 [ ; when no signals from neighboring cells sustain the derepression ...
    if any? D-here [ ; ...repression occurs, and ...
      ask D-here [ die ] sprout-R 1 [ define-R ] ;...D is replaced by R
      set repression repression + 1 ;repression event is counted
    ]
  ]
end

to attempt-horizontal-gene-transfer ;HORIZONTAL GENE TRANSFER
  ;is processed, if the focal bacterium carries a derepressed plasmid D
  let empty_neighbors count neighbors with [pcolor = white] ;check local availability of resources in terms of free space
  if empty_neighbors < 8 [ ;is there any bacterium to interact with?
    repeat random (1 + empty_neighbors) [ ;the more resources are available, the more transfer attempts can be performed (by chance)
      if random-float 1 <= gamma [ ;probability to perform a successful transfer attempt
        ask one-of neighbors with [pcolor != white] [
          if not any? D-here and not any? R-here [ ;check if the candidate is a potential recipient (not yet carrying R or D)
            sprout-D 1 [ define-D ] ;the cell receives D
            set pcolor P ;becomes a plasmid-bearer
            set transmission transmission + 1 ;the transmission-event is counted
          ]
        ]
      ]
    ]
  ]
end

to become-washed-out-or-die ;WASHOUT & MORTALITY
  ;is processed, if a random number lower than omega was drawn
  set pcolor white ask turtles-here [die]
end

to get-dispersed ;GLOBAL DISPERSAL
  if pcolor != white [ ;only bacteria become dispersed
    set target. one-of patches with [pcolor = white] ;select a random empty patch to move to
    ask target. [ set pcolor [pcolor] of myself ] ;change its color to the own ('F' or 'P')
    ask turtles-here [move-to target.] ;if the focal patch bears plasmids, they also move to the selected patch
    set pcolor white ;the focal patch becomes an empty one, indicated by white color
  ]
end

to recolonization ;VERTICAL GENE TRANSFER
  ;is processed, if the focal patch is empty
  if (any? other neighbors with [pcolor != white]) [ ;is there any bacterium in the local neighborhood that could perform fission?
    set target. one-of neighbors with [pcolor != white] ;select one at random
    let recipient. self ;helper variable to process fission and plasmid replication
    if random-float 1 <= (1 - sum [pb.] of turtles-on target.) [ ;the probability for fission is reduced by the plasmid-burden(s)
      ask turtles-on target. [ hatch 1 [ move-to recipient. ] ] ;perform plasmid replication
      set pcolor [pcolor] of target. ;perform fission
      set fission fission + 1 ;count this fission event
      ifelse incompatibility-exclusion?
      [ check-incompatibility ] ;consider missegregation due to incompatibility (in this computercode plasmids have been copied yet, but they might get erased again)
      [ check-segregation ]
    ]
  ]
  set success. FALSE ;this global variable is reset before the transfer attempt of the next bacterium is processed
end

to check-incompatibility ;INCOMPATIBILITY-EXCLUSION
  let Rh count R-here
  let Dh count D-here
  let Nh count N-here
  ifelse (Nh + Rh + Dh = 2) ;Does the bacterium carry both plasmid types?
  [
    ifelse random-float 1 <= kappa_RD ;by chance the transmissible plasmid (R or D) gets lost in one of the daughters
      [ ;there are either R or D present, which get lost in a random daughter
        ifelse any? R-here
        [ ;R instead of D are present; one daughter cell receives no R
          set R-exclusion R-exclusion + 1 ;count this event
          ifelse random 2 = 0 ;random selection of daughter
            [ ask R-here [die] set receiver. target.]
            [ ask R-on target. [die] set receiver. self]
        ]
        [ ;D instead of R are present; one daughter cell receives no D
          set D-exclusion D-exclusion + 1 ;count this event
          ifelse random 2 = 0 ;random selection of daughter
            [ ask D-here [die] set receiver. target.]
            [ ask D-on target. [die] set receiver. self]
        ]
        ;when R or D became excluded in one daughter, N might only become excluded in the other cell
        if random-float 1 <= kappa_N [;the cell that receives R or D probably receives no N
          set N-exclusion N-exclusion + 1 ;count this event
          ask N-on receiver. [die]
        ]
      ]
      [ ;even if R or D were not excluded, N-exclusion might happen to one of both daughters
        if random-float 1 <= kappa_N [;one daughter cell probably receives no N
          set N-exclusion N-exclusion + 1 ;count this event
          ifelse random 2 = 0 ;random selection of daughter
            [ ask N-here [die] ]
            [ ask N-on target. [die]]
        ]
      ]
  ]
  [
    check-segregation
  ]
end

to check-segregation ;either called after a check of incompatibility exclusion or directly, when incompatibility exclusion is not considered
  if consider-segregation? [
    let segregation-probability 10 ^ tau ;calculate the probability for segregation using the given estimate for tau
    if random-float 1 <= segregation-probability [ ;by chance one daughter cell will not receive a plasmid copy from the mother cell
      ifelse incompatibility-exclusion?
      [ ;if incompatibility exclusion is considered, only one plasmid can be present in the focal cell (see procedure 'check-incompatibility')
        ifelse random 2 = 0 ;random selection of daughter that looses its plasmid and becomes a plasmid-free cell F
        [ ask turtles-here [die] set pcolor F ]
        [ ask target. [ ask turtles-here [die] set pcolor F ]]
      ]
      [ ;if incompatibility exclusion is not considered, two plasmid types can be present in the focal cell
        if any? R-here [ifelse random 2 = 0 [ask R-here [die]] [ask target. [ask R-here [die]]]] ;one of the daughters looses its R plasmid
        if any? D-here [ifelse random 2 = 0 [ask D-here [die]] [ask target. [ask D-here [die]]]] ;one of the daughters looses its D plasmid
        if any? N-here [ifelse random 2 = 0 [ask N-here [die]] [ask target. [ask N-here [die]]]] ;one of the daughters looses its N plasmid
        if not any? turtles-here [ set pcolor F ] ;if no plasmid remains, this daughter becomes a plasmid-free cell F
        ask target. [ if not any? turtles-here [ set pcolor F ] ] ;if no plasmid remains, this daughter becomes a plasmid-free cell F
      ]
    ]
  ]
end

to make-video
  if record-video? [
    if vid-counter = 0
    [ vid:record-interface
      let num-generations 5 ; defines the number of generation after which a frame for the video is made; can be changed as required
      set vid-counter (1 / omega) * num-generations ] ; set vid-counter to 0 in order to make a frame for the video every tick
    set vid-counter vid-counter - 1
  ]
end

to save-recordings
  if record-video? [ vid:save-recording "out.mp4" ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; REPORTER AND PLOTTING PROCEDURES ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to measure-and-do-plots
  set Fc count patches with [pcolor = F]
  set Rc count R
  set Dc count D
  set Nc count N
  set NTc ((count patches with [pcolor = P]) - Rc - Dc - Nc) / -1 ;some cells harbor both N and R or D
  if ticks > 1 [set generations (ticks * (mean (growth-rate)))]
  set growth-rate lput (fission / (Fc + Rc + Dc + Nc - NTc)) growth-rate
  do-plots ;procedure to update plots in GUI (see below)
  ;reset the tick-based reporter variables:
  set fission 0
  set repression 0
  set derepression 0
  set transmission 0
  set D-exclusion 0
  set R-exclusion 0
  set N-exclusion 0
end

to do-plots

  set-current-plot "Dynamics"
  set-current-plot-pen "Plasmid-free"
    plot Fc / K
  set-current-plot-pen "Repressed"
    plot Rc / K
  set-current-plot-pen "Derepressed"
    plot Dc / K
  set-current-plot-pen "Non-transmissible"
    plot Nc / K
  set-current-plot-pen "NT"
    plot NTc / K

  set-current-plot "Transition-of-transfer-competence"
  set-current-plot-pen "R-derepression"
    ifelse (Rc > 0) [plot derepression / Rc] [plot 0]
  set-current-plot-pen "D-repression"
    ifelse (Dc > 0) [plot repression / Dc] [plot 0]

  set-current-plot "Plasmid-gain-loss"
  set-current-plot-pen "transmission"
    plot transmission
  set-current-plot-pen "R-exclusion"
    plot R-exclusion
  set-current-plot-pen "D-exclusion"
    plot D-exclusion
  set-current-plot-pen "N-exclusion"
    plot N-exclusion

  set-current-plot "F-T-plane"
  set-current-plot-pen "F-T"
    plotxy (Fc / K) ((Rc + Dc) / K)

  set-current-plot "F-N-plane"
  set-current-plot-pen "F-N"
    plotxy (Fc / K) (Nc / K)

  set-current-plot "T-N-plane"
  set-current-plot-pen "T-N"
    plotxy ((Rc + Dc) / K) (Nc / K)

end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
1018
819
-1
-1
4.0
1
10
1
1
1
0
1
1
1
0
199
0
199
1
1
1
ticks
8.0

BUTTON
155
43
210
76
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
156
76
211
109
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
1377
798
1549
831
omega
omega
0
1
0.2
0.01
1
NIL
HORIZONTAL

TEXTBOX
1378
597
1739
692
HORIZONTAL GENE TRANSFER\n============================================\nA bacterium that harbours a derepressed plasmid D may perform a series of transfer attempts. The number of attempts is influenced by the local resource availability (=empty patches). Each attempt comprises the selection of a random bacterium from the local neighborhood. If this is a recipient, it receives the plasmid D with the probability gamma.
10
0.0
1

SWITCH
17
706
184
739
global-dispersal?
global-dispersal?
1
1
-1000

SLIDER
1377
699
1549
732
gamma
gamma
0
1
0.3
0.01
1
NIL
HORIZONTAL

SLIDER
1376
208
1548
241
alpha_N
alpha_N
0
1
0.0625
0.01
1
NIL
HORIZONTAL

SLIDER
1376
170
1548
203
alpha_R
alpha_R
0
1
0.125
0.01
1
NIL
HORIZONTAL

SLIDER
1376
132
1548
165
alpha_D
alpha_D
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
1376
548
1548
581
rho_D
rho_D
-10
0
-3.0
1
1
NIL
HORIZONTAL

PLOT
1028
120
1367
264
Dynamics
time
proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Plasmid-free" 1.0 0 -11053225 true "" ""
"Repressed" 1.0 0 -955883 true "" ""
"Derepressed" 1.0 0 -13840069 true "" ""
"Non-transmissible" 1.0 0 -13345367 true "" ""
"NT" 1.0 0 -2064490 true "" ""

SLIDER
14
207
195
240
initial-N-proportion
initial-N-proportion
0
1
0.5
0.01
1
NIL
HORIZONTAL

SWITCH
14
474
183
507
incompatibility-exclusion?
incompatibility-exclusion?
0
1
-1000

SLIDER
1377
349
1549
382
kappa_RD
kappa_RD
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
1588
348
1760
381
kappa_N
kappa_N
0
1
0.0
0.01
1
NIL
HORIZONTAL

TEXTBOX
1377
317
1584
343
Probability that R or D is excluded, when N is present in the same cell:
10
0.0
1

TEXTBOX
1590
317
1779
343
Probability that N is excluded, when R or D is present in the same cell:
10
0.0
1

TEXTBOX
1380
16
1760
111
VERTICAL GENE TRANSFER\n============================================\nBurden of non-transmissible (N), repressed (R) and derepressed (D) plasmid types; Burden means that the probability to perform bacterial fission when an empty lattice cell can be occupied is reduced by x, resulting in a probability of 1 - x. You can either define the plasmid burden alpha_N and alpha_R directly (see left) or in relation to each other (with reference to alpha_D).
10
0.0
1

PLOT
1031
706
1353
826
Transition-of-transfer-competence
time
proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"R-derepression" 1.0 0 -13840069 true "" ""
"D-repression" 1.0 0 -955883 true "" ""

PLOT
1030
583
1350
703
Plasmid-gain-loss
time
frequency
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"transmission" 1.0 0 -7500403 true "" ""
"R-exclusion" 1.0 0 -955883 true "" ""
"D-exclusion" 1.0 0 -10899396 true "" ""
"N-exclusion" 1.0 0 -13345367 true "" ""

PLOT
1028
271
1188
426
F-T-plane
F
T
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"F-T" 1.0 0 -16777216 true "" ""

PLOT
1029
427
1189
579
F-N-plane
F
N
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"F-N" 1.0 0 -16777216 true "" ""

PLOT
1194
427
1354
579
T-N-plane
T
N
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"T-N" 1.0 0 -16777216 true "" ""

SWITCH
12
353
199
386
permanently-derepressed?
permanently-derepressed?
1
1
-1000

SWITCH
1561
117
1762
150
define-relative-plasmid-burden?
define-relative-plasmid-burden?
1
1
-1000

SLIDER
1561
190
1762
223
alpha_N-alpha_R-ratio
alpha_N-alpha_R-ratio
0
1.1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
1561
153
1763
186
alpha_R-alpha_D-ratio
alpha_R-alpha_D-ratio
0
1
0.25
0.01
1
NIL
HORIZONTAL

SLIDER
1558
548
1730
581
rho_R
rho_R
-10
0
-1.0
1
1
NIL
HORIZONTAL

TEXTBOX
1377
258
1742
308
INCOMPATIBILITY-EXCLUSION\n============================================\nBacteria carrying both plasmid types (either N and R or N and D) may loose one plasmid type during fission in one of the resulting daughter cells.
10
0.0
1

TEXTBOX
14
110
189
138
 MODEL ROBUSTNESS TESTS\n==================
11
0.0
1

TEXTBOX
1378
401
1734
524
SWITCH TRANSFER COMPETENCE\n============================================\nBacterial conjugation can be regulated by a switch of an individual cells transfer competence. This switch can be affected by signaling molecules as pheromones, that induce a repression or derepression of transfer genes. Define the sensitivity for derepression rho_D and repression rho_R to occur in dependence to the local proportion of bacteria that do not carry a transmissible plasmid (R or D). The lower the value, the less likely is a switch at low proportions (=signaling molecules).
10
0.0
1

TEXTBOX
1379
748
1746
791
WASHOUT & MORTALITY\n============================================\nA bacterium dies or is washed out with the probability omega.
10
0.0
1

INPUTBOX
1222
52
1359
112
simulation-period-in-generations
10000.0
1
0
Number

CHOOSER
12
286
181
331
initial-spatial-distribution
initial-spatial-distribution
"randomly-mixed" "portioned"
0

TEXTBOX
15
139
210
206
Initially, a fraction of 0.5 of all bacteria are plasmid-bearers, carrying repressed transmissible plasmids R or non- transmissible plasmids N. Define the fraction of plasmid-bearers carrying N:
10
0.0
1

TEXTBOX
14
255
203
283
Define how plasmids are spatially distributed in the initial model world:
10
0.0
1

SWITCH
17
572
184
605
consider-segregation?
consider-segregation?
0
1
-1000

SLIDER
17
656
119
689
tau
tau
-7
-1
-4.0
1
1
NIL
HORIZONTAL

TEXTBOX
1377
531
1527
549
Sensitivity for derepression
10
0.0
1

TEXTBOX
1558
531
1706
549
Sensitivity for repression
10
0.0
1

TEXTBOX
15
389
194
473
If ON, transmissible plasmids are always transfer competent. Thus, bacterial conjugation is not regulated by switching between the repressed and derepressed mode.\n-------------------------------------
10
0.0
1

TEXTBOX
18
515
198
567
If OFF, the co-occurrence of both plasmid types in the same cell does not affect their vertical transmission.\n-------------------------------------
10
0.0
1

TEXTBOX
19
608
206
651
If ON, a cell may looses its plasmid during fission in one of both daughters with a probability of 10 ^ tau.
10
0.0
1

TEXTBOX
18
691
168
709
-------------------------------------
11
0.0
1

TEXTBOX
18
743
200
884
If ON, any bacterium moves to a random non-occupied patch in the model world (after performing some other actions by chance), leaving the current patch non-occupied. As a consequence, bacteria still interact locally, but mixing continuously changes their local neighborhood.
10
0.0
1

MONITOR
1194
273
1328
318
executed generations
generations
0
1
11

TEXTBOX
15
335
165
353
-------------------------------------
11
0.0
1

TEXTBOX
1198
17
1366
46
Define the number of generations that should be executed:
10
0.0
1

SWITCH
1025
41
1186
74
stop_when_N=0_?
stop_when_N=0_?
0
1
-1000

SWITCH
1026
79
1203
112
stop_when_R+D=0_?
stop_when_R+D=0_?
0
1
-1000

MONITOR
1195
322
1318
367
current growth rate
last growth-rate
2
1
11

TEXTBOX
1028
10
1152
38
 STOP CONDITIONS\n=============
11
0.0
1

SWITCH
1555
699
1760
732
suppressed_HGT_of_NT?
suppressed_HGT_of_NT?
1
1
-1000

BUTTON
155
10
210
43
NIL
reset
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
15
42
149
75
record-video?
record-video?
1
1
-1000

BUTTON
13
77
150
110
NIL
save-recordings
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
15
10
140
41
Video? 1. \"setup\", 2. \"go\", 3. \"save\"
12
0.0
1

@#$#@#$#@
## WHAT IS IT?

This is a model about plasmid dynamics. It primarily explores the faith of plasmids with and without an interaction between non-transmissible and transmissible plasmids in an isogenic bacterial population. Furthermore, it demonstrates how local conditions triggering switches in the individual transfer competence of bacteria bearing transmissible plasmids affect population dynamics.

## HOW IT WORKS

Hosts with non-transmissible plasmids are able to outcompete hosts with transmissible plasmids (as the transmissible plasmids are more costly and more easily lost from a host with both plasmids). Hosts with transmissible plasmids can lead to decreases in the density of plasmid-free cells (through conjugation). Finally, plasmid-free cells outcompete hosts with non-transmissible plasmids (due to the cost of these plasmids). The upshot of this non-transitive dynamic is the maintenance of the non-transmissible plasmid.

## HOW TO USE IT

To start a default simulation click on "setup" and then "go" in the interface. Parameter can be changed live, that means while the simulation is running. In order to run a more distinct simulation, stop the model run by disabling "go", change the desired parameters and start a new simulation by clicking first on "setup" and then on "go". Please note that the plasmid-burden is either defined directly or as a relative measure to another plasmid burden. To disable the relative definition you have to switch of 'define-relative-plasmid-burden?'. If your model stops when a plasmid died out and you wish to run it for an extended period, you have to disable 'stop_when_N=0_?' and 'stop_when_R+D=0_?', respectively.

## THINGS TO NOTICE

The model performs a sort of warm up. That means although all bacteria may be initially mixed, they will separate into clusters in a period of some generations. During this time the model behaves different than in the time that follows after this warm up. Certain parameter settings may cause a collapse of the population (extinction of plasmid-bearing or plasmid-free bacteria) as the oscillations can be extraordinary pronounced in this warm up phase.

## THINGS TO TRY

Change sliders and switches at the left side of the interface (referring to model robustness tests). Start to change one at a time, before you undertake parallel adjustments. If the 'initial-N-proportion' is 0, you will only consider the presence of transmissible plasmids, while a value of 1 considers only non-transmissible plasmids. Both cause a distinct change in model behavior. But independent of a small or huge initial proportion of the one or the other plasmid type, the model stabilizes quite likely at the same level, i.e. performs comparable oscillations after warm up, for values >0 and <1.

## EXTENDING THE MODEL

A feasible model extension comprises the implementation of another sensing mechanism, that determines a bacteriums switch of transfer competence. It might also be feasible to investigate the outcome of an evolution of certain parameters, e.g. the transfer probability that affects conjugation. Apart from this, it can also be interesting to include a direct selection effect, e.g. through the presence of heavy metals or antibiotics. Given a certain additional cost for resistance, would transmissible or non-transmissible plasmids be doomed as sensitive plasmid-free bacteria?

## NETLOGO FEATURES

This model was implemented in NetLogo 6.0 and refurbished in NetLogo 6.2.2 including a 'reset'-procedure to set the model parameters to default values and including a video option using  no additional features.

## RELATED MODELS

By Martin Zwanzig (neé Werisch), more plasmid dynamics models have been developed an published. You can find them on GitHub: https://github.com/mzwanzig

A popular simulation model demonstrating the rock-paper-scissor mechanism for the instransitive dynamics in another system is published along with this article:
Kerr, B., Riley, M.A., Feldman, M.W., Bohannan, B.J.M., 2002. Local dispersal promotes biodiversity in a real-life game of rock-paper-scissors. Nature 418, 171–174. http://dx.doi.org/10.1038/nature00823

## CREDITS AND REFERENCES

Please indicate any use of this model that contributes to a publication with a reference to the following article:

Werisch, M., Berger, U. Berendonk, T. (2017): Conjugative plasmids enable the maintenance of low cost non-transmissible plasmids. Plasmid 91: 96-104.
https://doi.org/doi:10.1016/j.plasmid.2017.04.004
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="def-experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="world-size-experiment1" repetitions="100" runMetricsEveryStep="false">
    <setup>resize-world 0 49 0 49
setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="world-size-experiment2" repetitions="100" runMetricsEveryStep="false">
    <setup>resize-world 0 99 0 99
setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="world-size-experiment3" repetitions="100" runMetricsEveryStep="false">
    <setup>resize-world 0 149 0 149
setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="world-size-experiment4" repetitions="100" runMetricsEveryStep="false">
    <setup>resize-world 0 199 0 199
setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="world-size-experiment5" repetitions="100" runMetricsEveryStep="false">
    <setup>resize-world 0 249 0 249
setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="N-burden_experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <steppedValueSet variable="alpha_N" first="0" step="0.005" last="0.2"/>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gamma-omega_experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="gamma" first="0.1" step="0.1" last="1"/>
    <steppedValueSet variable="omega" first="0.1" step="0.1" last="0.5"/>
  </experiment>
  <experiment name="kappaRD-kappaN_experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <steppedValueSet variable="kappa_RD" first="0" step="0.1" last="1"/>
    <steppedValueSet variable="kappa_N" first="0" step="0.1" last="1"/>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="rhoD-rhoR_experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <steppedValueSet variable="rho_D" first="-9" step="1" last="-1"/>
    <steppedValueSet variable="rho_R" first="-9" step="1" last="-1"/>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="initial-N-proportion_experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.001"/>
      <value value="0.01"/>
      <value value="0.1"/>
      <value value="0.5"/>
      <value value="0.9"/>
      <value value="0.99"/>
      <value value="0.999"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="alphaR-alphaD_experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <steppedValueSet variable="alpha_D" first="0" step="0.01" last="0.5"/>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.125"/>
      <value value="0.25"/>
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="world-size-experiment6" repetitions="100" runMetricsEveryStep="false">
    <setup>resize-world 0 24 0 24
setup</setup>
    <go>go</go>
    <metric>generations</metric>
    <metric>Fc</metric>
    <metric>Rc</metric>
    <metric>Dc</metric>
    <metric>Nc</metric>
    <metric>NTc</metric>
    <metric>fission</metric>
    <metric>transmission</metric>
    <metric>repression</metric>
    <metric>derepression</metric>
    <metric>R-exclusion</metric>
    <metric>D-exclusion</metric>
    <metric>N-exclusion</metric>
    <enumeratedValueSet variable="initial-N-proportion">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-spatial-distribution">
      <value value="&quot;randomly-mixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="permanently-derepressed?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incompatibility-exclusion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-segregation?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tau">
      <value value="-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_N=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop_when_R+D=0_?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-period-in-generations">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N">
      <value value="0.0625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R">
      <value value="0.125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_D">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="define-relative-plasmid-burden?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_N-alpha_R-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_R-alpha_D-ratio">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_RD">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kappa_N">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_D">
      <value value="-3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rho_R">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gamma">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="omega">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
