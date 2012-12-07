(module stochastic-discrete *
(import chicken scheme extras)
(require-extension nondeterminism define-structure srfi-1 traversal)

(define  *gm-strategy* 'ac)

;; (define-macro fold-distribution
;;  (lambda (form expander)
;;   (unless (= (length form) 4)
;;    (error 'support "Improper FOLD-DISTRIBUTION: ~s" form))
;;   (expander `(fold-distribution-thunk
;; 	      ,(second form) ,(third form) (lambda () ,(fourth form)))
;; 	    expander)))

(define-syntax fold-distribution
 (syntax-rules ()
  ((_ f i thunk) (fold-distribution-thunk f i (lambda () thunk)))))

;; (define-macro support
;;  (lambda (form expander)
;;   (unless (= (length form) 2) (error 'support "Improper SUPPORT: ~s" form))
;;   (expander `(support-thunk (lambda () ,(second form))) expander)))

(define-syntax support
 (syntax-rules ()
  ((_ thunk) (support-thunk-thunk))))

(define-syntax supportq
 (syntax-rules ()
  ((_ thunk) (supportq-thunk-thunk))))

(define-syntax supportv
 (syntax-rules ()
  ((_ thunk) (supportv-thunk-thunk))))

(define-syntax supportp
 (syntax-rules ()
  ((_ thunk) (supportv-thunk-thunk))))

;; (define-macro probability
;;  (lambda (form expander)
;;   (unless (= (length form) 2)
;;    (error 'probability "Improper PROBABILITY: ~s" form))
;;   (expander `(probability-thunk (lambda () ,(second form))) expander)))

(define-syntax probability
 (syntax-rules ()
  ((_ thunk) (probability-thunk (lambda () thunk)))))

;; (define-macro expected-value
;;  (lambda (form expander)
;;   (unless (= (length form) 5)
;;    (error 'expected-value "Improper EXPECTED-VALUE: ~s" form))
;;   (expander `(expected-value-thunk
;; 	      ,(second form)
;; 	      ,(third form)
;; 	      ,(fourth form)
;; 	      (lambda () ,(fifth form)))
;; 	    expander)))

(define-syntax expected-value
 (syntax-rules ()
  ((_ plus times zero thunk) 
   (expected-value-thunk plus times zero (lambda () thunk)))))

;; (define-macro entropy
;;  (lambda (form expander)
;;   (unless (= (length form) 2) (error 'entropy "Improper ENTROPY: ~s" form))
;;   (expander `(entropy-thunk (lambda () ,(second form))) expander)))

(define-syntax entropy
 (syntax-rules ()
  ((_ thunk) (entropy-thunk (lambda () thunk)))))

;; (define-macro distribution
;;  (lambda (form expander)
;;   (unless (= (length form) 2)
;;    (error 'distribution "Improper DISTRIBUTION: ~s" form))
;;   (expander `(distribution-thunk (lambda () ,(second form))) expander)))

(define-syntax distribution
 (syntax-rules ()
  ((_ thunk) (distribution-thunk-thunk))))

(define-syntax distributionq
 (syntax-rules ()
  ((_ thunk) (distributionq-thunk-thunk))))

(define-syntax distributionv
 (syntax-rules ()
  ((_ thunk) (distributionv-thunk-thunk))))

(define-syntax distributionp
 (syntax-rules ()
  ((_ thunk) (distributionv-thunk-thunk))))

;;; TODO mutual information, K-L divergence, mean, median, mode, variance

(define flip (lambda (alpha) (error "Top-level flip")))

(define bottom (lambda () (error "Top-level bottom")))

(define current-probability (lambda () (error "Top-level current-probability")))

(define (fold-distribution-thunk f i thunk)
 (call-with-current-continuation
  (lambda (c)
   (let ((accumulation i) (p 1) (saved-flip flip)
	 (saved-bottom bottom) (saved-current-probability current-probability))
    (set! current-probability (lambda () p))
    (set! flip
	  (lambda (alpha)
	   (unless (<= 0 alpha 1) (error "Alpha not probability"))
	   (cond ((zero? alpha) #f)
		 ((= alpha 1) #t)
		 (else (call-with-current-continuation
			(lambda (c)
			 (let ((saved-p p) (saved-bottom bottom))
			  (set! p (* alpha p))
			  (set! bottom
				(lambda ()
				 (set! p (* (- 1 alpha) saved-p))
				 (set! bottom saved-bottom)
				 (c #f)))
			  #t)))))))
    (set! bottom
	  (lambda ()
	   (set! flip saved-flip)
	   (set! bottom saved-bottom)
	   (set! current-probability saved-current-probability)
	   (c accumulation)))
    (let ((value (thunk))) (set! accumulation (f value p accumulation)))
    (bottom)))))

(define (support-thunk thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (if (member value accumulation) accumulation (cons value accumulation)))
  '()
  thunk))

(define (supportq-thunk thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (if (memq value accumulation) accumulation (cons value accumulation)))
  '()
  thunk))

(define (supportv-thunk thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (if (memv value accumulation) accumulation (cons value accumulation)))
  '()
  thunk))

(define (supportp-thunk p? thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (if (memp p? value accumulation) accumulation (cons value accumulation)))
  '()
  thunk))

(define (probability-thunk thunk)
 (fold-distribution-thunk (lambda (value p accumulation)
			   (if (eq? value #t) (+ p accumulation) accumulation))
			  0
			  thunk))

(define (expected-value-thunk plus times zero thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation) (plus (times p value) accumulation))
  zero
  thunk))

(define (entropy-thunk thunk)
 (- 0
    (fold-distribution-thunk
     (lambda (value p accumulation) (+ (* p (log p)) accumulation))
     0
     thunk)))

(define (distribution-thunk thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (let ((pair (assoc value accumulation)))
    (if pair
	(cons (cons value (+ p (cdr pair))) (removeq pair accumulation))
	(cons (cons value p) accumulation))))
  '()
  thunk))

(define (distributionq-thunk thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (let ((pair (assq value accumulation)))
    (if pair
	(cons (cons value (+ p (cdr pair))) (removeq pair accumulation))
	(cons (cons value p) accumulation))))
  '()
  thunk))

(define (distributionv-thunk thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (let ((pair (assv value accumulation)))
    (if pair
	(cons (cons value (+ p (cdr pair))) (removeq pair accumulation))
	(cons (cons value p) accumulation))))
  '()
  thunk))

(define (distributionp-thunk p? thunk)
 (fold-distribution-thunk
  (lambda (value p accumulation)
   (let ((pair (assp p? value accumulation)))
    (if pair
	(cons (cons value (+ p (cdr pair))) (removeq pair accumulation))
	(cons (cons value p) accumulation))))
  '()
  thunk))

(define (most-likely-value distribution)
 (car (fold (lambda (pair best) (if (> (cdr pair) (cdr best)) pair best))
	    distribution
	    `(#f . -inf.0))))

(define (most-likely-pair distribution)
 (fold (lambda (pair best) (if (> (cdr pair) (cdr best)) pair best))
       distribution
       `(#f . -inf.0)))

(define (draw-pair distribution)
 (define (min x1 x2) (if (< x1 x2) x1 x2))
 (define (max x1 x2) (if (> x1 x2) x1 x2))
 (let loop ((distribution distribution) (p 1))
  (cond
   ((or (zero? p) (null? distribution)) (bottom))
   ((flip (min (/ (cdr (first distribution)) p) 1)) (first distribution))
   (else
    (loop (rest distribution) (max (- p (cdr (first distribution))) 0))))))

(define (draw distribution)
 (let loop ((p 1)
	    (pairs
	     ;; This is done not so much for efficiency but to prevent a
	     ;; situation where the last pair has zero probability and roundoff
	     ;; error is subtracting the earlier probabilities from one leads
	     ;; to flip being called with a value slightly greater than one.
	     (remove-if (lambda (pair) (zero? (cdr pair))) distribution)))
  (if (null? pairs)
      (bottom)
      (if (flip (/ (cdr (first pairs)) p))
	  (car (first pairs))
	  (loop (- p (cdr (first pairs))) (rest pairs))))))


;;; Graphical models

;;; domain -> distribution
;;; csp -> gm
;;; fail -> bottom

;;; This assume a model where there is no renormalization due to bottoms.

(define-structure distribution-variable distribution demons)

(define (upon-bottom-thunk thunk)
 (let ((saved-bottom bottom))
  (set! bottom (lambda () (thunk) (saved-bottom)))))

(define (create-distribution-variable distribution)
 (let ((distribution
	(remove-if (lambda (pair) (zero? (cdr pair))) distribution)))
  (when (null? distribution) (bottom))
  (make-distribution-variable distribution '())))

(define (restrict-distribution! x distribution)
 (define (local-set-distribution-variable-distribution!
	  distribution-variable distribution)
  (let ((distribution
	 (distribution-variable-distribution distribution-variable)))
   (upon-bottom-thunk
    (lambda ()
     (set-distribution-variable-distribution!
      distribution-variable distribution))))
  (set-distribution-variable-distribution! distribution-variable distribution))
 (when (null? distribution) (bottom))
 (when (< (length distribution)
	  (length (distribution-variable-distribution x)))
  (local-set-distribution-variable-distribution! x distribution)
  (for-each (lambda (demon) (demon)) (distribution-variable-demons x))))

(define (gm-bound? x) (null? (rest (distribution-variable-distribution x))))

(define (gm-binding x) (car (first (distribution-variable-distribution x))))

(define (gm-solution ds)
 (let loop ((ds ds) (xs '()))
  (if (null? ds)
      (reverse xs)
      (let ((pair
	     (draw-pair (distribution-variable-distribution (first ds)))))
       (restrict-distribution! (first ds) (list pair))
       (loop (rest ds) (cons (first pair) xs))))))

(define (gm-bb-solution ds)
 (let ((best 0))
  (let loop ((ds ds) (xs '()))
   (when (<= (current-probability) best) (bottom))
   (cond
    ((null? ds)
     (when (< best (current-probability)) (set! best (current-probability)))
     (reverse xs))
    (else
     (let ((pair
	    (draw-pair (distribution-variable-distribution (first ds)))))
      (restrict-distribution! (first ds) (list pair))
      (loop (rest ds) (cons (first pair) xs))))))))

(define (gm-bb-predict-solution ds)
 (let ((best 0))
  (let loop ((ds ds) (xs '()))
   (when (or (<= (current-probability) best)
	     (<= (foldl
		  (lambda (v d)
		   (* v
		      (map-reduce max
				  -inf.0
				  cdr
				  (distribution-variable-distribution d))))
		  ds
		  (current-probability))
	     	 best))
    (bottom))
   (cond
    ((null? ds)
     (when (< best (current-probability)) (set! best (current-probability)))
     (reverse xs))
    (else
     (let ((pair
	    (draw-pair (distribution-variable-distribution (first ds)))))
      (restrict-distribution! (first ds) (list pair))
      (loop (rest ds) (cons (first pair) xs))))))))

(define (gm-bb-predict-solution-with-start ds start)
 (let ((best start))
  (let loop ((ds ds) (xs '()))
   (when (or (<= (current-probability) best)
	     (<= (foldl
		  (lambda (v d)
		   (* v
		      (map-reduce max
				  -inf.0
				  cdr
				  (distribution-variable-distribution d))))
		  ds
		  (current-probability))
	     	 best))
    (bottom))
   (cond
    ((null? ds)
     (when (< best (current-probability)) (set! best (current-probability)))
     (reverse xs))
    (else
     (let ((pair
	    (draw-pair (distribution-variable-distribution (first ds)))))
      (restrict-distribution! (first ds) (list pair))
      (loop (rest ds) (cons (first pair) xs))))))))

(define (some-element p x)
 (some (lambda (x) (p (car x))) (distribution-variable-distribution x)))

(define (one-element p x)
 (one (lambda (x) (p (car x))) (distribution-variable-distribution x)))

(define (the-element p x)
 (list (find-if (lambda (x) (p (car x)))
		(distribution-variable-distribution x))))

(define (the-elements p x)
 (remove-if-not (lambda (x) (p (car x)))
		(distribution-variable-distribution x)))

(define (attach-demon! demon x)
 (define (local-set-distribution-variable-demons! distribution-variable demons)
  (let ((demons (distribution-variable-demons distribution-variable)))
   (upon-bottom-thunk
    (lambda ()
     (set-distribution-variable-demons! distribution-variable demons))))
  (set-distribution-variable-demons! distribution-variable demons))
 (local-set-distribution-variable-demons!
  x (cons demon (distribution-variable-demons x)))
 (demon))

(define (gm-assert-constraint-efd! constraint xs)
 (for-each
  (lambda (x)
   (attach-demon! (lambda ()
		   (when (every gm-bound? xs)
		    (unless (apply constraint (map gm-binding xs)) (bottom))))
		  x))
  xs))

(define (gm-assert-constraint-fc! constraint xs)
 (for-each
  (lambda (x)
   (attach-demon!
    (lambda ()
     (when (one (lambda (x) (not (gm-bound? x))) xs)
      (let* ((i (position-if (lambda (x) (not (gm-bound? x))) xs))
	     (x (list-ref xs i)))
       (unless (some-element
		(lambda (xe)
		 (apply
		  constraint
		  (map-indexed (lambda (x j) (if (= j i) xe (gm-binding x))) xs)))
		x)
	(bottom)))))
    x))
  xs))

(define (gm-assert-constraint-vp! constraint xs)
 (for-each
  (lambda (x)
   (attach-demon!
    (lambda ()
     (when (one (lambda (x) (not (gm-bound? x))) xs)
      (let* ((i (position-if (lambda (x) (not (gm-bound? x))) xs))
	     (x (list-ref xs i)))
       (when (one-element
	      (lambda (xe)
	       (apply
		constraint
		(map-indexed (lambda (x j) (if (= j i) xe (gm-binding x))) xs)))
	      x)
	(restrict-distribution!
	 x
	 (the-element
	  (lambda (xe)
	   (apply constraint
		  (map-indexed (lambda (x j) (if (= j i) xe (gm-binding x))) xs)))
	  x))))))
    x))
  xs))

(define (gm-assert-constraint-gfc! constraint xs)
 (for-each
  (lambda (x)
   (attach-demon!
    (lambda ()
     (when (one (lambda (x) (not (gm-bound? x))) xs)
      (let* ((i (position-if (lambda (x) (not (gm-bound? x))) xs))
	     (x (list-ref xs i)))
       (restrict-distribution!
	x
	(the-elements
	 (lambda (xe)
	  (apply constraint
		 (map-indexed (lambda (x j) (if (= j i) xe (gm-binding x))) xs)))
	 x)))))
    x))
  xs))

(define (gm-assert-constraint-ac! constraint ds)
 (for-each-indexed
  (lambda (d k)
   (attach-demon!
    (lambda ()
     (for-each-indexed
      (lambda (d i)
       (unless (= i k)
	(restrict-distribution!
	 d
	 (the-elements
	  (lambda (x)
	   (let loop ((ds ds) (xs '()) (j 0))
	    (if (null? ds)
		(apply constraint (reverse xs))
		(if (= j i)
		    (loop (rest ds) (cons x xs) (+ j 1))
		    (some (lambda (pair)
			   (loop (rest ds) (cons (car pair) xs) (+ j 1)))
			  (distribution-variable-distribution (first ds)))))))
	  d))))
      ds))
    d))
  ds))

(define (gm-assert-constraint! constraint . distribution-variables)
 (case *gm-strategy*
  ((efd) (gm-assert-constraint-efd! constraint distribution-variables))
  ((fc) (gm-assert-constraint-fc! constraint distribution-variables))
  ((vp) (gm-assert-constraint-fc! constraint distribution-variables)
   (gm-assert-constraint-vp! constraint distribution-variables))
  ((gfc) (gm-assert-constraint-gfc! constraint distribution-variables))
  ((ac) (gm-assert-constraint-ac! constraint distribution-variables))
  (else (error "Unrecognized strategy"))))
)
