(module stochastic-discrete *
(import chicken scheme)
(require-extension nondeterminism define-structure srfi-1 traversal)

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

)
