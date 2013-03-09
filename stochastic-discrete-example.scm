
(define (gradient-ascent-F f x0 n eta)
 (if (zero? n)
     (list x0 (f x0) ((gradient-F f) x0))
     (gradient-ascent-F
      f
      (map-vector (lambda (xi gi) (+ xi (* eta gi))) x0 ((gradient-F f) x0))
      (- n 1)
      eta)))

(define (gradient-ascent-R f x0 n eta)
 (if (zero? n)
     (list x0 (f x0) ((gradient-R f) x0))
     (gradient-ascent-R
      f
      (map-vector (lambda (xi gi) (+ xi (* eta gi))) x0 ((gradient-R f) x0))
      (- n 1)
      eta)))

(define (stochastic-scheme-example)
 (gradient-ascent-R
  (lambda (p)
   (reduce-vector
    *
    (map-vector (lambda (observation)
		 (probability (eqv? observation
				    (if (flip (vector-ref p 0))
					0
					(if (flip (vector-ref p 1)) 1 2)))))
		'#(0 1 2 2))
    1))
  '#(0.5 0.5)
  1000
  0.1))

(define (example1)
 (let* ((rain (flip 0.2))
        (sprinkler (if rain (flip 0.01) (flip 0.4)))
        (grass-wet (cond ((and (not sprinkler) (not rain)) (flip 0))
                         ((and (not sprinkler) rain)       (flip 0.8))
                         ((and sprinkler (not rain))       (flip 0.9))
                         ((and sprinkler rain)             (flip 0.99)))))
  (list rain sprinkler grass-wet)))

;; expect to see 0.288 and 0.16038, which normalized
;;  0.16038/(0.288*0.16038)=0.3576
;; agrees with wikipedia
(distribution
 (let ((r (example1)))
  (when (not (last r)) (bottom))
  (first r)))
