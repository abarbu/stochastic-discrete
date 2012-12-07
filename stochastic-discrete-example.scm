
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
