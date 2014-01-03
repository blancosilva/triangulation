(define p1  (make-point 0 0 'p1))  (define p8  (make-point 1 2 'p8))
(define p2  (make-point 1 -1 'p2)) (define p9  (make-point 0 2 'p9))
(define p3  (make-point 2 0 'p3))  (define p10 (make-point  -1  3.5 'p10))
(define p4  (make-point 2 2 'p4))  (define p11 (make-point  -3 3 'p11))
(define p5  (make-point 4 0 'p5))  (define p12 (make-point -2 2.5 'p12))
(define p6  (make-point 4 4 'p6))  (define p13 (make-point -4 0.5 'p13))
(define p7  (make-point 2 4 'p7))  (define p14 (make-point 5 2 'p14))
(define polygon (list->vector (list p1 p2 p3 p4 p5 p14 p6 p7 p8 p9 p13)))

(define r1 (make-point 2 2 'r1))
(define r2 (make-point 4 0 'r2))
(define r3 (make-point 6 0 'r3))
(define r4 (make-point 7 1 'r4))
(define r5 (make-point 9 9 'r5))
(define r6 (make-point 8 14 'r6))
(define r7 (make-point 7 6 'r7))
(define r8 (make-point 6 10 'r8))
(define r9 (make-point 7 15 'r9))
(define r10 (make-point 0 16 'r10))
(define r11 (make-point 0 9 'r11))
(define r12 (make-point 5 5 'r12))
(define r13 (make-point 5 4 'r13))
(define mesh (list->vector (list r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13)))

(define q1 (make-point 0 0 'q1))
(define q2 (make-point 15 0 'q2))
(define q3 (make-point 1 1 'q3))   
(define q4 (make-point 3 10 'q4))
(define q5 (make-point -13 4 'q5)) 
(define q6 (make-point -13 -13 'q6))
(define q7 (make-point 0 -20 'q7)) 
(define q7b (make-point -5 -15 'q7b))
(define q8 (make-point -7 -3 'q8))
(define q9 (make-point -7 -10 'q9))
(define q10 (make-point -10 -10 'q10))
(define q11 (make-point -10 0 'q11))
(define q12 (make-point -5 3 'q12))
(define q13 (make-point 0 3 'q13))
(define q14 (make-point -3 1 'q14))
(define devil (list->vector (list q1 q2 q3 q4 q5 q6 q7 q7b q8 q9 q10 q11 q12 q13 q14)))

(define s1 (make-point 15 -1 's1))
(define s2 (make-point 0 0 's2))
(define s3 (make-point 15 1 's3))
(define s4 (make-point 15 20 's4))
(define s5 (make-point 16 0 's5))
(define s6 (make-point 15 -20 's6))
(define teaser (list->vector (list s1 s6 s5 s4 s3 s2)))

(define t1 (make-point 1 1 't1))
(define t2 (make-point 11 2 't2))
(define t3 (make-point 21 0 't3))
(define t4 (make-point 11 -2 't4))
(define t5 (make-point 1 -1 't5))
(define t6 (make-point 2 -11 't6))
(define t7 (make-point 0 -21 't7))
(define t8 (make-point -2 -11 't8))
(define t9 (make-point -1 -1 't9))
(define t10 (make-point -11 -2 't10))
(define t11 (make-point -21 0 't11))
(define t12 (make-point -11 2 't12))
(define t13 (make-point -1 1 't13))
(define t14 (make-point -2 11 't14))
(define t15 (make-point 0 21 't15))
(define t16 (make-point 2 11 't16))
(define cross (list->vector (list t1  t2  t3  t4  t5  t6  t7  t8
		    t9  t10 t11 t12 t13 t14 t15 t16)))

(define (convex-polygon sides)
  (let loop ((position 0) (result '()))
    (if (= position sides)
	(list->vector result)
	(loop (+ position 1) 
	      (append result 
		      (list 
			(make-point (+ 0.1 
				       (* 10.0 
					  (cos 
					    (* 2 
					       pi 
					       position 
					       (/ 1.0 sides))))) 
				    (+ 0.1 
				       (* 10.0 
					  (sin 
					    (* 2 
					       pi 
					       position 
					       (/ 1.0 sides))))) 
				    position)))))))
