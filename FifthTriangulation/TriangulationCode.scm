(declare (standard-bindings) (extended-bindings) (block) (not safe))

(define-macro (FLOAT . rest) 
	      '(let ()
		 (declare (flonum))
		 ,@rest))

(define-macro (FIX . rest)
	      '(let ()
		 (declare (fixnum))
		 ,@rest))

(define-macro (GENERIC . rest)
	      '(let ()
		 (declare (generic))
		 ,@rest))

(define pi  3.14159265358969)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Definition of dynamic-vector (Dvector) and related basic functions
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-structure Dvector body valid-entries)

(define (build-Dvector n #!optional initial-value)
  (make-Dvector 
    (make-vector (one-third-more n) initial-value)
    0))

(define (Dvector-ref Dvector index)
  (declare (fixnum))
  (if (< index (Dvector-valid-entries Dvector))
    (vector-ref (Dvector-body Dvector) index)
    (error "out of range")))

(define (Dvector-set! Dvector index value)
  (declare (fixnum))
  (let ((size (Dvector-valid-entries Dvector)))
    (cond ((< index size)
	   (vector-set! (Dvector-body Dvector) index value))
	  ((< index (vector-length (Dvector-body Dvector)))
	   (vector-set! (Dvector-body Dvector) index value)
	   (Dvector-valid-entries-set! Dvector (+  size 1)))
	  (else
	    (let ((new-body
		    (make-vector 
		      (one-third-more (+ index 1))
		      #f))
		  (body (Dvector-body Dvector))
		  (size (Dvector-valid-entries Dvector)))
	      (do ((counter (- size 1) (- counter 1)))
		((< counter 0)
		 (Dvector-body-set! Dvector new-body)
		 (vector-set! new-body index value)
		 (Dvector-valid-entries-set! Dvector (+ size 1)))
		(vector-set! 
		  new-body counter (vector-ref body counter))))))))

(define Dvector-length Dvector-valid-entries)

(define (Dvector-append Dvector value)
  (Dvector-set! Dvector (Dvector-valid-entries Dvector) value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Geometric definitions and basic functions 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-structure point x y pointers triangles)

(define-structure triangle 
		  vertex-1 vertex-2 vertex-3
		  triangle-1-2 triangle-2-3 triangle-3-1)

(define no-point    #f)
(define no-triangle #f)

(define-structure polygon limits vertices holes triangulation)

(define (determinant a b c
		     d e f 
		     g h i)
  (declare (flonum))
  (- (+ (* a e i) (* b f g) (* c d h))
     (+ (* c e g) (* b d i) (* a f h))))

;; SOLVE-SYSTEM solves the following linear system of equations.
;;  a*x + b*y = c
;;  d*x + e*y = f
;; it assumes that the system has solution, since it will be called
;; by procedures that already know that a solution exists.
(define (solve-system a b c
		      d e f)
  (declare (flonum))
  (list (/ (- (* b f) (* c e))
	   (- (* b d) (* a e 1.0)))
	(/ (- (* c d) (* a f))
	   (- (* b d) (* a e 1.0)))))

(define (one-third-more x)
  (+ 2 (inexact->exact (round (* 1.5 x)))))

(define (add-squares x y)
  (declare (flonum))
  (+ (* x x) (* y y)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Functions of points
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (equal-points? point-1 point-2)
  (declare (flonum))
  (and (= (point-x point-1) (point-x point-2))
       (= (point-y point-1) (point-y point-2))))

(define (distance point-1 point-2)
  (declare (flonum))
  (let ((x1 (point-x point-1))
	(y1 (point-y point-1))
	(x2 (point-x point-2))
	(y2 (point-y point-2)))
  (sqrt (add-squares (- x1 x2) (- y1 y2)))))

(define (midpoint point-1 point-2)
  (declare (flonum))
  (let* ((x1 (point-x point-1))
	 (y1 (point-y point-1))
	 (x2 (point-x point-2))
	 (y2 (point-y point-2))
	 (point (make-point (/ (+ x1 x2) 2.0) (/ (+ y1 y2) 2.0) 
			    (build-Dvector 0 no-point) 
			    (build-Dvector 0 no-triangle))))
    (Dvector-append (point-pointers point) point-2)
    point))

(define (is-a-segment? point-1 point-2)
  
  (define (is-point-in-pointers? fixed-point query-point)
    (let* ((pointers (point-pointers fixed-point))
	   (size (Dvector-length pointers)))
      (let loop ((index 0))
	(if (= index size)
	  #f
	  (if (equal-points? query-point (Dvector-ref pointers index))
	    #t
	    (loop (+ index 1)))))))

  (or (is-point-in-pointers? point-1 point-2)
      (is-point-in-pointers? point-2 point-1)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Functions of triangles
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (equal-triangles? triangle-1 triangle-2)
  (declare (flonum))
  (cond ((and (not (triangle? triangle-1))
	      (not (triangle? triangle-2)))
	 #t)
	((or (and (triangle? triangle-1) (not (triangle? triangle-2)))
	     (and (triangle? triangle-2) (not (triangle? triangle-1))))
	 #f)
	(else
	  (let ((vertex-1-1 (triangle-vertex-1 triangle-1))
		(vertex-2-1 (triangle-vertex-2 triangle-1))
		(vertex-3-1 (triangle-vertex-3 triangle-1))
		(vertex-1-2 (triangle-vertex-1 triangle-2))
		(vertex-2-2 (triangle-vertex-2 triangle-2))
		(vertex-3-2 (triangle-vertex-3 triangle-2)))
	    (or (and (equal-points? vertex-1-1 vertex-1-2)
		     (equal-points? vertex-2-1 vertex-2-2)
		     (equal-points? vertex-3-1 vertex-3-2))
		(and (equal-points? vertex-1-1 vertex-1-2)
		     (equal-points? vertex-2-1 vertex-3-2)
		     (equal-points? vertex-3-1 vertex-2-2))
		(and (equal-points? vertex-1-1 vertex-2-2)
		     (equal-points? vertex-2-1 vertex-3-2)
		     (equal-points? vertex-3-1 vertex-1-2))
		(and (equal-points? vertex-1-1 vertex-2-2)
		     (equal-points? vertex-2-1 vertex-1-2)
		     (equal-points? vertex-3-1 vertex-3-2))
		(and (equal-points? vertex-1-1 vertex-3-2)
		     (equal-points? vertex-2-1 vertex-1-2)
		     (equal-points? vertex-3-1 vertex-2-2))
		(and (equal-points? vertex-1-1 vertex-3-2)
		     (equal-points? vertex-2-1 vertex-2-2)
		     (equal-points? vertex-3-1 vertex-1-2)))))))
					      
(define (change-pointers
	  object old-object new-object)
  (cond ((triangle? object)
	 (cond ((equal-triangles? old-object 
				  (triangle-triangle-1-2 object)) 
		(triangle-triangle-1-2-set! object new-object)) 
	       ((equal-triangles? old-object 
				  (triangle-triangle-2-3 object)) 
		(triangle-triangle-2-3-set! object new-object)) 
	       ((equal-triangles? old-object 
				  (triangle-triangle-3-1 object)) 
		(triangle-triangle-3-1-set! object new-object)) 
	       (else (error "input triangles not adjacent"))))
	((and (point? object) (triangle? old-object))
	 (let* ((triangles (point-triangles object))
		(size (Dvector-length triangles)))
	   (let loop ((index 0))
	     (if (= index size)
	       (error "input point not in triangle")
	       (if (equal-triangles? (Dvector-ref triangles index)
				     old-object)
		 (Dvector-set! triangles index new-object)
		 (loop (+ index 1)))))))
	((and (point? object) (point? old-object))
	 (let* ((segments (point-pointers object))
		(size (Dvector-length segments)))
	   (let loop ((index 0))
	     (if (= index size)
	       (error "input points not connected")
	       (if (equal-points? (Dvector-ref segments index) old-object)
		 (Dvector-set! segments index new-object)
		 (loop (+ index 1)))))))))
	       
(define (opposite-vertex-to-segment vertex-from vertex-to triangle)
  (declare (flonum))
  (let ((vertex-1 (triangle-vertex-1 triangle))
	(vertex-2 (triangle-vertex-2 triangle))
	(vertex-3 (triangle-vertex-3 triangle)))
  (cond ((or (and (equal-points? vertex-1 vertex-from)
		  (equal-points? vertex-2 vertex-to))
	     (and (equal-points? vertex-1 vertex-to)
		  (equal-points? vertex-2 vertex-from)))
	 vertex-3)
	((or (and (equal-points? vertex-2 vertex-from)
		  (equal-points? vertex-3 vertex-to))
	     (and (equal-points? vertex-2 vertex-to)
		  (equal-points? vertex-3 vertex-from)))
	 vertex-1)
	((or (and (equal-points? vertex-3 vertex-from)
		  (equal-points? vertex-1 vertex-to))
	     (and (equal-points? vertex-3 vertex-to)
		  (equal-points? vertex-1 vertex-from)))
	 vertex-2)
	(else (error "input vertice(s) not in triangle")))))

(define (triangle-associated-to-segment vertex-from vertex-to triangle)
  (declare (flonum))
  (let ((vertex-1 (triangle-vertex-1 triangle))
	(vertex-2 (triangle-vertex-2 triangle))
	(vertex-3 (triangle-vertex-3 triangle)))
  (cond ((or (and (equal-points? vertex-1 vertex-from)
		  (equal-points? vertex-2 vertex-to))
	     (and (equal-points? vertex-1 vertex-to)
		  (equal-points? vertex-2 vertex-from)))
	 (triangle-triangle-1-2 triangle))
	((or (and (equal-points? vertex-2 vertex-from)
		  (equal-points? vertex-3 vertex-to))
	     (and (equal-points? vertex-2 vertex-to)
		  (equal-points? vertex-3 vertex-from)))
	 (triangle-triangle-2-3 triangle))
	((or (and (equal-points? vertex-3 vertex-from)
		  (equal-points? vertex-1 vertex-to))
	     (and (equal-points? vertex-3 vertex-to)
		  (equal-points? vertex-1 vertex-from)))
	 (triangle-triangle-3-1 triangle))
	(else (error "input vertice(s) not in triangle")))))

(define (is-triangle-vertex? vertex triangle)
  (declare (flonum))
  (if (triangle? triangle)
    (or (equal-points? vertex (triangle-vertex-1 triangle))
	(equal-points? vertex (triangle-vertex-2 triangle))
	(equal-points? vertex (triangle-vertex-3 triangle)))
    #f))

(define (triangle-support point triangulation)
  (let ((size (Dvector-length triangulation)))
    (let loop ((index 0))
      (if (= index size)
	(error "point not in domain")
	(let* ((triangle (Dvector-ref triangulation index))
	       (vertex-1 (triangle-vertex-1 triangle))
	       (vertex-2 (triangle-vertex-2 triangle))
	       (vertex-3 (triangle-vertex-3 triangle)))
	  (if (or (point-inside-triangle? point vertex-1 vertex-2 vertex-3)
		  (point-inside-segment? point vertex-1 vertex-2)
		  (point-inside-segment? point vertex-2 vertex-3)
		  (point-inside-segment? point vertex-3 vertex-1))
	    triangle
	    (loop (+ index 1))))))))
				 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Functionals and applications to:
;;	point-inside-triangle
;;	point-inside-segment
;;	point-inside-circumcircle
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; SEGMENT-FUNCTIONAL:  Given the vertices of a segment, consider a
;; (canonical) equation of the line that goes through them (one in
;; which division is not used, please).  
;;
;;		       f(x,y) = 0
;; 
;; The functional associated to that equation gives us information of 
;; the relative location of a point of the plane with respect to that
;; line.
;; Proof: The way I write the functional can be read as the inner pro-
;;        duct of the direction of the line, and the perpendicular to 
;;        the vector defined by the first vertex and the point (x,y).
;;        This inner product is nothing by the product of the lengths
;;        of the vectors (positive) with the cosine of the convex angle
;;        defined by them.  The sign of the inner product then indicates
;;        that the point belongs to either half plane.  
;; Remark: The bigger the absolute value, the closer to ±90°, and for 
;;        absolute values of this functional very close to zero, we have
;;        almost collinearity.
;;
(define (segment-functional vertex-1 vertex-2)
  (declare (flonum))
  (lambda (x y)
    (+ (* (- (point-x vertex-2) (point-x vertex-1)) 
	  (- y (point-y vertex-1)))
       (* (- (point-x vertex-1) x) 
	  (- (point-y vertex-2) (point-y vertex-1))))))

;; POINT-INSIDE-TRIANGLE?: It checks the sign of the three segment
;; functionals defined by the triangle.  If all of them have the same 
;; sign, we are in; otherwise, out.
(define (point-inside-triangle? point vertex-1 vertex-2 vertex-3)
  (declare (flonum))
  (let ((x (point-x point))
	(y (point-y point))
	(functional-1 (segment-functional vertex-1 vertex-2))
	(functional-2 (segment-functional vertex-2 vertex-3))
	(functional-3 (segment-functional vertex-3 vertex-1)))
    (and (positive? (* (functional-1 x y) (functional-2 x y)))
	 (positive? (* (functional-2 x y) (functional-3 x y))))))

;;POINT-INSIDE-SEGMENT?: Looks for collinearity and encroaching.
;;
(define (point-inside-segment? point vertex-1 vertex-2)
  (declare (flonum))
  (let ((functional (segment-functional vertex-1 vertex-2))
	(center (midpoint vertex-1 vertex-2)))
    (and (= (functional (point-x point) (point-y point)) 0)
	 (< (distance point center) (distance vertex-1 center)))))

;; CIRCLE-FUNCTIONAL:  Given three non-collinear vertices (a triangle),
;; consider the following equation of the circumcircle (division and
;; square-root free!).
;;              |                       |
;;              |  1    x   y  x²+y²    |
;;              |                       |
;;              |  1   x1  y1  x1²+y1²  |
;;     g(x,y) = |                       | = 0
;;              |  1   x2  y2  x2²+y2²  |
;;              |                       |
;;              |  1   x3  y3  x3²+y3²  |
;;              |                       |
;;
(define (circle-functional vertex-1 vertex-2 vertex-3)
  (declare (flonum))
  (let ((x1 (point-x vertex-1)) (y1 (point-y vertex-1))
	(x2 (point-x vertex-2)) (y2 (point-y vertex-2))
	(x3 (point-x vertex-3)) (y3 (point-y vertex-3)))
    (lambda (x y)
      (+ (determinant x1 y1 (add-squares x1 y1)
		      x2 y2 (add-squares x2 y2)
		      x3 y3 (add-squares x3 y3))
	 (* (- x)
	    (determinant 1 y1 (add-squares x1 y1)
			 1 y2 (add-squares x2 y2)
			 1 y3 (add-squares x3 y3)))
	 (* y 
	    (determinant 1 x1 (add-squares x1 y1)
			 1 x2 (add-squares x2 y2)
			 1 x3 (add-squares x3 y3)))
	 (* (- (add-squares x y))
	    (determinant 1 x1 y1
			 1 x2 y2
			 1 x3 y3))))))
	    
;; POINT-INSIDE-CIRCUMCIRCLE?: The circle-functional is not enough by
;; itself to decide if a query point is inside, since the sign of the
;; determinant depends heavily on the orientation of the vertices; 
;; therefore, we must consider if the points are ordered clock or
;; counter-clockwise.
;;
(define (point-inside-circumcircle? point vertex-1 vertex-2 vertex-3)
  (declare (flonum))
  (let ((circle-test (circle-functional vertex-1 vertex-2 vertex-3))
        (orientation-test (segment-functional vertex-2 vertex-3)))
    (positive? 
      (* (circle-test (point-x point) 
		      (point-y point)) 
	 (orientation-test (point-x vertex-1) 
			   (point-y vertex-1))))))
