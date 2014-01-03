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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Definition of dynamic-vectors (Dvector) and basic functions
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
  (cond ((< index (Dvector-valid-entries Dvector))
	 (vector-set! (Dvector-body Dvector) index value))
	((< index (vector-length (Dvector-body Dvector)))
	 (vector-set! (Dvector-body Dvector) index value)
	 (Dvector-valid-entries-set! Dvector (+ index 1)))
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
		new-body counter (vector-ref body counter)))))))

(define Dvector-length Dvector-valid-entries)

(define (Dvector-append Dvector value)
  (Dvector-set! Dvector (Dvector-valid-entries Dvector) value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Basic definitions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-structure point x y)

(define-structure segment from to)

(define-structure triangle 
		  pvertex-1 pvertex-2 pvertex-3
		  pointer-1-2 pointer-2-3 pointer-3-1)

(define no-ptriangle -1)
		  
(define-structure polygon limits vertices segments holes triangulation)

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

(define pi 3.14159265358969)

(define (one-third-more x)
  (+ 2 (inexact->exact (round (* 1.5 x)))))

(define (add-squares x y)
  (declare (flonum))
  (+ (* x x) (* y y)))

(define (determinant a b c
		     d e f 
		     g h i)
  (declare (flonum))
  (- (+ (* a e i) (* b f g) (* c d h))
     (+ (* c e g) (* b d i) (* a f h))))

(define (distance point-1 point-2)
  (declare (flonum))
  (let ((x1 (point-x point-1))
	(y1 (point-y point-1))
	(x2 (point-x point-2))
	(y2 (point-y point-2)))
  (sqrt (add-squares (- x1 x2) (- y1 y2)))))

(define (midpoint point-1 point-2)
  (declare (flonum))
  (let ((x1 (point-x point-1))
	(y1 (point-y point-1))
	(x2 (point-x point-2))
	(y2 (point-y point-2)))
    (make-point (/ (+ x1 x2) 2) (/ (+ y1 y2) 2))))

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
    (and (zero? (functional (point-x point) (point-y point)))
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
  (let ((circle-test (circle-functional vertex-1 vertex-2 vertex-3))
        (orientation-test (segment-functional vertex-2 vertex-3)))
    (positive? 
      (* (circle-test (point-x point) 
		      (point-y point)) 
	 (orientation-test (point-x vertex-1) 
			   (point-y vertex-1))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Functions of segments 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (opposite-vertex-to-segment pvertex-from pvertex-to triangle)
  (declare (fixnum))
  (let ((pvertex-1 (triangle-pvertex-1 triangle))
	(pvertex-2 (triangle-pvertex-2 triangle))
	(pvertex-3 (triangle-pvertex-3 triangle)))
    (cond ((or (and (= pvertex-1 pvertex-from)
		    (= pvertex-2 pvertex-to))
	       (and (= pvertex-1 pvertex-to)
		    (= pvertex-2 pvertex-from)))
	   pvertex-3)
	  ((or (and (= pvertex-2 pvertex-from)
		    (= pvertex-3 pvertex-to))
	       (and (= pvertex-2 pvertex-to)
		    (= pvertex-3 pvertex-from)))
	   pvertex-1)
	  ((or (and (= pvertex-3 pvertex-from)
		    (= pvertex-1 pvertex-to))
	       (and (= pvertex-3 pvertex-to)
		    (= pvertex-1 pvertex-from)))
	   pvertex-2))))

(define (pointer-associated-to-segment 
	  pvertex-from pvertex-to triangle)
  (declare (fixnum))
  (let ((pvertex-1 (triangle-pvertex-1 triangle))
	(pvertex-2 (triangle-pvertex-2 triangle))
	(pvertex-3 (triangle-pvertex-3 triangle)))
    (cond ((or (and (= pvertex-1 pvertex-from)
		    (= pvertex-2 pvertex-to))
	       (and (= pvertex-1 pvertex-to)
		    (= pvertex-2 pvertex-from)))
	   (triangle-pointer-1-2 triangle))
	  ((or (and (= pvertex-2 pvertex-from)
		    (= pvertex-3 pvertex-to))
	       (and (= pvertex-2 pvertex-to)
		    (= pvertex-3 pvertex-from)))
	   (triangle-pointer-2-3 triangle))
	  ((or (and (= pvertex-3 pvertex-from)
		    (= pvertex-1 pvertex-to))
	       (and (= pvertex-3 pvertex-to)
		    (= pvertex-1 pvertex-from)))
	   (triangle-pointer-3-1 triangle)))))

(define (segment-support pvertex-from pvertex-to triangulation)
  (declare (fixnum))
  (let ((size (Dvector-length triangulation)))
    (let loop ((triangle-entry 0))
      (if (= triangle-entry size)
	(error "input edge not present in triangulation")
	(let* ((triangle (Dvector-ref triangulation triangle-entry))
	       (pvertex-1 (triangle-pvertex-1 triangle))
	       (pvertex-2 (triangle-pvertex-2 triangle))
	       (pvertex-3 (triangle-pvertex-3 triangle)))
	  (cond ((or (and (= pvertex-1 pvertex-from)
			  (= pvertex-2 pvertex-to))
		     (and (= pvertex-1 pvertex-to)
			  (= pvertex-2 pvertex-from)))
		 (list triangle-entry
		       (triangle-pointer-1-2 triangle)))
		((or (and (= pvertex-2 pvertex-from)
			  (= pvertex-3 pvertex-to))
		     (and (= pvertex-2 pvertex-to)
			  (= pvertex-3 pvertex-from)))
		 (list triangle-entry
		       (triangle-pointer-2-3 triangle)))
		((or (and (= pvertex-3 pvertex-from)
			  (= pvertex-1 pvertex-to))
		     (and (= pvertex-3 pvertex-to)
			  (= pvertex-1 pvertex-from)))
		 (list triangle-entry
		       (triangle-pointer-3-1 triangle)))
		(else (loop (+ triangle-entry 1)))))))))

(define (is-segment-in-segments? pvertex-from pvertex-to vertices)
  (declare (flonum))
  (or (= pvertex-from 
	 (point-pointer (Dvector-ref vertices pvertex-to)))
      (= pvertex-to 
	 (point-pointer (Dvector-ref vertices pvertex-from)))))

(define (change-pointers 
	  ptriangle old-pointer new-pointer triangulation)
  (declare (fixnum))
  (if (not (= ptriangle no-ptriangle))
    (let ((triangle (Dvector-ref triangulation ptriangle)))
      (cond ((= old-pointer (triangle-pointer-1-2 triangle))
	     (triangle-pointer-1-2-set! triangle new-pointer))
	    ((= old-pointer (triangle-pointer-2-3 triangle))
	     (triangle-pointer-2-3-set! triangle new-pointer))
	    ((= old-pointer (triangle-pointer-3-1 triangle))
	     (triangle-pointer-3-1-set! triangle new-pointer))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Angle functions and related.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (wedge-angle side-point-1 corner-point side-point-2)
  (declare (flonum))
  (let* ((direction-1 (+ (- (point-x side-point-1)
			    (point-x corner-point))
			 (* +i (- (point-y side-point-1)
				  (point-y corner-point)))))
	 (direction-2 (+ (- (point-x side-point-2)
			    (point-x corner-point))
			 (* +i (- (point-y side-point-2)
				  (point-y corner-point)))))
	 (pre-angle (* 180 (/ 1 pi) 
		       (abs (- (angle direction-2) 
			       (angle direction-1))))))
    (if (> pre-angle 180) (- 360 pre-angle) pre-angle)))
					  
(define (polygon-angles vertices triangulation)
  (declare (flonum))
  (let ((size (Dvector-length triangulation)))
    (let loop ((triangle-entry 0))
      (if (= triangle-entry size)
	(newline)
	(let* ((triangle (Dvector-ref triangulation triangle-entry))
	       (pvertex-1 (triangle-pvertex-1 triangle))
	       (pvertex-2 (triangle-pvertex-2 triangle))
	       (pvertex-3 (triangle-pvertex-3 triangle))
	       (vertex-1 (Dvector-ref vertices pvertex-1))
	       (vertex-2 (Dvector-ref vertices pvertex-2))
	       (vertex-3 (Dvector-ref vertices pvertex-3)))
	  (display "triangle #")
	  (display triangle-entry)
	  (display ":\t( ")
	  (display (wedge-angle vertex-3 vertex-1 vertex-2))
	  (display " ")
	  (display (wedge-angle vertex-1 vertex-2 vertex-3))
	  (display " ")
	  (display (wedge-angle vertex-2 vertex-3 vertex-1))
	  (display " )")
	  (newline)
	  (loop (+ 1 triangle-entry)))))))
	  
(define (polygon->file polygon file)

  (define (write-polygon polygon)
    (let* ((vertices (polygon-vertices polygon))
	   (triangulation (polygon-triangulation polygon))
	   (size (Dvector-length triangulation)))
      (let loop ((triangle-entry 0))
	(if (= triangle-entry size)
	  (newline)
	  (let* ((triangle (Dvector-ref triangulation triangle-entry))
		 (pvertex-1 (triangle-pvertex-1 triangle))
		 (pvertex-2 (triangle-pvertex-2 triangle))
		 (pvertex-3 (triangle-pvertex-3 triangle))
		 (vertex-1 (Dvector-ref vertices pvertex-1))
		 (vertex-2 (Dvector-ref vertices pvertex-2))
		 (vertex-3 (Dvector-ref vertices pvertex-3)))
	    (display (point-x vertex-1))
	    (display " ")
	    (display (point-y vertex-1))
	    (newline)
	    (display (point-x vertex-2))
	    (display " ")
	    (display (point-y vertex-2))
	    (newline)
	    (display (point-x vertex-3))
	    (display " ")
	    (display (point-y vertex-3))
	    (newline)
	    (display (point-x vertex-1))
	    (display " ")
	    (display (point-y vertex-1))
	    (newline)
	    (newline)
	    (loop (+ triangle-entry 1)))))))
  
  (with-output-to-file file (lambda () (write-polygon polygon))))

