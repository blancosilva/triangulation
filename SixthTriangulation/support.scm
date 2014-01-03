;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Geometric definitions and basic functions 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(define (add-squares x y)
  (declare (flonum))
  (+ (* x x) (* y y)))

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
	 (point (make-point (/ (+ x1 x2) 2.0) (/ (+ y1 y2) 2.0))))
    point))

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
(define (segment-functional segment)
  (declare (flonum))
  (let ((vertex-1 (segment-from segment))
	(vertex-2 (segment-to segment)))
    (lambda (x y)
      (+ (* (- (point-x vertex-2) (point-x vertex-1)) 
	    (- y (point-y vertex-1)))
	 (* (- (point-x vertex-1) x) 
	    (- (point-y vertex-2) (point-y vertex-1)))))))

;; POINT-IN-TRIANGLE-CLOSURE?: It checks the sign of the three segment
;; functionals defined by the triangle.  If either all of them have the same
;; sign, or one of them is zero, then we are in; otherwise, out.
(define (point-in-triangle-closure? point triangle)
  (declare (flonum))
  (let ((segment-1 (triangle-segment-1 triangle))
	(segment-2 (triangle-segment-2 triangle))
	(segment-3 (triangle-segment-3 triangle))
	(x (point-x point))
	(y (point-y point)))
    (let ((functional-1 (segment-functional segment-1))
	  (functional-2 (segment-functional segment-2))
	  (functional-3 (segment-functional segment-3)))
      (and (not (negative? (* (functional-1 x y) (functional-2 x y))))
	   (not (negative? (* (functional-2 x y) (functional-3 x y))))
	   (not (negative? (* (functional-3 x y) (functional-1 x y))))))))

;;POINT-INSIDE-SEGMENT?: Looks for collinearity and encroaching.
;;
;; (define (point-inside-segment? point vertex-1 vertex-2)
;;   (declare (flonum))
;;   (let ((functional (segment-functional vertex-1 vertex-2))
;; 	(center (midpoint vertex-1 vertex-2)))
;;     (and (= (functional (point-x point) (point-y point)) 0)
;; 	 (< (distance point center) (distance vertex-1 center)))))

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
(define (circle-functional triangle)
  (declare (flonum))
  (let ((vertex-1 (triangle-vertex-1 triangle))
	(vertex-2 (triangle-vertex-2 triangle))
	(vertex-3 (triangle-vertex-3 triangle)))
    (let ((x1 (point-x vertex-1)) 
	  (y1 (point-y vertex-1))
	  (x2 (point-x vertex-2)) 
	  (y2 (point-y vertex-2))
	  (x3 (point-x vertex-3)) 
	  (y3 (point-y vertex-3)))
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
			   1 x3 y3)))))))

;; POINT-INSIDE-CIRCUMCIRCLE?: The circle-functional is not enough by itself to
;; decide if a query point is inside, since the sign of the determinant depends
;; on the orientation of the vertices; therefore, we must consider if the
;; points are ordered clock or counter-clockwise.
;;
(define (point-inside-circumcircle? point triangle)
  (declare (flonum))
  (let ((x0 (point-x point))
	(y0 (point-y point))
	(x1 (point-x (triangle-vertex-1 triangle)))
	(y1 (point-y (triangle-vertex-1 triangle)))
	(circle-test (circle-functional triangle))
	(orientation-test 
	  (segment-functional (triangle-segment-1 triangle))))
    (positive? 
      (* (circle-test x0 y0) 
	 (orientation-test x1 y1)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Extra functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
(define (same-segments? segment-1 segment-2)
  (or (and (equal? (segment-from segment-1) (segment-from segment-2))
	   (equal? (segment-to segment-1) (segment-to segment-2)))
      (and (equal? (segment-from segment-1) (segment-to segment-2))
	   (equal? (segment-to segment-1) (segment-from segment-2)))))

(define (same-triangles? triangle-1 triangle-2)
  (if (and (triangle? triangle-1) (triangle? triangle-2))
    (and (equal? (triangle-vertex-1 triangle-1)
		 (triangle-vertex-1 triangle-2))
	 (equal? (triangle-vertex-2 triangle-1)
		 (triangle-vertex-2 triangle-2))
	 (equal? (triangle-vertex-3 triangle-1)
		 (triangle-vertex-3 triangle-2)))
    (eq? triangle-1 triangle-2)))
  
