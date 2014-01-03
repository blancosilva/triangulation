;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Basic functions and definitions.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (add-squares x y)
  (+ (* x x) (* y y)))

(define (subtract-squares x y)
  (- (* x x) (* y y)))

(define (determinant a b c
		     d e f 
		     g h i)
  (- (+ (* a e i) (* b f g) (* c d h))
     (+ (* c e g) (* b d i) (* a f h))))

(define (distance point-1 point-2)
  (sqrt (add-squares (- (point-x point-1) (point-x point-2))
		     (- (point-y point-1) (point-y point-2)))))

(define (vector-append vector element)
  (list->vector
    (append (vector->list vector)
	    (list element))))

;; SOLVE-SYSTEM solves the following linear system of equations.
;;  a*x + b*y = c
;;  d*x + e*y = f
;; it assumes that the system has solution, since it will be called
;; by procedures that already know that a solution exists.
(define (solve-system a b c
		      d e f)
  (list (/ (- (* b f) (* c e))
	   (- (* b d) (* a e 1.0)))
	(/ (- (* c d) (* a f))
	   (- (* b d) (* a e 1.0)))))

(define pi 3.14159265358969)

(define (integer-part dividend divisor)
  (/ (- dividend (remainder dividend divisor)) divisor))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Structures (with examples).
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; POLYGON-POINT is the main structure.  Such a point is defined by its
;; coordinates, and as it is meant to be an entry of a vector (polygon),
;; it carries information to a boundary pointer.
;; If boundary-pointer = 1, then the point is interior to the polygon.
;; If boundary-pointer = 0, then the point is part of a boundary edge.
;;
(define-structure point x y boundary-pointer)

;; SEGMENT carries two integers pointing to vertices, and a boolean
;; indicating if the segment belongs to the boundary box.
;;
(define-structure segment from to not-a-border?)

;; POLYGON is a structure that contains a vector of points, segments,
;; "holes" information and a triangulation of the vertices.
;;
(define-structure polygon vertices segments holes triangulation)

;; Some examples:
;;
;; A point by itself, no polygon to refer to.
;;
;; (define origin (make-point 0 0 -1))
;; 
;; A square (without triangulation):
;;
;; (define square-vertices
;;   (list->vector
;;     (list
;;       (make-point 0 0 0)
;;       (make-point 1 0 0)
;;       (make-point 1 1 0)
;;       (make-point 0 1 0))))
;;
;; (define square-segments
;;   (list->vector (list (make-segment 0 1 #t)
;; 		         (make-segment 1 2 #t)
;; 		         (make-segment 2 3 #t)
;; 		         (make-segment 3 0 #t))))
;; 
;; (define square-hole
;;   (list->vector '()))
;; 
;; (define square 
;;   (make-polygon 
;;     square-vertices square-segments 
;;     square-hole (list->vector '())))
;;
;;           v3 +----<----+ v2
;;              |/////////|
;;              |/////////|
;;              v/////////^
;;              |/////////|
;;              |/////////|
;;           v0 +---->----+ v1
;; 
;; A square with a "hole" inside (no triangulation):
;; 
;; (define holed-square-vertices
;;   (list->vector
;;     (list
;;       (make-point 0 0 0)
;;       (make-point 3 0 0)
;;       (make-point 3 3 0)
;;       (make-point 0 3 0)
;;       (make-point 1 1 0)
;;       (make-point 2 1 0)
;;       (make-point 2 2 0)
;;       (make-point 1 2 0))))
;;
;; (define holed-square-segments
;;   (list->vector (list (make-segment 0 1 #t)
;; 		         (make-segment 1 2 #t)
;; 		         (make-segment 2 3 #t)
;; 		         (make-segment 3 0 #t)
;; 		         (make-segment 4 5 #t)
;; 		         (make-segment 5 6 #t)
;; 		         (make-segment 6 7 #t)
;; 		         (make-segment 7 4 #t))))
;; 
;; (define holed-square-holes 
;;   (vector->list (list (make-point 1.5 1.5 -1))))
;; 		      
;; (define square-with-hole 
;;   (make-polygon
;;     holed-square-vertices holed-square-segments
;;     holed-square-holes (list->vector '())))
;;
;;         v3 +--------<--------+ v2
;;            |/////////////////|
;;            |///+----<----+///|  
;;            |///|v7     v6|///|
;;            |///|         |///|
;;            v///v         ^///^
;;            |///|         |///|
;;            |///|v4     v5|///|
;;            |///+---->----+///|
;;            |/////////////////|
;;         v0 +-------->--------+ v1

;; TRIANGLE is another structure meant to be part of a vector; thus,
;; besides vertices, it carries pointers to the adjacent triangles,
;; and a key to its status 
;; (=0 interior, =1 exterior).
;; A pointer to no triangle is a -1.
;;
(define-structure triangle 
		  vertex-1 vertex-2 vertex-3
		  pointer-1-2 pointer-2-3 pointer-3-1 status)
(define no-triangle -1)

;; Another example: Brad's simple triangulation.
;;
;; (define square-triangulation
;;   (list->vector
;;     (list
;;       (make-triangle 0 1 3 -1 1 -1 0)
;;       (make-triangle 1 2 3 -1 -1 0 0))))
;;
;; 			   -1
;; 		     v3 +-----+ v2
;; 			|\  1 |    
;; 			| \   |
;; 		     -1 |  \  | -1
;; 			| 0 \ |
;; 			|    \|
;; 		     v0 +-----+ v1
;; 			   -1

;; NAIVE-TRIANGULATION is simply a triangulation taking the first 
;; vertex of the polygon, and the four boundary vertices.
;;
(define (naive-triangulation size)
  (list->vector
    (list
      (make-triangle (- size 4) (- size 3) 0 -1 1 3 1)
      (make-triangle (- size 3) (- size 2) 0 -1 2 0 1)
      (make-triangle (- size 2) (- size 1) 0 -1 3 1 1)
      (make-triangle (- size 1) (- size 4) 0 -1 0 2 1))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Functionals and applications to point-inside-triangle/circumcircle
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
  (lambda (x y)
    (+ (* (- (point-x vertex-2) (point-x vertex-1)) 
	  (- y (point-y vertex-1)))
       (* (- (point-x vertex-1) x) 
	  (- (point-y vertex-2) (point-y vertex-1))))))

;; POINT-INSIDE-TRIANGLE?: It checks the sign of the three segment
;; functionals defined by the triangle.  If all of them have the same 
;; sign, we are in; otherwise, out.
(define (point-inside-triangle? point vertex-1 vertex-2 vertex-3)
  (let ((x (point-x point))
	(y (point-y point))
	(functional-1 (segment-functional vertex-1 vertex-2))
	(functional-2 (segment-functional vertex-2 vertex-3))
	(functional-3 (segment-functional vertex-3 vertex-1)))
    (and (positive? (* (functional-1 x y) (functional-2 x y)))
	 (positive? (* (functional-2 x y) (functional-3 x y))))))

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
;;  Segment functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; SEGMENT-MIDPOINT: 
;;
(define (segment-midpoint vertex-1 vertex-2)
  (make-point (/ (+ (point-x vertex-1) (point-x vertex-2)) 2.0)
	      (/ (+ (point-y vertex-1) (point-y vertex-2)) 2.0)
	      0))

(define (segment-encroached-upon-by? point vertex-1 vertex-2)
  (< (distance point    (segment-midpoint vertex-1 vertex-2))
     (distance vertex-1 (segment-midpoint vertex-1 vertex-2))))

(define (first-segment-encroached-upon-by point segments vertices)
  (let loop ((segment-entry 0))
    (if (= segment-entry (vector-length segments))
      -1 
      (let* ((segment (vector-ref segments segment-entry))
	     (vertex-1 (vector-ref vertices (segment-from segment)))
	     (vertex-2 (vector-ref vertices (segment-to segment))))
	(if (and (segment-encroached-upon-by? point vertex-1 vertex-2)
		 (segment-not-a-border? segment))
	  segment-entry
	  (loop (+ 1 segment-entry)))))))

(define (is-segment-in-triangulation? 
	  point-entry-1 point-entry-2 triangulation)
  (let loop ((triangle-entry 0) (result #f))
    (if (= triangle-entry (vector-length triangulation))
      result
      (let* ((triangle (vector-ref triangulation triangle-entry))
	     (pvertex-1 (triangle-vertex-1 triangle))
	     (pvertex-2 (triangle-vertex-2 triangle))
	     (pvertex-3 (triangle-vertex-3 triangle)))
	(loop (+ 1 triangle-entry)
	      (or result 
		  (and (= pvertex-1 point-entry-1)
		       (= pvertex-2 point-entry-2))
		  (and (= pvertex-1 point-entry-2)
		       (= pvertex-2 point-entry-1))
		  (and (= pvertex-2 point-entry-1)
		       (= pvertex-3 point-entry-2))
		  (and (= pvertex-2 point-entry-2)
		       (= pvertex-3 point-entry-1))
		  (and (= pvertex-3 point-entry-1)
		       (= pvertex-1 point-entry-2))
		  (and (= pvertex-3 point-entry-2)
		       (= pvertex-1 point-entry-1))))))))

(define (segment-with-vertices pvertex-1 pvertex-2 segments)
  (let loop ((counter 0))
    (if (= counter (vector-length segments))
      (error "No edge with input vertices")
      (let ((segment (vector-ref segments counter)))
	(if (or (and (= pvertex-1 (segment-from segment))
		     (= pvertex-2 (segment-to segment)))
		(and (= pvertex-1 (segment-to segment))
		     (= pvertex-2 (segment-from segment))))
	  counter
	  (loop (+ counter 1)))))))

(define (is-segment-in-segments? segment segments)
  (let loop ((counter 0) (result #f))
    (if (= counter (vector-length segments))
      result
      (let ((current-segment (vector-ref segments counter)))
	(loop (+ counter 1) 
	      (or result 
		  (and (= (segment-from segment) 
			  (segment-from current-segment))
		       (= (segment-to segment)
			  (segment-to current-segment)))
		  (and (= (segment-from segment)
			  (segment-to current-segment))
		       (= (segment-to segment)
			  (segment-from current-segment)))))))))

;; TRIANGLES-PER-EDGE finds the pointers to the two triangles that
;; share the input edge.
;;
(define (triangles-per-edge segment triangulation)
  (let loop ((triangle-entry 0))
    (if (= triangle-entry (vector-length triangulation))
      (error "edge not present in triangulation")
      (let* ((triangle (vector-ref triangulation triangle-entry))
	     (pvertex-1 (triangle-vertex-1 triangle))
	     (pvertex-2 (triangle-vertex-2 triangle))
	     (pvertex-3 (triangle-vertex-3 triangle)))
	(cond ((or (and (= pvertex-1 (segment-from segment))
			(= pvertex-2 (segment-to segment)))
		   (and (= pvertex-1 (segment-to segment))
			(= pvertex-2 (segment-from segment))))
	       (list triangle-entry
		     (triangle-pointer-1-2 triangle)))
	      ((or (and (= pvertex-2 (segment-from segment))
			(= pvertex-3 (segment-to segment)))
		   (and (= pvertex-2 (segment-to segment))
			(= pvertex-3 (segment-from segment))))
	       (list triangle-entry
		     (triangle-pointer-2-3 triangle)))
	      ((or (and (= pvertex-3 (segment-from segment))
			(= pvertex-1 (segment-to segment)))
		   (and (= pvertex-3 (segment-to segment))
			(= pvertex-1 (segment-from segment))))
	       (list triangle-entry
		     (triangle-pointer-3-1 triangle)))
	      (else (loop (+ 1 triangle-entry))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Triangle functions.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (interior-triangle? triangle)
  (= 0 (triangle-status triangle)))

(define (border-triangle? triangle segments)
  (let* ((vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (segment-1 (make-segment vertex-1 vertex-2 #t))
	 (segment-2 (make-segment vertex-2 vertex-3 #t))
	 (segment-3 (make-segment vertex-3 vertex-1 #t)))
    (or (and (is-segment-in-segments? segment-1 segments)
	     (is-segment-in-segments? segment-2 segments))
	(and (is-segment-in-segments? segment-2 segments)
	     (is-segment-in-segments? segment-3 segments))
	(and (is-segment-in-segments? segment-3 segments)
	     (is-segment-in-segments? segment-1 segments)))))

;; IS-INTERSECTION-BOUNDARY-EDGE? takes the pointers of two triangles
;; in a triangulation and first computes if some edge in triangle-2 is 
;; a border; if this is not the case, then the answer is #f; otherwise,
;; look in triangle-1 for the two vertices you've found. If they are
;; not there, then again #f. If they are there, #t.
;;
(define (is-intersection-boundary-edge? 
	  ptriangle-1 ptriangle-2 triangulation segments)

  (define (is-vertex-in-triangle? vertex triangle)
    (or (= vertex (triangle-vertex-1 triangle))
	(= vertex (triangle-vertex-2 triangle))
	(= vertex (triangle-vertex-3 triangle))))
  
  (if (or (= ptriangle-1 no-triangle) (= ptriangle-2 no-triangle))
#t
    (let* ((triangle-1 (vector-ref triangulation ptriangle-1))
	   (triangle-2 (vector-ref triangulation ptriangle-2))
	   (vertex-1-2 (triangle-vertex-1 triangle-2))
	   (vertex-2-2 (triangle-vertex-2 triangle-2))
	   (vertex-3-2 (triangle-vertex-3 triangle-2)))
      (or (and (is-segment-in-segments? (make-segment vertex-1-2
						      vertex-2-2 #t)
					segments)
	       (is-vertex-in-triangle? vertex-1-2 triangle-1)
	       (is-vertex-in-triangle? vertex-2-2 triangle-1))
	  (and (is-segment-in-segments? (make-segment vertex-1-2
						      vertex-3-2 #t)
					segments)
	       (is-vertex-in-triangle? vertex-1-2 triangle-1)
	       (is-vertex-in-triangle? vertex-3-2 triangle-1))
	  (and (is-segment-in-segments? (make-segment vertex-3-2
						      vertex-2-2 #t)
					segments)
	       (is-vertex-in-triangle? vertex-3-2 triangle-1)
	       (is-vertex-in-triangle? vertex-2-2 triangle-1))))))

(define (triangle-circumcenter triangle vertices)
  (let* ((vertex-1 (vector-ref vertices (triangle-vertex-1 triangle)))
	 (vertex-2 (vector-ref vertices (triangle-vertex-2 triangle)))
	 (vertex-3 (vector-ref vertices (triangle-vertex-3 triangle)))
	 (x1 (point-x vertex-1)) (y1 (point-y vertex-1))
	 (x2 (point-x vertex-2)) (y2 (point-y vertex-2))
	 (x3 (point-x vertex-3)) (y3 (point-y vertex-3))
	 (intersection 
	   (solve-system 
	     (* 2 (- x1 x2)) 
	     (* 2 (- y1 y2)) 
	     (- (subtract-squares x1 x2) (subtract-squares y2 y1))
	     (* 2 (- x1 x3)) 
	     (* 2 (- y1 y3)) 
	     (- (subtract-squares x1 x3) (subtract-squares y3 y1))))
	 (x (car intersection))
	 (y (cadr intersection)))
    (make-point x y -1)))

;; POINT-SUPPORT finds either the triangle or edge where the input
;; point lies in.  In the second case, it also appends the pointers
;; of the two triangles sharing the edge.
;;
(define (point-support point vertices triangulation)

  (define (point-inside-segment? point vertex-1 vertex-2)
    (and (zero? ((segment-functional vertex-1 vertex-2) 
		 (point-x point) (point-y point)))
	 (segment-encroached-upon-by? point vertex-1 vertex-2)))
  
  (let loop ((triangle-entry 0))
    (if (= triangle-entry (vector-length triangulation))
      (error "input point not in domain.")
      (let* ((triangle (vector-ref triangulation triangle-entry))
	     (pvertex-1 (triangle-vertex-1 triangle))
	     (vertex-1 (vector-ref vertices pvertex-1))
	     (pvertex-2 (triangle-vertex-2 triangle))
	     (vertex-2 (vector-ref vertices pvertex-2))
	     (pvertex-3 (triangle-vertex-3 triangle))
	     (vertex-3 (vector-ref vertices pvertex-3)))
	(if (point-inside-triangle? point vertex-1 vertex-2 vertex-3)
	  (list triangle triangle-entry)
	  (if (point-inside-segment? point vertex-1 vertex-2)
	    (list 
	      (make-segment pvertex-1 pvertex-2 #t)
	      triangle-entry
	      (triangle-pointer-1-2 triangle))
	    (if (point-inside-segment? point vertex-2 vertex-3)
	      (list
		(make-segment pvertex-2 pvertex-3 #t)
		triangle-entry
		(triangle-pointer-2-3 triangle))
	      (if (point-inside-segment? point vertex-3 vertex-1)
		(list 
		  (make-segment pvertex-3 pvertex-1 #t)
		  triangle-entry
		  (triangle-pointer-3-1 triangle))
		(loop (+ 1 triangle-entry))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Angle functions and related.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (wedge-angle side-point-1 corner-point side-point-2)
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
    (if (> 180 pre-angle) (- 360 pre-angle) pre-angle)))
					  
(define (corner-angle corner-triangle vertices segments)
  (let* ((pvertex-1 (triangle-vertex-1 corner-triangle))
	 (vertex-1 (vector-ref vertices pvertex-1))
         (pvertex-2 (triangle-vertex-2 corner-triangle))
	 (vertex-2 (vector-ref vertices pvertex-2))
         (pvertex-3 (triangle-vertex-3 corner-triangle))
	 (vertex-3 (vector-ref vertices pvertex-3))
	 (segment-1 (make-segment pvertex-1 pvertex-2 #t))
	 (segment-2 (make-segment pvertex-2 pvertex-3 #t))
	 (segment-3 (make-segment pvertex-3 pvertex-1 #t)))
    (cond ((and (is-segment-in-segments? segment-1 segments)
		(is-segment-in-segments? segment-2 segments))
	   (wedge-angle vertex-1 vertex-2 vertex-3))
	  ((and (is-segment-in-segments? segment-1 segments)
		(is-segment-in-segments? segment-2 segments))
	   (wedge-angle vertex-2 vertex-3 vertex-1))
	  (else (wedge-angle vertex-3 vertex-1 vertex-2)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Parsers:
;;	.poly-file -> polygon structure
;;	polygon structure -> gnuplot file
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Given a .poly file, obtain a polygon structure.
;;
(define (file->polygon file)
  (let* ((port (open-input-file file))
	 (size (read port))
	 (dim  (read port))
	 (dummy (read port))
	 (dummy (read port))
	 (vertices (make-vector size)))
    (let first-loop ((counter 0))
      (if (= counter size)

	(let ((size (read port))
	      (dummy (read port))
	      (segments (make-vector size)))
	  (let second-loop ((counter 0))
	    (if (= counter size)

	      (let* ((holes (read port))
		     (exterior-points (make-vector holes)))
		(let third-loop ((counter 0))
		  (if (= counter holes) 
		    (time 
		      (include-edges-in-triangulation 
			(make-polygon 
			  vertices 
			  segments 
			  exterior-points 
			  (time
			    (convex-hull-triangulation vertices)))))
		    (begin
		      (vector-set! exterior-points
				   (- (read port) 1)
				   (make-point (read port)
					       (read port)
					       -1))
		      (third-loop (+ counter 1))))))
	      
	      (begin
		(vector-set! segments (- (read port) 1)
			     (make-segment (- (read port) 1)
			     	           (- (read port) 1)
					   (< counter (- size 4))))
		(second-loop (+ counter 1))))))
			     
	(begin
	  (vector-set! vertices (- (read port) 1)
		       (make-point (read port) (read port) 0))
	  (first-loop (+ counter 1)))))))

;; Given a polygon structure, generate a gnuplot file with the trian-
;; gulation.  You must use (get-rid-of-exterior-triangles) BEFORE you
;; use polygon as input; otherwise, you'll get all triangles, exterior
;; and interior.
;;
(define (polygon->file polygon file)
  (define (write-triangulation triangulation polygon-vertices)
    (let loop ((counter 0))
      (if (= counter (vector-length triangulation))
	(newline)
	(if (interior-triangle? (vector-ref triangulation counter))
	  (let* ((triangle (vector-ref triangulation counter))
		 (vertex-1 (vector-ref polygon-vertices
				       (triangle-vertex-1 triangle)))
		 (vertex-2 (vector-ref polygon-vertices
				       (triangle-vertex-2 triangle)))
		 (vertex-3 (vector-ref polygon-vertices
				       (triangle-vertex-3 triangle))))
	    (display (point-x vertex-1))(display " ")
	    (display (point-y vertex-1))(newline)
	    (display (point-x vertex-2))(display " ")
	    (display (point-y vertex-2))(newline)
	    (display (point-x vertex-3))(display " ")
	    (display (point-y vertex-3))(newline)
	    (display (point-x vertex-1))(display " ")
	    (display (point-y vertex-1))(newline)(newline)
	    (loop (+ counter 1)))
	  (loop (+ counter 1))))))
  (with-output-to-file 
    file 
    (lambda () 
      (write-triangulation 
	(polygon-triangulation polygon)
	(polygon-vertices polygon)))))


;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; \\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///
;; ///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
