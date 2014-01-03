;; This is the file that handles all plane geometrical constructions
;; which will be used in the *Triangulation codes.  It is meant to be
;; understood as a set of legal operations in the classical geometry
;; sense; therefore, we will have basic definitions of the objects we
;; are to handle, basic manipulations of those objects, and complex
;; manipulations using the basic ones as building blocks.

(define-structure point x y name)
;; We assign a name to the points for three reasons: debugging, display
;; and as an aid to translation from this objects to Meroon objects.

(define origin (make-point 0.0 0.0 'o))

(define-structure segment point-1 point-2)

(define no-segment (make-segment origin origin))

(define-structure triangle vertex-1 vertex-2 vertex-3
  			   ptriangle-1-2 ptriangle-2-3 ptriangle-3-1)
;; A triangle is a three points structure with pointers to
;; the triangles adjacent to each of its edges.  Those pointers
;; are meant to be plain integers, indicating the entry of
;; the pointed triangle in a vector (the triangulation vector).

(define no-triangle -1)

(define-structure line base-point direction)
;; We should not define a line using its slope.  Slopes imply division
;; and divisions introduce huge roundoff errors in computations.  It
;; will mess all tests for parallel, intersection, etc. 
;; Instead, we will represent a line by its polar coordinates.
;; line = [ x : y ] and a point through it.

(define x-axis (make-line origin (make-point 1 0 'direction)))
(define y-axis (make-line origin (make-point 0 1 'direction)))

(define-structure rectangle ne se sw nw)

(define-structure super-polygon vertices segments triangulation)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define pi 3.14159265358969)

(define (sign x)
  (cond ((= x 0)  0)
	((> x 0)  1)
	((< x 0) -1)))

(define (square x) (* x x))

;; SAME-POLAR is a division-free way of checking if two directions
;; are parallel.
(define (same-polar-points? polar-point-1 polar-point-2)
  (= (* (point-x polar-point-1) (point-y polar-point-2))
     (* (point-x polar-point-2) (point-y polar-point-1))))

;; SOLVE-SYSTEM solves the following linear system of equations.
;;  a*x + b*y = c
;;  d*x + e*y = f
;; it assumes that the system has solution, since it will be called
;; by procedures that already know that a solution exists.
(define (solve-system a b c
		      d e f)
  (list (/ (- (* b f) (* c e))
	   (- (* b d) (* a e)))
	(/ (- (* c d) (* a f))
	   (- (* b d) (* a e)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  point functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (same-points? point-1 point-2)
  (and (= (point-x point-1) (point-x point-2))
       (= (point-y point-1) (point-y point-2))))

(define (distance point-1 point-2)
  (sqrt (+ (square (- (point-x point-1) (point-x point-2)))
	   (square (- (point-y point-1) (point-y point-2))))))

(define (line-through-points point-1 point-2)
  (make-line point-1 (subtract-points point-2 point-1)))

(define (add-points point-1 point-2)
  (make-point (+ (point-x point-1)
		 (point-x point-2))
	      (+ (point-y point-1)
		 (point-y point-2))
	      'point-addition))

(define (subtract-points point-1 point-2)
  (add-points point-1
	      (point-times-scalar point-2 -1)))

(define (point-times-scalar point scalar)
  (make-point (* scalar (point-x point))
	      (* scalar (point-y point))
	      'point-x-scalar))

(define (angle side-point-1 corner-point side-point-2)
  (let ((query-angle 
	  (- (segment-angle 
	       (make-segment corner-point side-point-2)) 
	     (segment-angle 
	       (make-segment corner-point side-point-1)))))
    (cond ((< query-angle -180)
           (+ 360 query-angle))
	  ((> query-angle 180)
	   (- query-angle 360))
	  (else query-angle))))

(define (collinear-points? point-1 point-2 point-3)
  (parallel-lines? (line-through-points point-1 point-2)
		   (line-through-points point-1 point-3)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  segment and point-segment functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (same-segments? segment-1 segment-2)
  (or (and (same-points? (segment-point-1 segment-1)
	                 (segment-point-1 segment-2))
	   (same-points? (segment-point-2 segment-1)
		         (segment-point-2 segment-2)))
      (and (same-points? (segment-point-1 segment-1)
	                 (segment-point-2 segment-2))
	   (same-points? (segment-point-2 segment-1)
		         (segment-point-1 segment-2)))))

(define (segment-length segment)
  (distance (segment-point-1 segment) (segment-point-2 segment)))

(define (segment-slope segment)
  (let* ((direction (line-direction (segment-line segment)))
	 (run   (point-x direction))
	 (raise (point-y direction)))
    (if (zero? run)
        (* raise +inf.)
	(/ raise run))))

(define (segment-direction segment)
  (line-direction (segment-line segment)))

(define (segment-angle segment)
  (if (zero? (segment-slope segment))
      (if (< (point-x (segment-point-2 segment))
	     (point-x (segment-point-1 segment))) 
	  180 
	  0)
      (let ((radians-angle (atan (segment-slope segment))))
      (if (< (point-x (segment-point-2 segment))
	     (point-x (segment-point-1 segment))) 
	  (+ 180 (* 180 (/ 1 pi) radians-angle))
	  (if (< radians-angle 0)
	      (+ 360 (* 180 (/ 1 pi) radians-angle))
	      (* 180 (/ 1 pi) radians-angle))))))

(define (segment-line segment)
  (line-through-points (segment-point-1 segment)
		       (segment-point-2 segment)))
	     
(define (segment-midpoint segment)
  (let* ((target 
	   (point-times-scalar 
	     (add-points (segment-point-1 segment)
			 (segment-point-2 segment)) 
	     (/ 1 2))))
    (point-name-set! target 'segment-midpoint)
    target))

(define (segment-bisector segment)
  (perpendicular-line-through-point (segment-midpoint segment)
				    (segment-line segment)))

(define (point-encroaches-segment? point segment)
  (< (distance point (segment-midpoint segment))
     (distance (segment-point-1 segment)
	       (segment-midpoint segment))))

(define (point-inside-segment? point segment)
  (and (collinear-points? point 
			  (segment-point-1 segment)
			  (segment-point-2 segment))
       (point-encroaches-segment? point segment)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  line and point-line functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (line-functional line)
  (let ((x0 (point-x (line-base-point line)))
	(y0 (point-y (line-base-point line)))
	(dx (point-x (line-direction line)))
	(dy (point-y (line-direction line))))
    (lambda (x y)
      (- (* dx (- y y0)) (* dy (- x x0))))))

(define (perpendicular-line-through-point point line)
  (let* ((direction (line-direction line))
	 (x (point-x direction))
	 (y (point-y direction)))
    (make-line point (make-point (* -1 y) x 'polar))))

(define (parallel-line-through-point point line)
  (make-line point (line-direction line)))

(define (line-intersection line-1 line-2 name)
  (if (same-polar-points? (line-direction line-1)
			  (line-direction line-2))
    (error "Parallel (or coincident) lines")
    (let* ((x1 (point-x (line-base-point line-1)))
	   (y1 (point-y (line-base-point line-1)))
	   (x2 (point-x (line-base-point line-2)))
	   (y2 (point-y (line-base-point line-2)))
	   (dx1 (point-x (line-direction line-1)))
	   (dy1 (point-y (line-direction line-1)))
	   (dx2 (point-x (line-direction line-2)))
	   (dy2 (point-y (line-direction line-2)))
	   (it1 (- (* dy1 x1) (* dx1 y1)))
	   (it2 (- (* dy2 x2) (* dx2 y2)))
	   (xy  (solve-system dy1 (- dx1) it1
			      dy2 (- dx2) it2))
	   (x   (car xy))
	   (y   (cadr xy)))
      (make-point x y name))))

(define (parallel-lines? line-1 line-2)
  (same-polar-points? (line-direction line-1) 
		      (line-direction line-2)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  triangle and point-triangle functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (triangle-circumcenter triangle)
  (let* ((segment-1 (make-segment (triangle-vertex-1 triangle)
				  (triangle-vertex-2 triangle)))
	 (segment-2 (make-segment (triangle-vertex-2 triangle)
				  (triangle-vertex-3 triangle)))
	 (segment-3 (make-segment (triangle-vertex-3 triangle)
				  (triangle-vertex-1 triangle))))
    (cond ((not (same-polar-points? (segment-direction segment-1)
				    (segment-direction segment-2)))
	   (line-intersection 
	     (perpendicular-line-through-point
	       (segment-midpoint segment-1)
	       (segment-line segment-1))
	     (perpendicular-line-through-point
	       (segment-midpoint segment-2)
	       (segment-line segment-2))
	     'circumcenter))
	  ((not (same-polar-points? (segment-direction segment-1)
				    (segment-direction segment-3)))
	   (line-intersection 
	     (perpendicular-line-through-point
	       (segment-midpoint segment-1)
	       (segment-line segment-1)) 
	     (perpendicular-line-through-point
	       (segment-midpoint segment-3)
	       (segment-line segment-3)) 
	     'circumcenter))
	  ((not (same-polar-points? (segment-direction segment-2)
				    (segment-direction segment-3)))
	   (line-intersection 
	     (perpendicular-line-through-point
	       (segment-midpoint segment-2)
	       (segment-line segment-2))
	     (perpendicular-line-through-point 
	       (segment-midpoint segment-3)
	       (segment-line segment-3))
	     'circumcenter)))))
;	  (else (display-triangle triangle)
;		(display-point (triangle-vertex-1 triangle))
;		(display-point (triangle-vertex-2 triangle))
;		(display-point (triangle-vertex-3 triangle))
;		(error "fake triangle!")))))

(define (triangle-incenter triangle)
  (point-times-scalar
    (add-points (triangle-vertex-1 triangle)
		(add-points (triangle-vertex-2 triangle)
			    (triangle-vertex-3 triangle)))
    (/ 1 3)))

(define (point-inside-circumcircle? point triangle)
  (if (point? (triangle-circumcenter triangle))
    (> (distance (triangle-circumcenter triangle) (triangle-vertex-1 triangle))
       (distance (triangle-circumcenter triangle) point))
    (= 1 0)))

(define (point-inside-triangle? point triangle)
  (let* ((segment-1 (make-segment (triangle-vertex-1 triangle)
				  (triangle-vertex-2 triangle)))
	 (segment-2 (make-segment (triangle-vertex-2 triangle)
				  (triangle-vertex-3 triangle)))
	 (segment-3 (make-segment (triangle-vertex-3 triangle)
				  (triangle-vertex-1 triangle)))
	 (function1 (line-functional (segment-line segment-1)))
	 (function2 (line-functional (segment-line segment-2)))
	 (function3 (line-functional (segment-line segment-3))))
    (= (sign (function1 (point-x point) (point-y point)) 0)
       (sign (function2 (point-x point) (point-y point)) 0)
       (sign (function3 (point-x point) (point-y point)) 0))))

(define (point-inside-triangle? point triangle)
  (let* ((vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (polygon  (list->vector (list vertex-1 vertex-2 vertex-3))))
    (point-inside-polygon? point polygon)))

(define (triangle-minimum-angle triangle)
  (min (abs (angle (triangle-vertex-1 triangle)
		   (triangle-vertex-2 triangle)
		   (triangle-vertex-3 triangle)))
       (abs (angle (triangle-vertex-2 triangle)
		   (triangle-vertex-3 triangle)
		   (triangle-vertex-1 triangle)))
       (abs (angle (triangle-vertex-3 triangle)
		   (triangle-vertex-1 triangle)
		   (triangle-vertex-2 triangle)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  segment-triangle and line-segment functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (is-segment-a-triangle-edge? segment triangle)
  (let ((vertex-1 (triangle-vertex-1 triangle))
	(vertex-2 (triangle-vertex-2 triangle))
	(vertex-3 (triangle-vertex-3 triangle)))
      (or (same-segments? segment (make-segment vertex-1 vertex-2))
	  (same-segments? segment (make-segment vertex-2 vertex-3))
	  (same-segments? segment (make-segment vertex-3 vertex-1)))))

(define (line-intersects-segment? line segment)
  (let ((cross (line-intersection line (segment-line segment) 'i)))
    (point-inside-segment? cross segment)))

(define (opposite-vertex-to-segment segment triangle)
  (let ((segment-1 (make-segment (triangle-vertex-1 triangle)
				 (triangle-vertex-2 triangle)))
	(segment-2 (make-segment (triangle-vertex-2 triangle)
				 (triangle-vertex-3 triangle)))
	(segment-3 (make-segment (triangle-vertex-3 triangle)
				 (triangle-vertex-1 triangle))))
    (cond ((same-segments? segment-1 segment)
	   (triangle-vertex-3 triangle))
	  ((same-segments? segment-2 segment)
	   (triangle-vertex-1 triangle))
	  ((same-segments? segment-3 segment)
	   (triangle-vertex-2 triangle))
	  (else 
	    (display-triangle triangle)
	    (display-segment segment)
	    (error "input segment not edge of triangle")))))

(define (pointer-associated-to-segment segment triangle)
  (let ((segment-1 (make-segment (triangle-vertex-1 triangle)
				 (triangle-vertex-2 triangle)))
	(segment-2 (make-segment (triangle-vertex-2 triangle)
				 (triangle-vertex-3 triangle)))
	(segment-3 (make-segment (triangle-vertex-3 triangle)
				 (triangle-vertex-1 triangle))))
    (cond ((same-segments? segment-1 segment)
	   (triangle-ptriangle-1-2 triangle))
	  ((same-segments? segment-2 segment)
	   (triangle-ptriangle-2-3 triangle))
	  ((same-segments? segment-3 segment)
	   (triangle-ptriangle-3-1 triangle))
	  (else (error "input segment not edge of triangle")))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  segment-triangulation functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (is-segment-in-triangulation? segment triangulation)
  (let loop ((position 0) (result #f))
    (if (= position (vector-length triangulation))
      result
      (loop (+ position 1) 
	    (or 
	      result 
	      (is-segment-a-triangle-edge? segment 
				       (vector-ref triangulation 
						   position)))))))

(define (is-segment-a-boundary-edge? segment segments)
  (let loop ((position 0) (result #f))
    (if (= position (vector-length segments))
      result
      (loop (+ position 1) 
	    (or result 
		(same-segments? segment 
				(vector-ref segments position)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; polygon functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; WARNING!! Polygon is not just a set of vertices.  We define a polygon
;; as a vector of points that delimite the polygonal region on its left
;; hand side... Messy? an example:
;; Let p0=(0,0), p1=(2,0) and p2=(1,1).
;; (list->vector (list p0 p1 p2)) is a polygon, but
;; (list->vector (list p0 p2 p1)) is not.
;; Of course, in the definition is implicit the fact that no interior
;; points are allowed, anly vertices joined by boundary edges.

(define (segments->polygon segments)
  (let loop ((position 0) (result '()))
    (if (= position (vector-length segments))
      (list->vector result)
      (loop (+ position 1)
	    (append result
		    (list (segment-point-1
			    (vector-ref segments position))))))))

(define (polygon->segments polygon)
  (let ((first-vertex (vector-ref polygon 0)))
    (let loop ((position 1) (result '()))
      (if (= position (vector-length polygon))
	(list->vector
	  (append result
		  (list 
		    (make-segment (vector-ref polygon (- position 1))
				  first-vertex))))
	(loop (+ position 1)
	      (append result
		      (list
			(make-segment
			  (vector-ref polygon (- position 1))
			  (vector-ref polygon position)))))))))

;; BOUNDARY-BOX creates a rectangle with dimensions three times that
;; of the smallest horizontal rectangle containing the given polygon.
;; This assures us that the circumcenters of the triangles are within
;; reach.
(define (boundary-box polygon)
  (let* ((limits (polygon-limits polygon))
	 (xmin (vector-ref limits 0))
	 (xmax (vector-ref limits 1))
	 (ymin (vector-ref limits 2))
	 (ymax (vector-ref limits 3))) 
    (make-rectangle 
      (make-point (- (* 2 xmax) xmin) (- (* 2 ymax) ymin) 'NE) 
      (make-point (- (* 2 xmax) xmin) (- (* 2 ymin) ymax) 'SE) 
      (make-point (- (* 2 xmin) xmax) (- (* 2 ymin) ymax) 'SW) 
      (make-point (- (* 2 xmin) xmax) (- (* 2 ymax) ymin) 'NW))))

(define (polygon-limits polygon)
  (let loop ((position 0) (xmin 0) (xmax 0) (ymin 0) (ymax 0))
    (if (= position (vector-length polygon))
      (list->vector (list xmin xmax ymin ymax))
      (let ((vertex (vector-ref polygon position))) 
	(loop (+ position 1) 
	      (min xmin (point-x vertex))
	      (max xmax (point-x vertex))
	      (min ymin (point-y vertex))
	      (max ymax (point-y vertex)))))))

(define (point-inside-polygon? point polygon)
  (let* ((length (vector-length polygon))
	 (last-vertex (vector-ref polygon (- length 1)))
	 (first-vertex (vector-ref polygon 0))) 
    (let loop ((position 1) 
	       (result (angle last-vertex 
			      point
			      first-vertex))) 
      (if (= position length)
	(> (abs result) 359)
	(loop (+ position 1) 
	      (+ result 
		 (angle 
		   (vector-ref polygon (- position 1)) 
		   point
		   (vector-ref polygon position))))))))

(define (convex-polygon-corner? entry polygon)
  (cond ((= entry 0)
	 (point-inside-polygon? 
	   (triangle-incenter 
	     (make-triangle (vector-ref polygon 0)
			    (vector-ref polygon 1)
			    (vector-ref polygon 
					(- (vector-length polygon) 
					   1))
			    -1 -1 -1))
	   polygon))
	((= entry (- (vector-length polygon) 1))
	 (point-inside-polygon?
	   (triangle-incenter
	     (make-triangle (vector-ref polygon entry)
			    (vector-ref polygon (- entry 1))
			    (vector-ref polygon 0)
			    -1 -1 -1))
	   polygon))
	(else
	  (point-inside-polygon?
	    (triangle-incenter
	      (make-triangle (vector-ref polygon entry)
			     (vector-ref polygon (- entry 1))
			     (vector-ref polygon (+ entry 1))
			     -1 -1 -1))
	      polygon))))


(define (is-inner-triangle? triangle polygon)
  (point-inside-polygon? (triangle-incenter triangle) polygon))

(define (polygon-point-entry point polygon)
  (let loop ((position 0))
    (if (= position (vector-length polygon))
      (error "point is not vertex of polygon")
      (if (same-points? point (vector-ref polygon position))
	position
	(loop (+ position 1))))))

(define (is-corner-triangle? triangle segments)
  ;;
  ;; At least two of the edges of this triangles must be boundary-edges
  ;; 
  (let* ((segment-1 (make-segment (triangle-vertex-1 triangle) 
				  (triangle-vertex-2 triangle)))
	 (segment-2 (make-segment (triangle-vertex-2 triangle)
				  (triangle-vertex-3 triangle)))
	 (segment-3 (make-segment (triangle-vertex-3 triangle)
				  (triangle-vertex-1 triangle)))
	 (premise-1 (is-segment-a-boundary-edge?
		      segment-1 segments))
	 (premise-2 (is-segment-a-boundary-edge?
		      segment-2 segments))
	 (premise-3 (is-segment-a-boundary-edge?
		      segment-3 segments)))
    (or (and premise-1 premise-2)
	(and premise-1 premise-3)
	(and premise-2 premise-3))))
    
(define (triangle-corner-angle-in-polygon corner-triangle polygon)
  (let ((top (- (vector-length polygon) 1))
	(entry-1 (polygon-point-entry
		   (triangle-vertex-1 corner-triangle) polygon))
	(entry-2 (polygon-point-entry
		   (triangle-vertex-2 corner-triangle) polygon))
	(entry-3 (polygon-point-entry
		   (triangle-vertex-3 corner-triangle) polygon)))
    (cond ((or (= 2   (- (max entry-1 entry-3) 
	                 (min entry-1 entry-3)))
	       (= top (- (max entry-1 entry-3)
		         (min entry-1 entry-3))))
	   (abs (angle (triangle-vertex-1 corner-triangle)
		       (triangle-vertex-2 corner-triangle)
		       (triangle-vertex-3 corner-triangle)))) 
	  ((or (= 2   (- (max entry-1 entry-2) 
	                 (min entry-1 entry-2)))
	       (= top (- (max entry-1 entry-2)
		         (min entry-1 entry-2))))
	   (abs (angle (triangle-vertex-1 corner-triangle)
		       (triangle-vertex-3 corner-triangle)
		       (triangle-vertex-2 corner-triangle)))) 
	  ((or (= 2   (- (max entry-2 entry-3) 
	                 (min entry-2 entry-3)))
	       (= top (- (max entry-2 entry-3)
		         (min entry-2 entry-3))))
	   (abs (angle (triangle-vertex-2 corner-triangle)
		       (triangle-vertex-1 corner-triangle)
		       (triangle-vertex-3 corner-triangle)))))))
