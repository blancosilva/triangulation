;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
;;;;;;   basic geometric constructions
;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (midpoint point-1 point-2)
  (make-point (/ (+ (point-x point-1)
		    (point-x point-2)) 
		 2.0)
	      (/ (+ (point-y point-1)
		    (point-y point-2)) 
		 2.0)
	      'SS))

(define (circumcenter triangle)
  (let* ((x1 (point-x (triangle-vertex-1 triangle)))
	 (y1 (point-y (triangle-vertex-1 triangle)))
	 (x2 (point-x (triangle-vertex-2 triangle)))
	 (y2 (point-y (triangle-vertex-2 triangle)))
	 (x3 (point-x (triangle-vertex-3 triangle)))
	 (y3 (point-y (triangle-vertex-3 triangle)))
	 (center-x (/ (- (+ (* x1 x1 y2)
			    (* y1 y3 y3)
			    (* y1 x3 x3)
			    (* y1 y1 y2)
			    (* y3 y2 y2)
			    (* y3 x2 x2))
			 (+ (* x1 x1 y3)
			    (* x3 x3 y2)
			    (* y1 y2 y2)
			    (* y1 x2 x2)
			    (* y3 y3 y2)
			    (* y1 y1 y3)))
		      (- (+ (* x2 y3) 
			    (* x1 y2)
			    (* x3 y1))
			 (+ (* x1 y3)
			    (* x2 y1)
			    (* x3 y2)))
		      2.0))
	 (center-y (/ (- (+ (* x1 x3 x3) 
			    (* x1 y3 y3) 
			    (* x1 x1 x2) 
			    (* x2 y1 y1) 
			    (* x2 x2 x3) 
			    (* y2 y2 x3)) 
			 (+ (* x1 x1 x3) 
			    (* x2 x3 x3) 
			    (* x2 y3 y3) 
			    (* x2 x2 x1) 
			    (* y1 y1 x3) 
			    (* y2 y2 x1))) 
		      (- (+ (* x2 y3) 
			    (* x1 y2) 
			    (* x3 y1)) 
			 (+ (* x1 y3) 
			    (* x2 y1) 
			    (* x3 y2))) 
		      -2.0)))
    (make-point center-x center-y 'C)))

(define (angle point-1 point-2 point-3)
  (let* ((x1 (point-x point-1))
	 (y1 (point-y point-1))
	 (x2 (point-x point-2))
	 (y2 (point-y point-2))
	 (x3 (point-x point-3))
	 (y3 (point-y point-3))
	 (inner-product (+ (* (- x1 x2) (- x3 x2)) (* (- y1 y2) (- y3 y2))))
	 (dist2a (+ (* (- x1 x2) (- x1 x2)) (* (- y1 y2) (- y1 y2))))
	 (dist2b (+ (* (- x3 x2) (- x3 x2)) (* (- y3 y2) (- y3 y2))))
	 (cosine (/ inner-product (sqrt (* dist2a dist2b))))
	 (discr (discriminant x1 y1 x2 y2 x3 y3)))
	(/ (* 180 (sign discr) (acos cosine)) pi)))
    
(define (square-distance point-1 point-2)
  (let ((x1 (point-x point-1))
	(y1 (point-y point-1))
	(x2 (point-x point-2))
	(y2 (point-y point-2)))
    (+ (* (- x2 x1) (- x2 x1))
       (* (- y2 y1) (- y2 y1)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
;;;;;;   Basic list-vector manipulations
;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (insert-in-list element after-element items-list)
  (let* ((items-vector (list->vector items-list))
	 (augmented-vector (make-vector (+ 1 (vector-length items-vector)))))
    (let loop ((position 0) (pointer 0))
      (if (= (vector-length items-vector) position)
	  (vector->list augmented-vector)
	  (if (equal? (vector-ref items-vector position) after-element)
	      (begin
		(vector-set! augmented-vector position after-element)
		(vector-set! augmented-vector (+ position 1) element)
		(loop (+ position 1) (+ pointer 2)))
	      (begin
		(vector-set! augmented-vector pointer (vector-ref items-vector position))
		(loop (+ position 1) (+ pointer 1))))))))

(define (locate-entry-in-vector equality-function entry vector)
  (let loop ((position 0))
    (if (= (vector-length vector) position)
	(error "unable to locate entry in vector")
	(if (equality-function entry (vector-ref vector position))
	    position
	    (loop (+ position 1))))))

(define (replace-entry-in-vector equality-function old-entry new-entry vector)
  (vector-set! vector (locate-entry-in-vector equality-function old-entry vector) new-entry)
  vector)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;
;;;;;;;;;;;  Abstractions abstractions...
;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (get-segments polygon)  ;;;; be careful with this one. Use it only ONCE
  (let ((poly-vector (list->vector polygon))
	(segments-vector (make-vector (length polygon))))
    (let loop ((position 1))
      (if (= (length polygon) position)
	  (begin
	    (vector-set! segments-vector (- position 1) (list (vector-ref poly-vector (- position 1))
							      (vector-ref poly-vector 0)))
	    segments-vector)
	  (begin
	    (vector-set! segments-vector 
			 (- position 1) 
			 (list (vector-ref poly-vector (- position 1))
			       (vector-ref poly-vector position)))
	    (loop (+ position 1)))))))

(define-structure super-polygon vertices segments triangulation)
;; A polygon was a list of vertices.  A super-polygon carries some extra vertices,
;; and the set of segments is either the set of original segments of the original
;; polygon, or maybe each segment on that one is substituted by some pieces.  We
;; also include a triangulation of the set of vertices, but such that every segment
;; must be an edge of the triangulation.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
;;;;;;  Introduction of segments and basic utilities
;;;;;;  (a segment is a pair of unordered vertices)
;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (same-segments? segment-1 segment-2)
  (let ((vertex-1-1 (car segment-1))
	(vertex-2-1 (cadr segment-1))
	(vertex-1-2 (car segment-2))
	(vertex-2-2 (cadr segment-2)))
    (or (and (equal? vertex-1-1 vertex-1-2)
	     (equal? vertex-2-1 vertex-2-2))
	(and (equal? vertex-1-1 vertex-2-2)
	     (equal? vertex-2-1 vertex-1-2)))))

(define (is-segment-in-triangulation? segment triangulation)
  (define (is-segment-in-triangle? triangle)
    (let ((vertex-1 (triangle-vertex-1 triangle))
	  (vertex-2 (triangle-vertex-2 triangle))
	  (vertex-3 (triangle-vertex-3 triangle)))
      (or (same-segments? segment (list vertex-1 vertex-2))
	  (same-segments? segment (list vertex-2 vertex-3))
	  (same-segments? segment (list vertex-3 vertex-1)))))
  (let loop ((position 0))
    (if (= (vector-length triangulation) position)
	(= 0 1) ;;; xemacs messes parenthesis matching if I use #f or #t
	(if (is-segment-in-triangle? (vector-ref triangulation position))
	    (= 0 0)  ;; same story
	    (loop (+ position 1))))))

(define (include-segments-in-triangulation super-polygon)
  (let loop ((position 0) (vertices (super-polygon-vertices super-polygon))
	     (segments (super-polygon-segments super-polygon))
	     (triangulation (super-polygon-triangulation super-polygon)))
    (if (= (vector-length segments) position)
	(make-super-polygon vertices segments triangulation)
	(let* ((segment (vector-ref segments position))
	       (vertex-1 (car segment))
	       (vertex-2 (cadr segment))
	       (middle (midpoint vertex-1 vertex-2)))
	  (if (is-segment-in-triangulation? segment triangulation)
	      (loop (+ position 1) vertices segments triangulation)
	      (begin
		(vector-set! segments position (list vertex-1 middle))
		(let* ((new-vertices (insert-in-list middle vertex-1 vertices))
		       (new-segments (list->vector (insert-in-list (list middle vertex-2)
								   (list vertex-1 middle)
								   (vector->list segments))))
		       (new-triangulation (split-something middle triangulation)))
		  (loop 0 new-vertices new-segments new-triangulation))))))))

(define (first-segment-encroached-upon-by-triangle-circumcenter triangle segments)
  (let loop ((position 0))
    (if (= position (vector-length segments))
	(list origen origen)
	(let* ((segment (vector-ref segments position))
	       (point (circumcenter triangle))
	       (vertex-1 (car segment))
	       (vertex-2 (cadr segment))
	       (middle (midpoint vertex-1 vertex-2)))
	  (if (< (square-distance point middle) (square-distance vertex-1 middle))
	      segment
	      (loop (+ position 1)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;
;;;;;;;;;;;;;  Introduction of angle handling (notice that we no longer have polygons 
;;;;;;;;;;;;;  instead we have "vertices", since we do not consider the edges connec-
;;;;;;;;;;;;;  ting the interior vertices we introduce as forced segments in our tri-
;;;;;;;;;;;;;  angulation)
;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (get-rid-of-angles-less-than-threshold threshold-angle super-polygon xmin xmax ymin ymax)
  (define (helper position super-polygon)
    (define (encroaching-step triangle)
      (let* ((vertices (super-polygon-vertices super-polygon))
	     (segments (super-polygon-segments super-polygon))
	     (triangulation (super-polygon-triangulation super-polygon))
	     (bad-segment (first-segment-encroached-upon-by-triangle-circumcenter triangle segments)))
	(if (equal? bad-segment no-segment)
	    (let ((target-point (circumcenter triangle)))
	      (if (and (< xmin (point-x target-point))
		       (> xmax (point-x target-point))
		       (< ymin (point-y target-point))
		       (> ymax (point-y target-point)))
		  (begin
					;(display "\nIncluding in structure triangle-circumcenter. Triangle:\n")
					;(display-triangle triangle)
		    (let* ((new-triangulation (split-something target-point triangulation))
			   (new-vertices (append vertices (list target-point)))
			   (new-super-polygon (make-super-polygon new-vertices segments new-triangulation))
			   (newest-super-polygon (include-segments-in-triangulation new-super-polygon)))
		      (helper 0 newest-super-polygon)))
		  (helper (+ position 1) super-polygon)))
	    (let* ((vertex-1 (car bad-segment))
		   (vertex-2 (cadr bad-segment))
		   (middle (midpoint vertex-1 vertex-2)))
	      (replace-entry-in-vector same-segments? bad-segment (list vertex-1 middle) segments)
	      (let* ((new-segments (list->vector (insert-in-list (list middle vertex-2)
								 (list vertex-1 middle)
								 (vector->list segments))))
		     (new-vertices (insert-in-list middle vertex-1 vertices))
		     (new-triangulation (split-something middle triangulation))
		     (new-super-polygon (make-super-polygon new-vertices new-segments new-triangulation))
		     (newest-super-polygon (include-segments-in-triangulation new-super-polygon)))
		(helper 0 newest-super-polygon))))))
    (if (= position (vector-length (super-polygon-triangulation super-polygon)))
	super-polygon
	(let* ((segments (super-polygon-segments super-polygon))
	       (triangulation (super-polygon-triangulation super-polygon))
	       (triangle (vector-ref triangulation position)))
	  (cond ((or (>= (minimum-angle triangle) threshold-angle)
		     (not (is-inner-triangle? triangle segments)))
		 (helper (+ position 1) super-polygon))
		((is-corner-triangle? triangle segments)
		 (if (= (minimum-angle triangle) (corner-angle triangle segments))
		     (helper (+ position 1) super-polygon)
		     (encroaching-step triangle)))
		(else (encroaching-step triangle))))))
  (helper 0 super-polygon))

(define (delaunay-triangulation threshold-angle polygon)
  (let* ((box (polygon-boundary polygon))
	 (N (vector-ref box 0))
	 (E (vector-ref box 1))
	 (S (vector-ref box 2))
	 (W (vector-ref box 3))
	 (xmin (point-x W))
	 (xmax (point-x N))
	 (ymin (point-y S))
	 (ymax (point-y N))
	 (initial-triangulation (polygon-triangulation polygon))
	 (initial-segments      (get-segments polygon))
	 (initial-super-polygon (make-super-polygon polygon initial-segments initial-triangulation))
	 (second-super-polygon  (include-segments-in-triangulation initial-super-polygon)))
    (get-rid-of-angles-less-than-threshold threshold-angle second-super-polygon xmin xmax ymin ymax)))
