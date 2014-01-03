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

;; A point is a complex number:
(define point-x real-part)
(define point-y imag-part)
(define (make-point x y) (+ x (* +i y)))
;; An edge is defined by two points in no particular order, pointers to the
;; adjacent triangles (in numerical order), and two extra fields: 
;;     boundary?: Set to #t if the edge is (part of) an original segment of the
;;     PSGL.  
;;     useful?:   Set to #t if it is needed in the triangulation.
(define-structure edge from to 
		  triangle-pointer-1 triangle-pointer-2 
		  boundary? useful?)
;; A triangle is defined by its three vertices, and pointers to the adjacent
;; edges and triangles.
(define-structure triangle 
		  vertex-1 vertex-2 vertex-3
		  ;; triangle-k and segment-k are opposite to vertex-k
		  ;; a pointer to a no-triangle has value -1
		  edge-pointer-1 edge-pointer-2 edge-pointer-3
		  triangle-pointer-1 triangle-pointer-2 triangle-pointer-3)

;; Collections:  
(include "Dvector.scm")
(define (collection . elements)
  (let ((set (build-Dvector 0)))
    (do ((counter (- (length elements) 1) (- counter 1)))
      ((< counter 0) set)
      (Dvector-set! set counter (list-ref elements counter)))))
(define collection-length Dvector-length)
(define collection-ref Dvector-ref)
(define collection-set! Dvector-set!)
(define (add-element-to-collection set element)
  (Dvector-append set element))
(define (add-elements-to-collection set . elements)
  (do ((counter (- (length elements) 1) (- counter 1)))
    ((< counter 0) void)
    (add-element-to-collection set (list-ref elements counter))))
(define (replace-element-in-collection
	  set old-element-pointer . new-elements)
  (let ((size (length new-elements)))
    (do ((counter 1 (+ counter 1)))
      ((= counter size)
       (collection-set! set old-element-pointer (car new-elements)))
      (add-element-to-collection set (list-ref new-elements counter)))))

;; A parser for .poly files.
(define (read-file file-name)

  (define (read-points port points-size)
    (let ((result (collection)))
      (do ((counter 0 (+ counter 1)))
	((= counter points-size) result)
	(let ((index (- (read port) 1))
	      (x-coordinate (read port))
	      (y-coordinate (read port)))
	  (collection-set! 
	    result index (make-point x-coordinate y-coordinate))))))

  (define (read-edges port edges-size points)
    (let ((result (collection)))
      (do ((counter 0 (+ counter 1)))
	((= counter edges-size) result)
	(let ((index (- (read port) 1))
	      (from-index (- (read port) 1))
	      (to-index (- (read port) 1)))
	  (collection-set! result index
		       (make-edge (collection-ref points from-index)
				     (collection-ref points to-index)
				     -1 -1 #t #t))))))

  (let* ((port (open-input-file file-name))
	 (vertices-size (read port))
	 (ignore (read port))
	 (ignore (read port))
	 (ignore (read port))
	 (vertices (read-points port vertices-size))
	 (edges-size (read port))
	 (ignore (read port))
	 (segments (read-edges port edges-size vertices))
	 (holes-size (read port))
	 (holes (read-points port holes-size)))
    (list vertices segments holes)))

#|   XXXXX             XX                                           
      X   X             X                                           
      X    X            X                                           
      X    X  XXXXX     X     XXXX   XX  XX  XX XX    XXXX   XXX XXX
      X    X X     X    X         X   X   X   XX  X       X   X   X 
      X    X XXXXXXX    X     XXXXX   X   X   X   X   XXXXX   X   X 
      X    X X          X    X    X   X   X   X   X  X    X    X X  
      X   X  X     X    X    X    X   X  XX   X   X  X    X    X X  
     XXXXX    XXXXX   XXXXX   XXXX X   XX XX XXX XXX  XXXX X    X   
								X   
							      XX    |#

(define (delaunay-triangulation parsed-poly-file)

  ;; First step: Construction of a naive triangulation of a region big
  ;; enough to contain all the vertices of the given poly file.

  ;; Second step: Introduction of the vertices of the poly file in the previous
  ;; triangulation; one at a time.  Each vertex will be lying in either a
  ;; triangle or an edge of the current triangulation.  That object must be
  ;; split, creating new triangles and edges.  The old triangle(s)/edge(s)
  ;; supportng the new point must be erased from the triangulation.  And
  ;; finally, a check must be performed to ensure all edges are legal in the
  ;; sense of Delaunay (no vertex lies in the interior of the circumcircle of
  ;; any other three vertices).
  (let* ((vertices (car parsed-poly-file))
	 (holes (caddr parsed-poly-file))
	 (initial-setup (initial-setup vertices holes))
	 (edges (cdr initial-setup))
	 (triangulation (car initial-setup)))
    (do ((index (- (collection-length vertices) 1) (- index 1)))
      ((< index 0) 
       (cons triangulation edges))
      (begin
	(insert-point-in-triangulation
	  (collection-ref vertices index) edges triangulation)))))

#|                          X                           
           X                                            
           X                                            
  XXXXX   XXXX    XXXXX   XXX    XX XX    XXXXX  XXX XX 
 X     X   X     X     X    X     XX  X  X     X   XX  X
  XXX      X     XXXXXXX    X     X   X  XXXXXXX   X    
     XX    X     X          X     X   X  X         X    
 X     X   X  X  X     X    X     X   X  X     X   X    
  XXXXX     XX    XXXXX   XXXXX  XXX XXX  XXXXX  XXXXX  |#
                                                        
(define (constrained-steiner-delaunay-triangulation parsed-poly-file)
  ;; Repeat until each segment is present in the triangulation as the edge of
  ;; at least one triangle: If a given segment is not an edge in the
  ;; triangulation, split the segment in two, replace the old segment with
  ;; these two pieces, and include the midpoint in the triangulation. Start
  ;; all over with the new edges, triangles and segments.
  ;;
  ;; As a side-effect, we update the "boundary?" field of each edge when
  ;; needed. 
  (declare (flonum))
  (let* ((vertices (car parsed-poly-file))
	 (segments (cadr parsed-poly-file))
	 (delaunay-triangulation 
	   (time (delaunay-triangulation parsed-poly-file)))
	 (triangulation (car delaunay-triangulation))
	 (edges (cdr delaunay-triangulation)))
    (let loop ((counter (- (collection-length segments) 1)))
      (if (< counter 0) 
	(cons triangulation edges)
	(let ((segment (collection-ref segments counter)))
	  (if (not (is-segment-in-edges? segment edges))
	    (let* ((vertex-1 (edge-from segment))
		   (vertex-2 (edge-to segment))
		   (midpoint (midpoint vertex-1 vertex-2))
		   (segment-1 (make-edge vertex-1 midpoint -1 -1 #t #t))
		   (segment-2 (make-edge midpoint vertex-2 -1 -1 #t #t)))
	      (replace-element-in-collection
		segments counter segment-1 segment-2)
	      (add-element-to-collection vertices midpoint)
	      (insert-point-in-triangulation midpoint edges triangulation)
	      (loop (- (collection-length segments) 1)))
	    (loop (- counter 1))))))))

#|				       XX           
					X           
					X           
	      XXXX   XX XX    XXXXXX    X     XXXXX 
		  X   XX  X  X    X     X    X     X
	      XXXXX   X   X  X    X     X    XXXXXXX
	     X    X   X   X  X    X     X    X      
	     X    X   X   X   XXXXX     X    X     X
	      XXXX X XXX XXX      X   XXXXX   XXXXX 
				  X                 
			      XXXX                  |#

(define (minimum-angle-constrained-steiner-delaunay-triangulation
	  threshold-angle parsed-poly-file)
  (let* ((vertices (car parsed-poly-file))
	 (segments (cadr parsed-poly-file))
	 (holes (caddr parsed-poly-file))
	 (constrained-steiner-delaunay-triangulation
	   (time (constrained-steiner-delaunay-triangulation parsed-poly-file)))
	 (triangulation (car constrained-steiner-delaunay-triangulation))
	 (edges (cdr constrained-steiner-delaunay-triangulation)))
    (remove-exterior-triangles-from-triangulation holes edges triangulation)
;;     (modify-triangles-with-minimum-angle-greater-than-threshold
;;       threshold-angle segments edges triangulation)
    (cons triangulation edges)))



;; Auxiliary functions for DELAUNAY-TRIANGULATION.

(define (initial-setup vertices holes)

  (define (box-limits points)
    (let ((min-x +inf.) (max-x -inf.) (min-y +inf.) (max-y -inf.))
      (let loop ((counter (- (collection-length points) 1))
		 (min-x +inf.) (max-x -inf.)
		 (min-y +inf.) (max-y -inf.))
	(if (< counter 0) 
	  (list min-x max-x min-y max-y)
	  (let ((point (collection-ref points counter)))
	    (loop (- counter 1)
		  (min min-x (point-x point))
		  (max max-x (point-x point))
		  (min min-y (point-y point))
		  (max max-y (point-y point))))))))

  (let* ((limits (box-limits vertices))
	 (min-x (car limits))
	 (max-x (cadr limits))
	 (min-y (caddr limits))
	 (max-y (cadddr limits)))
    (let* ((NE (make-point (- (* 2 max-x) min-x) (- (* 2 max-y) min-y)))
	   (SE (make-point (- (* 2 max-x) min-x) (- (* 2 min-y) max-y)))
	   (SW (make-point (- (* 2 min-x) max-x) (- (* 2 min-y) max-y)))
	   (NW (make-point (- (* 2 min-x) max-x) (- (* 2 max-y) min-y)))
	   (north (make-edge NW NE -1 1 #f #f))    ;; edge 0
	   (east  (make-edge NE SE -1 0 #f #f))    ;; edge 1
	   (south (make-edge SE SW -1 0 #f #f))    ;; edge 2
	   (west  (make-edge SW NW -1 1 #f #f))    ;; edge 3
	   (diagonal (make-edge NE SW 0 1 #f #f))) ;; edge 4
      (let ((triangle-1 (make-triangle SW SE NE 1 4 2 -1 1 -1))
	    (triangle-2 (make-triangle SW NW NE 0 4 3 -1 0 -1)))
	(add-elements-to-collection
	  holes
	  (make-point (- (point-x NE) 0.001) (- (point-y NE) 0.001)))
	(cons 
	  (collection triangle-1 triangle-2)
	  (collection north east south west diagonal))))))

(define (change-triangle-pointers triangle old-pointer new-pointer)
  (cond ((= old-pointer (triangle-triangle-pointer-1 triangle))
	 (triangle-triangle-pointer-1-set! triangle new-pointer))
	((= old-pointer (triangle-triangle-pointer-2 triangle))
	 (triangle-triangle-pointer-2-set! triangle new-pointer))
	((= old-pointer (triangle-triangle-pointer-3 triangle))
	 (triangle-triangle-pointer-3-set! triangle new-pointer))))

(define (change-edge-pointers edge old-pointer new-pointer)
  (cond ((= old-pointer (edge-triangle-pointer-1 edge))
	 (let ((other-pointer (edge-triangle-pointer-2 edge)))
	   (edge-triangle-pointer-1-set! 
	     edge (min other-pointer new-pointer))
	   (edge-triangle-pointer-2-set! 
	     edge (max other-pointer new-pointer))))
	((= old-pointer (edge-triangle-pointer-2 edge))
	 (let ((other-pointer (edge-triangle-pointer-1 edge)))
	   (edge-triangle-pointer-1-set! 
	     edge (min other-pointer new-pointer))
	   (edge-triangle-pointer-2-set! 
	     edge (max other-pointer new-pointer))))
	(else (error "bad pointer"))))
	     
(define (the-other-vertex vertex-1 vertex-2 triangle)
  (let* ((sorted-vertices (lex-sort 
			    (triangle-vertex-1 triangle)
			    (triangle-vertex-2 triangle)
			    (triangle-vertex-3 triangle)))
	 (triangle-vertex-1 (car sorted-vertices))
	 (triangle-vertex-2 (cadr sorted-vertices))
	 (triangle-vertex-3 (caddr sorted-vertices))
	 (segment-vertex-1 (lex-min vertex-1 vertex-2))
	 (segment-vertex-2 (lex-max vertex-1 vertex-2)))
    (cond ((and (equal? triangle-vertex-1 segment-vertex-1)
		(equal? triangle-vertex-2 segment-vertex-2))
	   triangle-vertex-3)
	  ((and (equal? triangle-vertex-1 segment-vertex-1)
		(equal? triangle-vertex-3 segment-vertex-2))
	   triangle-vertex-2)
	  ((and (equal? triangle-vertex-2 segment-vertex-1)
		(equal? triangle-vertex-3 segment-vertex-2))
	   triangle-vertex-1)
	  (else (error "input points not triangle-vertices")))))

(define (opposite-triangle-pointer-to-vertex vertex triangle)
  (cond ((equal? vertex (triangle-vertex-1 triangle))
	 (triangle-triangle-pointer-1 triangle))
	((equal? vertex (triangle-vertex-2 triangle))
	 (triangle-triangle-pointer-2 triangle))
	((equal? vertex (triangle-vertex-3 triangle))
	 (triangle-triangle-pointer-3 triangle))
	(else (error "input point not triangle-vertex"))))

(define (opposite-edge-pointer-to-vertex vertex triangle)
  (cond ((equal? vertex (triangle-vertex-1 triangle))
	 (triangle-edge-pointer-1 triangle))
	((equal? vertex (triangle-vertex-2 triangle))
	 (triangle-edge-pointer-2 triangle))
	((equal? vertex (triangle-vertex-3 triangle))
	 (triangle-edge-pointer-3 triangle))
	(else (error "input point not triangle-vertex"))))

(define (point-support point triangulation)
  (let loop ((counter (- (collection-length triangulation) 1)))
    (if (< counter 0)
      -1
      (let ((triangle (collection-ref triangulation counter)))
	(if triangle
	  (let ((vertex-1 (triangle-vertex-1 triangle))
		(vertex-2 (triangle-vertex-2 triangle))
		(vertex-3 (triangle-vertex-3 triangle)))
	    (cond ((point-inside-triangle? point vertex-1 vertex-2 vertex-3)
		   counter)
		  ((point-inside-edge? point vertex-1 vertex-2)
		   (make-edge
		     vertex-1 vertex-2
		     (min counter (triangle-triangle-pointer-3 triangle))
		     (max counter (triangle-triangle-pointer-3 triangle))
		     #f #t))
		  ((point-inside-edge? point vertex-2 vertex-3)
		   (make-edge
		     vertex-2 vertex-3
		     (min counter (triangle-triangle-pointer-1 triangle))
		     (max counter (triangle-triangle-pointer-1 triangle))
		     #f #t))
		  ((point-inside-edge? point vertex-3 vertex-1)
		   (make-edge
		     vertex-3 vertex-1
		     (min counter (triangle-triangle-pointer-2 triangle))
		     (max counter (triangle-triangle-pointer-2 triangle))
		     #f #t))
		  (else (loop (- counter 1)))))
	  (loop (- counter 1)))))))

(define (split-triangle-by-point point support edges triangulation)
  (let* ((old-triangle-pointer support)
	 (triangle (collection-ref triangulation support))
	 (vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (new-triangle-pointer-1 support)
	 (new-triangle-pointer-2 (collection-length triangulation))
	 (new-triangle-pointer-3 (+ new-triangle-pointer-2 1))
	 (triangle-pointer-1 (triangle-triangle-pointer-1 triangle))
	 (triangle-pointer-2 (triangle-triangle-pointer-2 triangle))
	 (triangle-pointer-3 (triangle-triangle-pointer-3 triangle))
	 (edge-pointer-1 (triangle-edge-pointer-1 triangle))
	 (edge-pointer-2 (triangle-edge-pointer-2 triangle))
	 (edge-pointer-3 (triangle-edge-pointer-3 triangle))
	 (edge-pointer-1p (collection-length edges))
	 (edge-pointer-2p (+ edge-pointer-1p 1))
	 (edge-pointer-3p (+ edge-pointer-2p 1))
	 (edge-1p (make-edge
		    vertex-1 point
		    (min new-triangle-pointer-1 new-triangle-pointer-3)
		    (max new-triangle-pointer-1 new-triangle-pointer-3) 
		    #f #t))
	 (edge-2p (make-edge
		    vertex-2 point
		    (min new-triangle-pointer-1 new-triangle-pointer-2)
		    (max new-triangle-pointer-1 new-triangle-pointer-2) 
		    #f #t))
	 (edge-3p (make-edge
		    vertex-3 point
		    (min new-triangle-pointer-2 new-triangle-pointer-3)
		    (max new-triangle-pointer-2 new-triangle-pointer-3) 
		    #f #t))
	 (new-triangle-1 (make-triangle vertex-1 vertex-2 point
					edge-pointer-2p 
					edge-pointer-1p 
					edge-pointer-3
					new-triangle-pointer-2 
					new-triangle-pointer-3
					triangle-pointer-3))
	 (new-triangle-2 (make-triangle vertex-2 vertex-3 point
					edge-pointer-3p 
					edge-pointer-2p 
					edge-pointer-1
					new-triangle-pointer-3 
					new-triangle-pointer-1
					triangle-pointer-1))
	 (new-triangle-3 (make-triangle vertex-3 vertex-1 point
					edge-pointer-1p 
					edge-pointer-3p 
					edge-pointer-2
					new-triangle-pointer-1 
					new-triangle-pointer-2
					triangle-pointer-2)))
    ;; First, a few changes in the pointers of edges and triangles.
    (if (> triangle-pointer-1 -1)
      (change-triangle-pointers
	(collection-ref triangulation triangle-pointer-1)
	old-triangle-pointer new-triangle-pointer-2))
    (if (> triangle-pointer-2 -1)
      (change-triangle-pointers 
	(collection-ref triangulation triangle-pointer-2)
	old-triangle-pointer new-triangle-pointer-3))
    (change-edge-pointers
      (collection-ref edges edge-pointer-1)
      old-triangle-pointer new-triangle-pointer-2)
    (change-edge-pointers
      (collection-ref edges edge-pointer-2)
      old-triangle-pointer new-triangle-pointer-3)
    (change-edge-pointers
      (collection-ref edges edge-pointer-3)
      old-triangle-pointer new-triangle-pointer-1)
    ;; And now, get rid of the old triangle and append the three new ones...
    (replace-element-in-collection
      triangulation old-triangle-pointer 
      new-triangle-1 new-triangle-2 new-triangle-3)
    ;; ...and append the three new edges.
    (add-elements-to-collection edges edge-3p edge-2p edge-1p)
    (legalize-edge vertex-1 vertex-2 
		   new-triangle-pointer-1 triangle-pointer-3
		   edges triangulation)
    (legalize-edge vertex-2 vertex-3 
		   new-triangle-pointer-2 triangle-pointer-1
		   edges triangulation)
    (legalize-edge vertex-3 vertex-1
		   new-triangle-pointer-3 triangle-pointer-2
		   edges triangulation)))

(define (split-edge-by-point point support edges triangulation)
  (let* ((vertex-1 (edge-from support))
	 (vertex-3 (edge-to support))
	 (triangle-pointer-up (edge-triangle-pointer-1 support))
	 (triangle-pointer-down (edge-triangle-pointer-2 support))
	 (triangle-up (collection-ref triangulation triangle-pointer-up))
	 (triangle-down (collection-ref triangulation triangle-pointer-down))
	 (vertex-2 (the-other-vertex vertex-1 vertex-3 triangle-up))
	 (vertex-4 (the-other-vertex vertex-1 vertex-3 triangle-down))
	 (triangle-pointer-1 
	   (opposite-triangle-pointer-to-vertex vertex-3 triangle-up))
	 (triangle-pointer-2
	   (opposite-triangle-pointer-to-vertex vertex-1 triangle-up))
	 (triangle-pointer-3
	   (opposite-triangle-pointer-to-vertex vertex-1 triangle-down))
	 (triangle-pointer-4
	   (opposite-triangle-pointer-to-vertex vertex-3 triangle-down))
	 (new-triangle-pointer-1 triangle-pointer-up)
	 (new-triangle-pointer-2 (collection-length triangulation))
	 (new-triangle-pointer-3 (+ new-triangle-pointer-2 1))
	 (new-triangle-pointer-4 triangle-pointer-down)
	 (edge-index 
	   (opposite-edge-pointer-to-vertex vertex-2 triangle-up))
	 (edge-pointer-1 
	   (opposite-edge-pointer-to-vertex vertex-3 triangle-up))
	 (edge-pointer-2 
	   (opposite-edge-pointer-to-vertex vertex-1 triangle-up))
	 (edge-pointer-3 
	   (opposite-edge-pointer-to-vertex vertex-1 triangle-down))
	 (edge-pointer-4
	   (opposite-edge-pointer-to-vertex vertex-3 triangle-down))
	 (edge-pointer-1p edge-index)
	 (edge-pointer-2p (collection-length edges))
	 (edge-pointer-3p (+ edge-pointer-2p 1))
	 (edge-pointer-4p (+ edge-pointer-2p 2))
	 (edge-1p (make-edge
		    vertex-1 point
		    (min new-triangle-pointer-1 new-triangle-pointer-4)
		    (max new-triangle-pointer-1 new-triangle-pointer-4) 
		    #f #t))
	 (edge-2p (make-edge
		    vertex-2 point
		    (min new-triangle-pointer-2 new-triangle-pointer-1)
		    (max new-triangle-pointer-2 new-triangle-pointer-1) 
		    #f #t))
	 (edge-3p (make-edge
		    vertex-3 point
		    (min new-triangle-pointer-3 new-triangle-pointer-2)
		    (max new-triangle-pointer-3 new-triangle-pointer-2) 
		    #f #t))
	 (edge-4p (make-edge
		    vertex-4 point
		    (min new-triangle-pointer-3 new-triangle-pointer-4)
		    (max new-triangle-pointer-3 new-triangle-pointer-4) 
		    #f #t))
	 (new-triangle-1 
	   (make-triangle vertex-1 vertex-2 point
			  edge-pointer-2p 
			  edge-pointer-1p 
			  edge-pointer-1
			  new-triangle-pointer-2
			  new-triangle-pointer-4
			  triangle-pointer-1))
	 (new-triangle-2 
	   (make-triangle vertex-2 vertex-3 point
			  edge-pointer-3p 
			  edge-pointer-2p 
			  edge-pointer-2
			  new-triangle-pointer-3
			  new-triangle-pointer-1
			  triangle-pointer-2))
	 (new-triangle-3 
	   (make-triangle vertex-3 vertex-4 point
			  edge-pointer-4p 
			  edge-pointer-3p 
			  edge-pointer-3
			  new-triangle-pointer-4
			  new-triangle-pointer-2
			  triangle-pointer-3))
	 (new-triangle-4 
	   (make-triangle vertex-4 vertex-1 point
			  edge-pointer-1p 
			  edge-pointer-4p 
			  edge-pointer-4
			  new-triangle-pointer-1
			  new-triangle-pointer-3
			  triangle-pointer-4)))
    ;; First, a few changes in the pointers of edges and triangles.
    (if (> triangle-pointer-2 -1)
      (change-triangle-pointers 
	(collection-ref triangulation triangle-pointer-2)
	triangle-pointer-up new-triangle-pointer-2))
    (if (> triangle-pointer-3 -1)
      (change-triangle-pointers 
	(collection-ref triangulation triangle-pointer-3)
	triangle-pointer-down new-triangle-pointer-3))
    (change-edge-pointers
      (collection-ref edges edge-pointer-1)
      triangle-pointer-up new-triangle-pointer-1)
    (change-edge-pointers
      (collection-ref edges edge-pointer-2)
      triangle-pointer-up new-triangle-pointer-2)
    (change-edge-pointers
      (collection-ref edges edge-pointer-3)
      triangle-pointer-down new-triangle-pointer-3)
    (change-edge-pointers
      (collection-ref edges edge-pointer-4)
      triangle-pointer-down new-triangle-pointer-4)
    ;; And now, get rid of old triangles and append the four new ones...
    (replace-element-in-collection
      triangulation triangle-pointer-up new-triangle-1 new-triangle-2)
    (replace-element-in-collection
      triangulation triangle-pointer-down new-triangle-4 new-triangle-3)
    ;; ...and get rid of the old edge and append the four new ones.
    (replace-element-in-collection
      edges edge-index edge-1p edge-2p edge-3p edge-4p)
    (legalize-edge vertex-1 vertex-2 
		   new-triangle-pointer-1 triangle-pointer-1
		   edges triangulation)
    (legalize-edge vertex-2 vertex-3 
		   new-triangle-pointer-2 triangle-pointer-2
		   edges triangulation)
    (legalize-edge vertex-3 vertex-4 
		   new-triangle-pointer-3 triangle-pointer-3
		   edges triangulation)
    (legalize-edge vertex-4 vertex-1 
		   new-triangle-pointer-4 triangle-pointer-4
		   edges triangulation)))

(define (legalize-edge vertex-1 vertex-3 
		       triangle-pointer-up triangle-pointer-down
		       edges triangulation)
  (if (> triangle-pointer-down -1)
    (let* ((triangle-up (collection-ref triangulation triangle-pointer-up))
	   (vertex-2 (the-other-vertex vertex-1 vertex-3 triangle-up))
	   (triangle-down 
	     (collection-ref triangulation triangle-pointer-down))
	   (vertex-4 (the-other-vertex vertex-1 vertex-3 triangle-down)))
      (if (point-inside-circumcircle? vertex-4 vertex-1 vertex-2 vertex-3)
	(let* ((triangle-pointer-1
		 (opposite-triangle-pointer-to-vertex vertex-3 triangle-up))
	       (triangle-pointer-2
		 (opposite-triangle-pointer-to-vertex vertex-1 triangle-up))
	       (triangle-pointer-3
		 (opposite-triangle-pointer-to-vertex vertex-1 triangle-down))
	       (triangle-pointer-4
		 (opposite-triangle-pointer-to-vertex vertex-3 triangle-down))
	       (triangle-pointer-left triangle-pointer-up)
	       (triangle-pointer-right triangle-pointer-down)
	       (edge-index 
		 (opposite-edge-pointer-to-vertex vertex-2 triangle-up))
	       (other-edge-index
		 (opposite-edge-pointer-to-vertex vertex-4 triangle-down))
	       (edge-pointer-1 
		 (opposite-edge-pointer-to-vertex vertex-3 triangle-up))
	       (edge-pointer-2
		 (opposite-edge-pointer-to-vertex vertex-1 triangle-up))
	       (edge-pointer-3
		 (opposite-edge-pointer-to-vertex vertex-1 triangle-down))
	       (edge-pointer-4
		 (opposite-edge-pointer-to-vertex vertex-3 triangle-down))
	       (edge-24 
		 (make-edge
		   vertex-2 vertex-4
		   (min triangle-pointer-left triangle-pointer-right)
		   (max triangle-pointer-left triangle-pointer-right) 
		   #f #t))
	       (triangle-left
		 (make-triangle vertex-1 vertex-2 vertex-4
				edge-index
				edge-pointer-4
				edge-pointer-1
				triangle-pointer-right
				triangle-pointer-4
				triangle-pointer-1))
	       (triangle-right
		 (make-triangle vertex-2 vertex-3 vertex-4
				edge-pointer-3 
				edge-index
				edge-pointer-2
				triangle-pointer-3
				triangle-pointer-left
				triangle-pointer-2)))
	  ;; First, a few changes in the pointers of edges and triangles.
	  (if (> triangle-pointer-2 -1)
	    (change-triangle-pointers
	      (collection-ref triangulation triangle-pointer-2)
	      triangle-pointer-up triangle-pointer-right))
	  (if (> triangle-pointer-4 -1)
	    (change-triangle-pointers
	      (collection-ref triangulation triangle-pointer-4)
	      triangle-pointer-down triangle-pointer-left))
	  (change-edge-pointers
	    (collection-ref edges edge-pointer-2)
	    triangle-pointer-up triangle-pointer-right)
	  (change-edge-pointers
	    (collection-ref edges edge-pointer-4)
	    triangle-pointer-down triangle-pointer-left)
	  ;; Get rid of the old two triangles, append the two new ones...
	  (replace-element-in-collection
	    triangulation triangle-pointer-up triangle-left)
	  (replace-element-in-collection
	    triangulation triangle-pointer-down triangle-right)
	  ;; ... and get rid of the illegal edge, and append the new one.
	  (replace-element-in-collection edges edge-index edge-24)
	  (legalize-edge vertex-1 vertex-4
			 triangle-pointer-left triangle-pointer-4
			 edges triangulation)
	  (legalize-edge vertex-3 vertex-4
			 triangle-pointer-right triangle-pointer-3
			 edges triangulation))))))

(define (insert-point-in-triangulation point edges triangulation)
  (let ((support (point-support point triangulation)))
    (if (edge? support)
      (split-edge-by-point point support edges triangulation)
      (if (> support -1)
	(split-triangle-by-point point support edges triangulation)
	(error "the point seems to be nowhere!")))))



;; Auxiliary functions for CONSTRAINED-STEINER-DELAUNAY-TRIANGULATION:

(define (is-segment-in-edges? segment set)
  (let ((vertex-1 (edge-from segment))
	(vertex-2 (edge-to segment)))
    (let loop ((counter (- (collection-length set) 1)))
      (if (< counter 0) #f
	(let ((edge (collection-ref set counter)))
	  (if edge
	    (let ((vertex-from (edge-from edge))
		  (vertex-to (edge-to edge)))
	      (if (or (and (equal? vertex-from vertex-1)
			   (equal? vertex-to vertex-2))
		      (and (equal? vertex-to vertex-1)
			   (equal? vertex-from vertex-2)))
		(edge-boundary?-set! edge #t)
		(loop (- counter 1))))
	    (loop (- counter 1))))))))
 


;; Auxiliary functions for Angle Handling:

(define (remove-exterior-triangles-from-triangulation 
	  holes edges triangulation)

  (define (remove-adjacent-exterior-triangles 
	    triangle-pointer edges triangulation)
    (if (> triangle-pointer -1)
      (let ((triangle (collection-ref triangulation triangle-pointer)))
	(if triangle
	  (let* ((edge-pointer-1 (triangle-edge-pointer-1 triangle))
		 (edge-pointer-2 (triangle-edge-pointer-2 triangle))
		 (edge-pointer-3 (triangle-edge-pointer-3 triangle))
		 (triangle-pointer-1 (triangle-triangle-pointer-1 triangle))
		 (triangle-pointer-2 (triangle-triangle-pointer-2 triangle))
		 (triangle-pointer-3 (triangle-triangle-pointer-3 triangle))
		 (edge-1 (collection-ref edges edge-pointer-1))
		 (edge-2 (collection-ref edges edge-pointer-2))
		 (edge-3 (collection-ref edges edge-pointer-3)))
	    (replace-element-in-collection triangulation triangle-pointer #f)
	    (change-edge-pointers edge-1 triangle-pointer -1)
	    (change-edge-pointers edge-2 triangle-pointer -1)
	    (change-edge-pointers edge-3 triangle-pointer -1)
	    (if (> triangle-pointer-1 -1)
	      (change-triangle-pointers 
		(collection-ref triangulation triangle-pointer-1)
		triangle-pointer -1))
	    (if (> triangle-pointer-2 -1)
	      (change-triangle-pointers
		(collection-ref triangulation triangle-pointer-2)
		triangle-pointer -1))
	    (if (> triangle-pointer-3 -1)
	      (change-triangle-pointers
		(collection-ref triangulation triangle-pointer-3)
		triangle-pointer -1))
	    (if (not (edge-boundary? edge-1))
	      (begin 
		(edge-useful?-set! edge-1 #f)
		(remove-adjacent-exterior-triangles 
		  triangle-pointer-1 edges triangulation)))
	    (if (not (edge-boundary? edge-2))
	      (begin
		(edge-useful?-set! edge-2 #f)
		(remove-adjacent-exterior-triangles 
		  triangle-pointer-2 edges triangulation)))
	    (if (not (edge-boundary? edge-3))
	      (begin
		(edge-useful?-set! edge-3 #f)
		(remove-adjacent-exterior-triangles 
		  triangle-pointer-3 edges triangulation))))))))

  (do ((counter (- (collection-length holes) 1) (- counter 1)))
    ((< counter 0) triangulation)
    (let* ((hole (collection-ref holes counter))
	   (support (point-support hole triangulation)))

      (if (edge? support)
	(remove-adjacent-exterior-triangles
	  (max (edge-triangle-pointer-2 support) 
	       (edge-triangle-pointer-1 support))
	  edges triangulation)
	(remove-adjacent-exterior-triangles support edges triangulation)))))




(define (modify-triangles-with-minimum-angle-greater-than-threshold
	  threshold-angle segments edges triangulation)
  (let ((thinner-second-class-triangle 
	  (triangulation-thinner-second-class-triangle 
	    threshold-angle edges triangulation)))
    (if thinner-second-class-triangle
      (let* ((circumcenter 
	       (triangle-circumcenter thinner-second-class-triangle))
	     (first-encroached-edge-pointer
	       (first-pointer-encroached-upon-by-point 
		 circumcenter edges)))
	(display "Circumcenter: ")
	(pp circumcenter)
	(if (> first-encroached-edge-pointer -1)
	  (let* ((first-encroached-edge 
		   (collection-ref edges first-encroached-edge-pointer))
		 (first-encroached-segment-pointer
		   (first-pointer-encroached-upon-by-point 
		     circumcenter segments))
		 (first-encroached-segment 
		   (collection-ref segments first-encroached-segment-pointer))
		 (vertex-1 (edge-from first-encroached-segment))
		 (vertex-2 (edge-to first-encroached-segment))
		 (midpoint (midpoint vertex-1 vertex-2))
		 (new-segment-1 (make-edge vertex-1 midpoint #f #f #t #t))
		 (new-segment-2 (make-edge midpoint vertex-2 #f #f #t #t)))
	    (replace-element-in-collection
	      segments first-encroached-segment-pointer 
	      new-segment-1 new-segment-2)
	    (if (or (> (edge-triangle-pointer-1 first-encroached-edge) -1)
		    (> (edge-triangle-pointer-2 first-encroached-edge) -1))
	      (begin
		(split-edge-by-point
		  midpoint (collection-ref edges first-encroached-edge-pointer)
		  edges triangulation)
		(include-segments-in-triangulation segments edges triangulation)
		(modify-triangles-with-minimum-angle-greater-than-threshold
		  threshold-angle segments edges triangulation))
	      (modify-triangles-with-minimum-angle-greater-than-threshold
		threshold-angle segments edges triangulation)))
	  (let ((support (point-support circumcenter)))
	    (if (edge? support)
	      (begin
		(split-edge-by-point circumcenter support edges triangulation)
		(include-segments-in-triangulation segments edges triangulation)
		(modify-triangles-with-minimum-angle-greater-than-threshold
		  threshold-angle segments edges triangulation))
	      (begin
		(split-triangle-by-point 
		  circumcenter support edges triangulation)
		(include-segments-in-triangulation segments edges triangulation)
		(modify-triangles-with-minimum-angle-greater-than-threshold
		  threshold-angle segments edges triangulation))))))
      (cons edges triangulation))))

(define (triangulation-thinner-second-class-triangle 
	  threshold-angle edges triangulation)
  ;; We define a "second class" triangle in a triangulation with respect to a
  ;; threshold angle ß as either any non-corner triangle with minimum angle
  ;; less than ß, or any corner-triangle with minimum angle less than ß, but
  ;; corner angle greater than or egual to 2ß.
  (let loop ((counter (- (collection-length triangulation) 1))
	     (thinner-triangle #f)
	     (thinner-angle +inf.))
    (if (< counter 0)
      (begin
	(display "\nThinner triangle:\n")
	(pp thinner-triangle)
	(display "Thinner angle: ")
	(pp thinner-angle)
	thinner-triangle)
      (let ((triangle (collection-ref triangulation counter)))
	(if triangle
	  (let ((minimum-angle (triangle-minimum-angle triangle)))
	    (if (> minimum-angle thinner-angle)
	      (loop (- counter 1) thinner-triangle thinner-angle)
	      (if (is-corner-triangle? triangle edges)
		(if (>= (triangle-corner-angle triangle edges) 
			(* 2 threshold-angle))
		  (loop (- counter 1) triangle minimum-angle)
		  (loop (- counter 1) thinner-triangle thinner-angle))
		(loop (- counter 1) triangle minimum-angle))))
	  (loop (- counter 1) thinner-triangle thinner-angle))))))

(define (first-pointer-encroached-upon-by-point point edge-collection)
  (let loop ((counter (- (collection-length edge-collection) 1)))
    (if (< counter 0)
      #f
      (let* ((edge (collection-ref edge-collection counter))
	     (vertex-1 (edge-from edge))
	     (vertex-2 (edge-to edge))
	     (midpoint (midpoint vertex-1 vertex-2)))
	(if (< (distance midpoint point) (distance midpoint vertex-1))
	  counter
	  (loop (- counter 1)))))))
	



;; Basic functions: 

(define (lex-min point-1 point-2)
  (declare (flonum))
  (if (or (< (point-x point-1) (point-x point-2))
	  (and (= (point-x point-1) (point-x point-2))
	       (< (point-y point-1) (point-y point-2))))
    point-1 point-2))

(define (lex-max point-1 point-2)
  (declare (flonum))
  (if (equal? point-1 (lex-min point-1 point-2))
    point-2 point-1))

(define (lex-sort point-1 point-2 point-3)
  (declare (flonum))
  (cond ((equal? point-2 (lex-min point-1 point-2))
	 (lex-sort point-2 point-1 point-3))
	((equal? point-3 (lex-min point-2 point-3))
	 (lex-sort point-1 point-3 point-2))
	(else (list point-1 point-2 point-3))))
	  
(define (segment-functional vertex-1 vertex-2)
  (declare (flonum))
  (lambda (x y)
    (+ (* (- (point-x vertex-2) (point-x vertex-1)) 
	  (- y (point-y vertex-1)))
       (* (- (point-x vertex-1) x) 
	  (- (point-y vertex-2) (point-y vertex-1))))))

(define (circle-functional vertex-1 vertex-2 vertex-3)
  (declare (flonum))
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
			 1 x3 y3))))))

(define (point-inside-triangle? point vertex-1 vertex-2 vertex-3)
  (declare (flonum))
  (let ((functional-1 (segment-functional vertex-1 vertex-2))
	(functional-2 (segment-functional vertex-2 vertex-3))
	(functional-3 (segment-functional vertex-3 vertex-1))
	(x (point-x point))
	(y (point-y point)))
    (let ((test-1 (functional-1 x y))
	  (test-2 (functional-2 x y))
	  (test-3 (functional-3 x y)))
      (if (positive? test-1)
	(and (positive? test-2) (positive? test-3))
	(and (negative? test-1) (negative? test-2) (negative? test-3))))))

(define (point-inside-edge? point vertex-1 vertex-2)
  (declare (flonum))
  (let ((functional (segment-functional vertex-1 vertex-2))
	 (x (point-x point))
	 (y (point-y point)))
    (if (zero? (functional x y))
      (let ((midpoint (midpoint vertex-1 vertex-2)))
	(< (distance point midpoint) (distance vertex-1 midpoint)))
      #f)))

(define (point-inside-circumcircle? point vertex-1 vertex-2 vertex-3)
  (declare (flonum))
  (let ((x0 (point-x point))
	(y0 (point-y point))
	(x1 (point-x vertex-1))
	(y1 (point-y vertex-1))
	(circle-test (circle-functional vertex-1 vertex-2 vertex-3))
	(orientation-test 
	  (segment-functional vertex-2 vertex-3)))
    (positive? 
      (* (circle-test x0 y0) 
	 (orientation-test x1 y1)))))

(define (triangle-circumcenter triangle)
  (let* ((vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (x1 (point-x vertex-1))
	 (y1 (point-y vertex-1))
	 (x2 (point-x vertex-2))
	 (y2 (point-y vertex-2))
	 (x3 (point-x vertex-3))
	 (y3 (point-y vertex-3))
	 (solution (solve-system 
		     (- x1 x2) (- y1 y2) (+ (* (/ (+ x1 x2) 2.0)
					       (- x1 x2))
					    (* (/ (+ y1 y2) 2.0)
					       (- y1 y2)))
		     (- x2 x3) (- y2 y3) (+ (* (/ (+ x2 x3) 2.0)
					       (- x2 x3))
					    (* (/ (+ y2 y3) 2.0)
					       (- y2 y3))))))
    (make-point (car solution) (cdr solution))))

(define (triangle-minimum-angle triangle)
  (if triangle
    (let ((vertex-1 (triangle-vertex-1 triangle))
	  (vertex-2 (triangle-vertex-2 triangle))
	  (vertex-3 (triangle-vertex-3 triangle)))
      (min (wedge-angle vertex-3 vertex-1 vertex-2)
	   (wedge-angle vertex-1 vertex-2 vertex-3)
	   (wedge-angle vertex-2 vertex-3 vertex-1)))))

(define (is-corner-triangle? triangle edges)
  (let ((edge-1 (collection-ref edges (triangle-edge-pointer-1 triangle)))
	(edge-2 (collection-ref edges (triangle-edge-pointer-2 triangle)))
	(edge-3 (collection-ref edges (triangle-edge-pointer-3 triangle))))
    (or (and (edge-boundary? edge-1)
	     (edge-boundary? edge-2))
	(and (edge-boundary? edge-1)
	     (edge-boundary? edge-3))
	(and (edge-boundary? edge-2)
	     (edge-boundary? edge-3)))))

(define (triangle-corner-angle triangle edges)
  (let ((vertex-1 (triangle-vertex-1 triangle))
	(vertex-2 (triangle-vertex-2 triangle))
	(vertex-3 (triangle-vertex-3 triangle))
	(edge-1 (collection-ref edges (triangle-edge-pointer-1 triangle)))
	(edge-2 (collection-ref edges (triangle-edge-pointer-2 triangle)))
	(edge-3 (collection-ref edges (triangle-edge-pointer-3 triangle))))
    (cond ((and (edge-boundary? edge-1) (edge-boundary? edge-2))
	   (wedge-angle vertex-2 vertex-3 vertex-1))
	  ((and (edge-boundary? edge-1) (edge-boundary? edge-3))
	   (wedge-angle vertex-1 vertex-2 vertex-3))
	  ((and (edge-boundary? edge-2) (edge-boundary? edge-3))
	   (wedge-angle vertex-3 vertex-1 vertex-2)))))

(define (wedge-angle side-point-1 corner-point side-point-2)
  (declare (flonum))
  (let ((candidate (abs (* 180 (- (angle (- side-point-2 corner-point)) 
				  (angle (- side-point-1 corner-point)))
			   (/ 1 pi)))))
    (if (> candidate 180)
      (- 360 candidate)
      candidate)))

(define (distance point-1 point-2)
  (magnitude (- point-1 point-2)))

(define  (add-squares x y) (declare (flonum)) (+ (* x x) (* y y)))

(define (determinant a b c
		     d e f 
		     g h i)
  ;;    |         |    
  ;;    | a  b  c |     |      |     |      |     |      |  
  ;;    |         |     | e  f |     | d  f |     | d  e |
  ;;    | d  e  f | = a |      | - b |      | + c |      |
  ;;    |         |     | h  i |     | g  i |     | g  h | 
  ;;    | g  h  i |     |      |     |      |     |      |  
  ;;    |         |    
  (declare (flonum))
  (- (+ (* a e i) (* b f g) (* c d h))
     (+ (* c e g) (* b d i) (* a f h))))

(define (solve-system a00 a01 b0
		      a10 a11 b1)
  ;;      |         | |   |   |    |
  ;;      | a00 a01 | | x |   | b0 |
  ;;      |         | |   | = |    |
  ;;      | a10 a11 | | y |   | b1 |
  ;;      |         | |   |   |    |
  (declare (flonum))
  (let ((denominator (- (* a00 a11) (* a10 a01)))
	(numerator-1 (- (* b0  a11) (* b1  a01)))
	(numerator-2 (- (* a00 b1 ) (* a10 b0 ))))
    (cons (/ numerator-1 denominator 1.0) 
	  (/ numerator-2 denominator 1.0))))

(define pi 3.14159265358979)

;; Functions for Displaying and debugging.

(define (write-triangulation edges file-name)

  (define (write-edge edge)
    (if (edge-useful? edge)
      (let ((x1 (point-x (edge-from edge)))
	    (y1 (point-y (edge-from edge)))
	    (x2 (point-x (edge-to edge)))
	    (y2 (point-y (edge-to edge))))
	(display x1)(display " ")(display y1)(newline)
	(display x2)(display " ")(display y2)(newline)(newline))))

  (display "output: ")(display file-name)(newline)
  (with-output-to-file
    file-name
    (lambda ()
      (do ((counter (- (collection-length edges) 1) (- counter 1)))
	((< counter 0) (newline))
	(let ((edge (collection-ref edges counter)))
	  (if edge
	    (write-edge (collection-ref edges counter))))))))

(define (gnu-triangulation triangulation file-name)
  
  (define (write-triangle triangle)
    (let* ((x1 (point-x (triangle-vertex-1 triangle)))
	   (y1 (point-y (triangle-vertex-1 triangle)))
	   (x2 (point-x (triangle-vertex-2 triangle)))
	   (y2 (point-y (triangle-vertex-2 triangle)))
	   (x3 (point-x (triangle-vertex-3 triangle)))
	   (y3 (point-y (triangle-vertex-3 triangle))))
      (display x1)(display " ")(display y1)(newline)
      (display x2)(display " ")(display y2)(newline)
      (display x3)(display " ")(display y3)(newline)
      (display x1)(display " ")(display y1)(newline)(newline)))

  (display "output: ")(display file-name)(newline)
  (with-output-to-file
    file-name
    (lambda ()
      (do ((counter (- (collection-length triangulation) 1) (- counter 1)))
	((< counter 0) (newline))
	(let ((triangle (collection-ref triangulation counter)))
	  (if triangle (write-triangle triangle)))))))

(define (midpoint vertex-1 vertex-2)
  (/ (+ vertex-1 vertex-2) 2))

(define (triangle-angles triangle)
  (if triangle
    (let ((vertex-1 (triangle-vertex-1 triangle))
	  (vertex-2 (triangle-vertex-2 triangle))
	  (vertex-3 (triangle-vertex-3 triangle)))
      (display 
	(list (wedge-angle vertex-3 vertex-1 vertex-2)
	      (wedge-angle vertex-1 vertex-2 vertex-3)
	      (wedge-angle vertex-2 vertex-3 vertex-1)))
      (newline))))

;; Misc:

(define DT delaunay-triangulation)
(define CSDT constrained-steiner-delaunay-triangulation)
(define MACSDT minimum-angle-constrained-steiner-delaunay-triangulation)

(define a (read-file "PolyFiles/nasa.poly"))
(define b 
  (time (minimum-angle-constrained-steiner-delaunay-triangulation 30 a)))
(write-triangulation (cdr b) "nasa.gnu")
