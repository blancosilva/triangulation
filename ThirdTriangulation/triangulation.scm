;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  First Step: Delaunay Triangulation of the
;;               convex hull of input polygon.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; CONVEX-HULL-TRIANGULATION creates a naive triangulation with the
;; boundary vertices (the last four) and the first vertex of the polygon
;; It also inserts the rest of the vertices one at a time, updating the
;; triangulation.
;;
(define (convex-hull-triangulation polygon-vertices)
  (let loop 
    ((vertex-entry 1) 
     (triangulation 
       (naive-triangulation (vector-length polygon-vertices))))
      (if (= vertex-entry (- (vector-length polygon-vertices) 4))
	triangulation
	(loop (+ 1 vertex-entry) 
	      (introduce-point-in-triangulation 
		vertex-entry polygon-vertices triangulation)))))

;; INTRODUCE-POINT-IN-TRIANGULATION is the main function.  It has to be
;; independent of (convex-hull-triangulation) since it is used also as
;; building block in (introduce-segments-in-triangulation).  
;; Given the input point, looks for a triangle where it is contained.
;; If it finds one, performs a call to (split-triangle); otherwise, it
;; assumes it is in the interior of an edge, and calls (split-edge).
;; Either way, the output is a delaunay-triangulation.
;;
(define (introduce-point-in-triangulation 
	  point-entry polygon-vertices triangulation)
  (let loop ((triangle-entry 0))
    (if (= triangle-entry (vector-length triangulation))
      (split-edge point-entry polygon-vertices triangulation)
      (let* ((triangle
	       (vector-ref triangulation triangle-entry))
	     (vertex-1 
	       (vector-ref polygon-vertices 
			   (triangle-vertex-1 triangle)))
	     (vertex-2
	       (vector-ref polygon-vertices
			   (triangle-vertex-2 triangle)))
	     (vertex-3
	       (vector-ref polygon-vertices
			   (triangle-vertex-3 triangle)))
	     (point (vector-ref polygon-vertices point-entry)))
	(if (point-inside-triangle?
	      point vertex-1 vertex-2 vertex-3)
	  (split-triangle point-entry triangle-entry
			  triangulation polygon-vertices)
	  (loop (+ 1 triangle-entry)))))))
	     
;; SPLIT-TRIANGLE replaces triangle in position triangle-entry with one
;; of the three new triangles, and appends the other two at the end of 
;; the triangulation vector.  It legalizes edges if necessary, and 
;; returns the updated triangulation.
;;                          
;;               _ _ _ _ _ _v2 _ _ _ _ _ _
;;                          /\
;;              |          /..\          |
;;                        / .. \
;;              |  pT12  /  ..  \   pT23 |
;;                      /   ..   \
;;              |      / pT ..  s \      |
;;                    /     .@.    \
;;              |    /    .. p ..   \    |
;;                  /   ..       ..  \ 
;;              |  /  ..           .. \  |
;;                / ..      s+1      ..\ 
;;              |/.                    .\|
;;           v1 *------------------------*v3
;;                +                    +
;;                   +    pT31      + 
;;                      +        +
;;                         +  + 
;;                            
;;
(define (split-triangle point-entry triangle-entry 
			triangulation polygon-vertices)
  (let* ((triangle (vector-ref triangulation triangle-entry))
	 (vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (pointer-1-2 (triangle-pointer-1-2 triangle))
	 (pointer-2-3 (triangle-pointer-2-3 triangle))
	 (pointer-3-1 (triangle-pointer-3-1 triangle))
	 (size (vector-length triangulation))
	 (new-triangle-1 
	   (make-triangle vertex-1 vertex-2 point-entry
			  pointer-1-2 size (+ size 1) 0))
	 (new-triangle-2
	   (make-triangle vertex-2 vertex-3 point-entry
			  pointer-2-3 (+ size 1) triangle-entry 0))
	 (new-triangle-3 
	   (make-triangle vertex-3 vertex-1 point-entry
			  pointer-3-1 triangle-entry size 0)))
    (vector-set! triangulation triangle-entry new-triangle-1)
    (change-pointers 
      pointer-2-3 triangle-entry size triangulation)
    (change-pointers 
      pointer-3-1 triangle-entry (+ size 1) triangulation)
    (let ((new-triangulation
	    (list->vector
	      (append (vector->list triangulation)
		      (list new-triangle-2 new-triangle-3)))))
      (legalize-edge vertex-1 vertex-2 triangle-entry 
		     new-triangulation polygon-vertices)
      (legalize-edge vertex-2 vertex-3 size 
		     new-triangulation polygon-vertices)
      (legalize-edge vertex-3 vertex-1 (+ size 1) 
		     new-triangulation polygon-vertices)
      new-triangulation)))
	   
;; SPLIT-EDGE "knows" that there is at least one edge in the triangula-
;; tion containing the input point, but it doesn't know which one it is.
;; So first, it looks for that edge and the two adjacent triangles.
;; Then replaces the two triangles sharing the edge that must be split
;; with the new four triangles.  It also legalizes edges, and returns
;; the updated triangulation.
;; And a little bit of cheating: as I am not happy with tests for
;; collinearity, what I do is keeping always an edge and the distance
;; from the input point to the segment.  I loop over all edges then, 
;; and will substitute edge and distance if the new edge is encroached
;; upon by the input point, and the distance is smaller than previous.
;; As "distance", I use the "collinearity test" that offers the segment
;; functional.
;;
;;                ______v2______         ______v2______
;;               |      /\      |       |      /\      |
;;               | pT1 /..\ pT2 |       | pT1 /  \ pT2 |
;;               |    / .. \    |       |    /    \    |
;;               |   /  ..  \   |       |   /      \   |
;;               |  /   ..   \  |       |  /        \  |
;;               | /    .. s  \ |       | /   oT1    \ |
;;               |/ nT1 ..     \|       |/            \|
;;            v1 +======@@======+ v3 v1 +====segment===+ v3
;;               |\ nT2 p.     /|       |\            /|
;;               | \    .. s+1/ |       | \   oT2    / |
;;               |  \   ..   /  |       |  \        /  |
;;               |   \  ..  /   |       |   \      /   |
;;               |    \ .. /    |       |    \    /    |
;;               | pT4 \../ pT3 |       | pT4 \  / pT3 |
;;               |______\/______|       |______\/______|
;;                      v4                     v4
;;
(define (split-edge point-entry polygon-vertices triangulation)
  (let ((size (vector-length triangulation)))
    (let loop ((triangle-entry 0) 
	       (edge-vertex-1 0) (edge-vertex-2 0) 
	       (triangle-1 -1) (triangle-2 -1)
	       (threshold-collinearity 0.1))
      (if (= triangle-entry (vector-length triangulation))
	(if (= triangle-1 no-triangle)
	  triangulation
	  (let* ((vertex-1 edge-vertex-1)
		 (vertex-2 
		   (opposite-vertex-to-edge 
		     edge-vertex-1 
		     edge-vertex-2 
		     (vector-ref triangulation triangle-1)))
		 (vertex-3 edge-vertex-2)
		 (vertex-4 
		   (opposite-vertex-to-edge 
		     edge-vertex-1 
		     edge-vertex-2 
		     (vector-ref triangulation triangle-2)))
		 (pointer-1 
		   (pointer-associated-to-edge
		     vertex-1 
		     vertex-2 
		     (vector-ref triangulation triangle-1)))
		 (pointer-2
		   (pointer-associated-to-edge
		     vertex-2 
		     vertex-3 
		     (vector-ref triangulation triangle-1)))
		 (pointer-3 
		   (pointer-associated-to-edge
		     vertex-3 
		     vertex-4 
		     (vector-ref triangulation triangle-2)))
		 (pointer-4
		   (pointer-associated-to-edge
		     vertex-4 
		     vertex-1 
		     (vector-ref triangulation triangle-2)))
		 (new-triangle-1
		   (make-triangle vertex-1 vertex-2 point-entry
				  pointer-1 size triangle-2 0))
		 (new-triangle-2
		   (make-triangle vertex-1 point-entry vertex-4
				  triangle-1 (+ size 1) pointer-4 0))
		 (new-triangle-3
		   (make-triangle vertex-4 vertex-3 point-entry
				  pointer-3 size triangle-2 0))
		 (new-triangle-4
		   (make-triangle vertex-3 vertex-2 point-entry
				  pointer-2 triangle-1 (+ size 1) 0)))
	    (vector-set! triangulation triangle-1 new-triangle-1)
	    (vector-set! triangulation triangle-2 new-triangle-2)
	    (change-pointers 
	      pointer-2 triangle-1 size triangulation)
	    (change-pointers 
	      pointer-3 triangle-2 (+ size 1) triangulation)
	    (let ((new-triangulation
		    (list->vector
		      (append
			(vector->list triangulation)
			(list new-triangle-4 new-triangle-3)))))
	      (legalize-edge 
		vertex-1 vertex-2 triangle-1 
		new-triangulation polygon-vertices)
	      (legalize-edge
		vertex-2 vertex-3 size 
		new-triangulation polygon-vertices)
	      (legalize-edge
		vertex-3 vertex-4 (+ size 1) 
		new-triangulation polygon-vertices)
	      (legalize-edge
		vertex-4 vertex-1 triangle-2 
		new-triangulation polygon-vertices)
	      new-triangulation)))
	  
	  (let* ((triangle (vector-ref triangulation triangle-entry))
		 (point (vector-ref polygon-vertices point-entry))
		 (x (point-x point))
		 (y (point-y point))
		 (vertex-1 (triangle-vertex-1 triangle))
		 (vertex-2 (triangle-vertex-2 triangle))
		 (vertex-3 (triangle-vertex-3 triangle))
		 (functional-1 
		   (segment-functional 
		     (vector-ref polygon-vertices vertex-1)
		     (vector-ref polygon-vertices vertex-2)))
		 (collinearity-1 (abs (functional-1 x y)))
		 (encroached?-1 
		   (segment-encroached-upon-by? 
		     point 
		     (vector-ref polygon-vertices vertex-1)
		     (vector-ref polygon-vertices vertex-2)))
		 (functional-2 
		   (segment-functional 
		     (vector-ref polygon-vertices vertex-2)
		     (vector-ref polygon-vertices vertex-3)))
		 (collinearity-2 (abs (functional-2 x y)))
		 (encroached?-2
		   (segment-encroached-upon-by? 
		     point 
		     (vector-ref polygon-vertices vertex-2)
		     (vector-ref polygon-vertices vertex-3)))
		 (functional-3 
		   (segment-functional 
		     (vector-ref polygon-vertices vertex-3)
		     (vector-ref polygon-vertices vertex-1)))
		 (collinearity-3 (abs (functional-3 x y)))
		 (encroached?-3
		   (segment-encroached-upon-by? 
		     point 
		     (vector-ref polygon-vertices vertex-3)
		     (vector-ref polygon-vertices vertex-1)))
		 (test-collinearity 
		   (min collinearity-1 collinearity-2 collinearity-3))
		 (test-encroached-upon
		   (or encroached?-1 encroached?-2 encroached?-3)))
	    (if (and (<= test-collinearity threshold-collinearity)
		     test-encroached-upon)
	      (cond ((and (= collinearity-1 test-collinearity)
			  encroached?-1
			  (not (= (triangle-pointer-1-2 triangle)
				  no-triangle)))
		     (loop (+ 1 triangle-entry)
			   vertex-1 vertex-2
			   triangle-entry 
			   (triangle-pointer-1-2 triangle)
			   test-collinearity))
		    ((and (= collinearity-2 test-collinearity)
			  encroached?-2
			  (not (= (triangle-pointer-2-3 triangle)
				  no-triangle)))
		     (loop (+ 1 triangle-entry)
			   vertex-2 vertex-3
			   triangle-entry 
			   (triangle-pointer-2-3 triangle)
			   test-collinearity))
		    ((and (= collinearity-3 test-collinearity)
			  encroached?-3
			  (not (= (triangle-pointer-3-1 triangle)
				  no-triangle)))
		     (loop (+ 1 triangle-entry)
			   vertex-3 vertex-1
			   triangle-entry 
			   (triangle-pointer-3-1 triangle)
			   test-collinearity))
		    (else (loop (+ 1 triangle-entry)
				edge-vertex-1 edge-vertex-2
				triangle-1 triangle-2
				threshold-collinearity)))
	      (loop (+ 1 triangle-entry)
		    edge-vertex-1 edge-vertex-2
		    triangle-1 triangle-2
		    threshold-collinearity)))))))

;; CHANGE-POINTERS... no need to explain
;;
(define (change-pointers 
	  triangle-entry old-pointer new-pointer triangulation)
  (if (not (= triangle-entry no-triangle))
    (let ((triangle (vector-ref triangulation triangle-entry)))
      (cond ((= old-pointer (triangle-pointer-1-2 triangle))
	     (triangle-pointer-1-2-set! triangle new-pointer))
	    ((= old-pointer (triangle-pointer-2-3 triangle))
	     (triangle-pointer-2-3-set! triangle new-pointer))
	    ((= old-pointer (triangle-pointer-3-1 triangle))
	     (triangle-pointer-3-1-set! triangle new-pointer))))))

;; LEGALIZE-EDGE is recursive: it checks if the input edge is legal; if
;; it is not, then it changes stuff, and checks if new edges are legal.
;; But fear not, it has been proved that this recursion must end.
;; Notice there are two cases:
;; Case 1: edge originated from splitting of an interior triangle or
;;         the edge of an interior triangle.  We change edge by segment
;;         [v2-v4], update all triangles, and check legality in the new
;;         triangles.
;;                    
;;                ____v2____  
;;               |    /\    |
;;               |pT1/  \pT2|
;;               |  /    \  |
;;               | /   <------------- triangle-entry
;;               |/        \|
;;            v1 +===edge===+ v3
;;               |\        /|
;;               | \      / |
;;               |  \  <------------- opposite-ptriangle
;;               |pT4\  /pT3|
;;               |____\/____|
;;                    v4       
;;         
;;         If edge is not legal, we change to this:
;;
;;                ____v2____
;;               |    /\    |
;;               |pT1/||\pT2|
;;               |  / || \  |
;;               | /  ||  \ |                             
;;               |/ <---------------- triangle-entry
;;            v1 +    ||    + v3
;;               |\   || <----------- opposite-ptriangle
;;               | \  ||  / |
;;               |  \ || /  |                   
;;               |pT4\||/pT3|
;;               |____\/____|
;;                    v4
;;
;;         And then check if [v1-v4] is legal in pT4.
;;                        if [v3-v4] is legal in pT3.
;;
;; Case 2: edge originated from splitting of an exterior triangle or 
;;         the edge of an exterior triangle.  This is legal.  Return.
;;
;;                ____v2____
;;               |    /\    |
;;               |pT1/  \pT2|
;;               |  /    \  |
;;               | /      \ |
;;               |/        \|
;;            v1 +===edge===+ v3
;;                    || 
;;                    \/
;;                no-triangle
;;
(define (legalize-edge 
	  edge-vertex-1 edge-vertex-2 
	  triangle-entry 
	  triangulation polygon-vertices)
  (if (not (= triangle-entry no-triangle))
    (let ((opposite-ptriangle 
	    (pointer-associated-to-edge
	      edge-vertex-1 
	      edge-vertex-2 
	      (vector-ref triangulation triangle-entry))))
      (if (not (= opposite-ptriangle no-triangle))
	(if (= (triangle-status 
		 (vector-ref triangulation triangle-entry))
	       (triangle-status
		 (vector-ref triangulation opposite-ptriangle)))
	  (let* ((triangle (vector-ref triangulation triangle-entry))
		 (opposite-triangle
		   (vector-ref triangulation opposite-ptriangle))
		 (vertex-1 edge-vertex-1)
		 (vertex-2 
		   (opposite-vertex-to-edge 
		     edge-vertex-1 edge-vertex-2 
		     triangle))
		 (vertex-3 edge-vertex-2)
		 (vertex-4
		   (opposite-vertex-to-edge
		     edge-vertex-1 edge-vertex-2 
		     opposite-triangle))
		 (pointer-1
		   (pointer-associated-to-edge
		     vertex-1 vertex-2 
		     triangle))
		 (pointer-2
		   (pointer-associated-to-edge
		     vertex-2 vertex-3
		     triangle))
		 (pointer-3
		   (pointer-associated-to-edge
		     vertex-3 vertex-4
		     opposite-triangle))
		 (pointer-4
		   (pointer-associated-to-edge
		     vertex-4 vertex-1
		     opposite-triangle)))
	    (if (point-inside-circumcircle? 
		  (vector-ref polygon-vertices vertex-4)
		  (vector-ref polygon-vertices vertex-1)
		  (vector-ref polygon-vertices vertex-2)
		  (vector-ref polygon-vertices vertex-3))
	      (begin
		(vector-set! triangulation triangle-entry
			     (make-triangle 
			       vertex-1 vertex-2 vertex-4 
			       pointer-1 opposite-ptriangle pointer-4 
			       (triangle-status triangle)))
		(vector-set! triangulation opposite-ptriangle 
			     (make-triangle 
			       vertex-2 vertex-3 vertex-4 
			       pointer-2 pointer-3 triangle-entry
			       (triangle-status opposite-triangle)))
		(change-pointers 
		  pointer-2 triangle-entry opposite-ptriangle 
		  triangulation)
		(change-pointers 
		  pointer-4 opposite-ptriangle triangle-entry
		  triangulation)
		(legalize-edge 
		  vertex-1 vertex-4 pointer-4 
		  triangulation polygon-vertices)
		(legalize-edge 
		  vertex-3 vertex-4 pointer-3 
		  triangulation polygon-vertices)))))))))

;; The following two functions, OPPOSITE-VERTEX and POINTER-ASSOCIATED
;; are used in both split-edge and legalize-edge
;;
(define (opposite-vertex-to-edge vertex-1 vertex-2 triangle)
  (let ((point-1 (triangle-vertex-1 triangle))
	(point-2 (triangle-vertex-2 triangle))
	(point-3 (triangle-vertex-3 triangle)))
    (cond ((or (and (= point-1 vertex-1)
		    (= point-2 vertex-2))
	       (and (= point-1 vertex-2)
		    (= point-2 vertex-1)))
	   point-3)
	  ((or (and (= point-2 vertex-1)
		    (= point-3 vertex-2))
	       (and (= point-2 vertex-2)
		    (= point-3 vertex-1)))
	   point-1)
	  ((or (and (= point-3 vertex-1)
		    (= point-1 vertex-2))
	       (and (= point-3 vertex-2)
		    (= point-1 vertex-1)))
	   point-2))))
  
(define (pointer-associated-to-edge vertex-1 vertex-2 triangle)
  (let ((point-1 (triangle-vertex-1 triangle))
	(point-2 (triangle-vertex-2 triangle))
	(point-3 (triangle-vertex-3 triangle)))
    (cond ((or (and (= point-1 vertex-1)
		    (= point-2 vertex-2))
	       (and (= point-1 vertex-2)
		    (= point-2 vertex-1)))
	   (triangle-pointer-1-2 triangle))
	  ((or (and (= point-2 vertex-1)
		    (= point-3 vertex-2))
	       (and (= point-2 vertex-2)
		    (= point-3 vertex-1)))
	   (triangle-pointer-2-3 triangle))
	  ((or (and (= point-3 vertex-1)
		    (= point-1 vertex-2))
	       (and (= point-3 vertex-2)
		    (= point-1 vertex-1)))
	   (triangle-pointer-3-1 triangle))
	   (else no-triangle))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Second Step: Inclusion of the boundary edges.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Now we sweep over the points in the polygon-vertices, and focus on
;; the boundary-pointer field.  If that pointer is not -1, then we check
;; if the segment so defined is an edge in the triangulation.
;; If it is not there, we introduce the midpoint of that segment,
;; updating where necessary.
;;
(define (include-edges-in-triangulation polygon)
  (let ((vertices  (polygon-vertices polygon))
	(segments  (polygon-segments polygon))
	(holes     (polygon-holes polygon))
	(triangulation (polygon-triangulation polygon)))
    (let loop ((counter 0) 
	       (vertices vertices) 
	       (segments segments)
	       (triangulation triangulation))
      (if (= counter (vector-length segments))
	(make-polygon vertices segments holes triangulation)
	(let ((from-vertex (segment-from (vector-ref segments counter)))
	      (to-vertex (segment-to (vector-ref segments counter))))
	  (if (not (is-segment-in-triangulation?
		     from-vertex to-vertex triangulation))
	    (let* ((size-vertices (vector-length vertices))
		   (size-segments (vector-length segments))
		   (new-vertices 
		     (vector-append 
		       vertices 
		       (segment-midpoint 
			 (vector-ref vertices from-vertex)
			 (vector-ref vertices to-vertex))))
		   (new-segments
		     (vector-append 
		       segments 
		       (make-segment size-vertices to-vertex #t)))
		   (new-triangulation 
		     (introduce-point-in-triangulation
		       size-vertices new-vertices triangulation)))
	      (segment-to-set! (vector-ref new-segments counter)
			       size-vertices)
	      (loop 0 new-vertices new-segments new-triangulation))
	    (loop (+ counter 1) vertices segments triangulation)))))))
		
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Third Step: Get rid of exterior triangles
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; GET-RID-OF-EXTERIOR-TRIANGLES: we use the holes information of the 
;; polygon structure.  For each point-in-hole, we locate the triangle 
;; or edge it belongs to, and eat adjacent triangles until we find 
;; either an interior or a triangle previously eaten.
;; WARNING!!!!  This procedure is time-consuming: For each exterior
;; triangle it checks three new, and start the recursion on them.
;;
(define (get-rid-of-exterior-triangles polygon)
  (let ((vertices (polygon-vertices polygon))
	(segments (polygon-segments polygon))
	(holes    (polygon-holes polygon))
	(triangulation (polygon-triangulation polygon)))
    (let loop ((hole-entry 0))
      (if (= hole-entry (vector-length holes))
	polygon
	(begin
	  (eat-adjacent-triangles 
	    (vector-ref holes hole-entry)
	    vertices segments triangulation)
	  (loop (+ hole-entry 1)))))))

;; EAT-ADJACENT-TRIANGLES-TO-TRIANGLE: When we know that a triangle
;; is exterior, we can perform this operation of eating it and eating
;; all adjacent triangles to it until a border is reached.
;;
(define (eat-adjacent-triangles-to-triangle 
	  triangle-entry triangulation segments)
  (if (and (not (= triangle-entry no-triangle))
	   (interior-triangle? 
	     (vector-ref triangulation triangle-entry)))
    (let* ((triangle (vector-ref triangulation triangle-entry))
	   (adjacent-1 (triangle-pointer-1-2 triangle))
	   (border?-1 
	     (is-intersection-boundary-edge? 
	       triangle-entry adjacent-1 triangulation segments))
	   (adjacent-2 (triangle-pointer-2-3 triangle))
	   (border?-2 
	     (is-intersection-boundary-edge? 
	       triangle-entry adjacent-2 triangulation segments))
	   (adjacent-3 (triangle-pointer-3-1 triangle))
	   (border?-3 
	     (is-intersection-boundary-edge? 
	       triangle-entry adjacent-3 triangulation segments)))
      (triangle-status-set! triangle 1)
      (if (not border?-1)
	(eat-adjacent-triangles-to-triangle 
	  adjacent-1 triangulation segments))
      (if (not border?-2)
	(eat-adjacent-triangles-to-triangle
	  adjacent-2 triangulation segments))
      (if (not border?-3)
	(eat-adjacent-triangles-to-triangle
	  adjacent-3 triangulation segments)))))

;; EAT-ADJACENT-TRIANGLES-TO-EDGE:  As in (split-edge), this function
;; "knows" the input point belongs to an edge, but does not know which
;; one.  So it looks for it in the same collinearity-&-encroaching way,
;; and when it has it, picks one of the two adjacent triangles and
;; performs (eat-adjacent-triangles-to-triangle)
;;
(define (eat-adjacent-triangles-to-edge 
	  point triangulation vertices segments)
  (let loop ((triangle-entry 0) 
	     (triangle-to-be-eaten no-triangle)
	     (threshold-collinearity 0.1))
    (if (= triangle-entry (vector-length triangulation))
      (eat-adjacent-triangles-to-triangle
	triangle-to-be-eaten triangulation segments)
      (let* ((triangle (vector-ref triangulation triangle-entry))
	     (x (point-x point))
	     (y (point-y point))
	     (vertex-1 (triangle-vertex-1 triangle))
	     (vertex-2 (triangle-vertex-2 triangle))
	     (vertex-3 (triangle-vertex-3 triangle))
	     (functional-1 
	       (segment-functional 
		 (vector-ref vertices vertex-1)
		 (vector-ref vertices vertex-2)))
	     (collinearity-1 (abs (functional-1 x y)))
	     (encroached?-1 
	       (segment-encroached-upon-by? 
		 point 
		 (vector-ref vertices vertex-1)
		 (vector-ref vertices vertex-2)))
	     (functional-2 
	       (segment-functional 
		 (vector-ref vertices vertex-2)
		 (vector-ref vertices vertex-3)))
	     (collinearity-2 (abs (functional-2 x y)))
	     (encroached?-2
	       (segment-encroached-upon-by? 
		 point 
		 (vector-ref vertices vertex-2)
		 (vector-ref vertices vertex-3)))
	     (functional-3 
	       (segment-functional 
		 (vector-ref vertices vertex-3)
		 (vector-ref vertices vertex-1)))
	     (collinearity-3 (abs (functional-3 x y)))
	     (encroached?-3
	       (segment-encroached-upon-by? 
		 point 
		 (vector-ref vertices vertex-3)
		 (vector-ref vertices vertex-1)))
	     (test-collinearity 
	       (min collinearity-1 collinearity-2 collinearity-3))
	     (test-encroached-upon
	       (or encroached?-1 encroached?-2 encroached?-3)))
	(if (and (< test-collinearity threshold-collinearity)
		 test-encroached-upon)
	  (cond ((and (= collinearity-1 test-collinearity)
		      encroached?-1)
		 (loop (+ 1 triangle-entry)
		       triangle-entry 
		       test-collinearity))
		((and (= collinearity-2 test-collinearity)
		      encroached?-2)
		 (loop (+ 1 triangle-entry)
		       triangle-entry 
		       test-collinearity))
		((and (= collinearity-3 test-collinearity)
		      encroached?-3)
		 (loop (+ 1 triangle-entry)
		       triangle-entry
		       test-collinearity))
		(else (loop (+ 1 triangle-entry)
			    triangle-to-be-eaten
			    threshold-collinearity)))
	  (loop (+ 1 triangle-entry)
		triangle-to-be-eaten
		threshold-collinearity))))))
      
;; EAT-ADJACENT-TRIANGLES: Takes as input a point we know exterior to
;; the polygon, locates either edge or triangle it belongs to, and
;; performs eating operations until borders are reached.
;;
(define (eat-adjacent-triangles point vertices segments triangulation)
  (let loop ((triangle-entry 0))
    (if (= triangle-entry (vector-length triangulation))
      (eat-adjacent-triangles-to-edge 
	point triangulation vertices segments)
      (let* ((triangle (vector-ref triangulation triangle-entry))
	     (vertex-1 
	       (vector-ref vertices (triangle-vertex-1 triangle)))
	     (vertex-2
	       (vector-ref vertices (triangle-vertex-2 triangle)))
	     (vertex-3
	       (vector-ref vertices (triangle-vertex-3 triangle))))
	(if (point-inside-triangle?
	      point vertex-1 vertex-2 vertex-3)
	  (eat-adjacent-triangles-to-triangle
	    triangle-entry triangulation segments)
	  (loop (+ 1 triangle-entry)))))))
		  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Fourth Step: Handling angles
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; GET-RID-OF-ANGLES: It sweeps over the triangulation vector and deci-
;; des if each triangle has a right to be improved.  A triangle with all
;; its angles greater than the threshold angle or an exterior triangle
;; to the polygon (including holes as exterior) are examples of lowest
;; class: they cannot induce a modification in the triangulation.  
;; For any other triangle (interior with some angle less than threshold)
;; we might modify or not. Let us examine the possible cases:
;; (1) Triangle with at most one boundary edge: 
;;              induce modification.
;; (2) Corner triangle (two boundary edges) with corner angle less than
;;     threshold and the other angles greater: 
;;              do not induce modification.
;; (3) Corner triangle with two angles less than threshold, one of them
;;     being the corner angle: 
;;              induce modification.
;; (4) Notice that it is possible to have a corner triangle with corner
;;     angle greater than threshold and the other two less, and still 
;;     being legal in the sense of Delaunay (no other vertex lies inside
;;     its circumcircle).  But these are extreme cases (for instance, 
;;     Brad's triangulation of a square): 
;;              do not induce modification.
;; 
;; Just one warning; it does not matter if the input polygon has its   
;; exterior triangles marked.  This code first set them all interior
;; and then marks the exterior ones.
;;
;(define (get-rid-of-angles-less-than-threshold 
;	  threshold-angle polygon)
;  
;  (define (helper triangle-entry polygon)
;
;    (define (encroaching-step triangle polygon position)
;      (let* ((vertices (polygon-vertices polygon))
;	     (segments (polygon-segments polygon))
;	     (holes (polygon-holes polygon))
;	     (triangulation (polygon-triangulation polygon))
;	     (point (triangle-circumcenter triangle vertices))
;	     (bad-segment 
;	       (first-segment-encroached-upon-by 
;		 point segments vertices)))
;	(if (= bad-segment -1) 
;	  (let* ((vertices-size (vector-length vertices))
;		 (step-vertices (vector-append vertices point))
;		 (step-triangulation 
;		   (introduce-point-in-triangulation
;		     vertices-size step-vertices triangulation))
;		 (new-polygon
;		   (include-edges-in-triangulation
;		     (make-polygon
;		       step-vertices 
;		       segments 
;		       holes 
;		       step-triangulation))))
;	      (helper 0 new-polygon))
;	    (helper 0 (split-this-edge bad-segment polygon)))))
;
;    (let ((vertices (polygon-vertices polygon))
;	  (segments (polygon-segments polygon))
;	  (holes (polygon-holes polygon))
;	  (triangulation (polygon-triangulation polygon)))
;      (if (= triangle-entry (vector-length triangulation))
;	polygon
;	(let* ((triangle (vector-ref triangulation triangle-entry))
;	       (pvertex-1 (triangle-vertex-1 triangle))
;	       (vertex-1 (vector-ref vertices pvertex-1))
;	       (pvertex-2 (triangle-vertex-2 triangle))
;	       (vertex-2 (vector-ref vertices pvertex-2))
;	       (pvertex-3 (triangle-vertex-3 triangle))
;	       (vertex-3 (vector-ref vertices pvertex-3))
;	       (angle-1 (wedge-angle vertex-3 vertex-1 vertex-2))
;	       (angle-2 (wedge-angle vertex-1 vertex-2 vertex-3))
;	       (angle-3 (wedge-angle vertex-2 vertex-3 vertex-1))
;	       (minimum-angle (min angle-1 angle-2 angle-3))
;	       (maximum-angle (max angle-1 angle-2 angle-3)))
;	  (if (>= minimum-angle threshold-angle)
;	    (helper (+ 1 triangle-entry) polygon)
;	    (if (border-triangle? triangle segments)
;	      (let ((corner 
;		      (corner-angle triangle vertices segments)))
;		(if (< corner threshold-angle)
;		  (if (>= (- 180 threshold-angle)
;			  (+ corner maximum-angle))
;		    (helper (+ 1 triangle-entry) polygon)
;		    (encroaching-step
;		      triangle polygon (+ triangle-entry 1)))
;		  (encroaching-step 
;		    triangle polygon (+ triangle-entry 1))))
;	      (encroaching-step 
;		triangle polygon (+ triangle-entry 1))))))))
;
;  (helper 0 polygon))
;
;(define (split-this-edge segment-entry polygon)
;  (let* ((vertices (polygon-vertices polygon))
;	 (segments (polygon-segments polygon))
;	 (holes (polygon-holes polygon))
;	 (triangulation (polygon-triangulation polygon))
;	 (segment (vector-ref segments segment-entry))
;	 (vertices-size (vector-length vertices))
;	 (triangulation-size (vector-length triangulation))
;	 (pvertex-1 (segment-from segment))
;	 (vertex-1 (vector-ref vertices pvertex-1))
;	 (pvertex-3 (segment-to segment))
;	 (vertex-3 (vector-ref vertices pvertex-3))
;	 (step-vertices 
;	   (vector-append vertices 
;			  (segment-midpoint vertex-1 vertex-3)))
;	 (step-segments
;	   (vector-append 
;	     segments
;	     (make-segment vertices-size (segment-to segment) #t)))
;	 (segment-info (triangles-per-edge segment triangulation))
;	 (old-ptriangle-1 (car segment-info))
;	 (old-triangle-1 (vector-ref triangulation old-ptriangle-1))
;	 (old-ptriangle-2 (cadr segment-info))
;	 (old-triangle-2 (vector-ref triangulation old-ptriangle-2))
;	 (pvertex-2 (opposite-vertex-to-edge 
;		      pvertex-1 pvertex-3 old-triangle-1))
;	 (pvertex-4 (opposite-vertex-to-edge
;		      pvertex-1 pvertex-3 old-triangle-2))
;	 (pointer-1 (pointer-associated-to-edge
;		      pvertex-1 pvertex-2 old-triangle-1))
;	 (pointer-2 (pointer-associated-to-edge
;		      pvertex-2 pvertex-3 old-triangle-1))
;	 (pointer-3 (pointer-associated-to-edge
;		      pvertex-3 pvertex-4 old-triangle-2))
;	 (pointer-4 (pointer-associated-to-edge
;		      pvertex-4 pvertex-1 old-triangle-2))
;	 (new-triangle-1
;	   (make-triangle
;	     pvertex-1 pvertex-2 vertices-size
;	     pointer-1 triangulation-size old-ptriangle-2
;	     (triangle-status old-triangle-1)))
;	 (new-triangle-2
;	   (make-triangle
;	     pvertex-2 pvertex-3 vertices-size
;	     pointer-2 (+ triangulation-size 1) old-ptriangle-1
;	     (triangle-status old-triangle-1)))
;	 (new-triangle-3
;	   (make-triangle
;	     pvertex-3 pvertex-4 vertices-size
;	     pointer-3 old-ptriangle-2 triangulation-size
;	     (triangle-status old-triangle-2)))
;	 (new-triangle-4
;	   (make-triangle
;	     pvertex-4 pvertex-1 vertices-size
;	     pointer-4 old-ptriangle-1 (+ triangulation-size 1)
;	     (triangle-status old-triangle-2)))
;	 (step-triangulation
;	   (vector-append
;	     (vector-append triangulation
;			    new-triangle-2)
;	     new-triangle-3)))
;    (segment-to-set! segment vertices-size)
;    (vector-set! step-triangulation
;		 old-ptriangle-1
;		 new-triangle-1)
;    (vector-set! step-triangulation
;		 old-ptriangle-2
;		 new-triangle-4)
;    (change-pointers 
;      pointer-2 old-ptriangle-1 
;      triangulation-size step-triangulation)
;    (change-pointers
;      pointer-3 old-ptriangle-2 
;      (+ 1 triangulation-size) step-triangulation)
;    (legalize-edge pvertex-1 pvertex-2 old-ptriangle-1
;		   step-triangulation step-vertices)
;    (legalize-edge pvertex-2 pvertex-3 triangulation-size
;		   step-triangulation step-vertices)
;    (legalize-edge pvertex-3 pvertex-4 (+ triangulation-size 1)
;		   step-triangulation step-vertices)
;    (legalize-edge pvertex-4 pvertex-1 old-ptriangle-2
;		   step-triangulation step-vertices)
;    (let* ((step-polygon 
;	     (make-polygon 
;	       step-vertices step-segments 
;	       holes step-triangulation))) 
;      (include-edges-in-triangulation step-polygon))))
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Another Fourth Step
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; It is assumed that the code knows at all times which triangles are 
;; exterior.  So it happens that after some insertion of points, the 
;; original delaunay triangulation loses track of that status for all
;; new triangles.  A possible solution is not to legalize edges after
;; new insertions.  This makes the code faster in some sense, but per-
;; haps time-consuming in the long run.
;;
(define (get-rid-of-angles-less-than-threshold threshold-angle polygon)

  (define (helper triangle-entry polygon)

    (define (encroaching-step triangle polygon position)

      (define (split-this-edge segment-entry polygon position)
	(let* ((vertices (polygon-vertices polygon))
	       (segments (polygon-segments polygon))
	       (holes (polygon-holes polygon))
	       (triangulation (polygon-triangulation polygon))
	       (vertices-size (vector-length vertices))
	       (triangulation-size (vector-length triangulation))
	       (segment (vector-ref segments segment-entry))
	       (triangles (triangles-per-edge segment triangulation))
	       (old-ptriangle-1 (car triangles))
	       (old-ptriangle-2 (cadr triangles)))
	  (if (or (negative? old-ptriangle-1)
		  (negative? old-ptriangle-2))
	    (helper position polygon)
	    (let* ((pvertex-1 (segment-from segment))
		   (pvertex-3 (segment-to segment))
		   (vertex-1 (vector-ref vertices pvertex-1))
		   (vertex-3 (vector-ref vertices pvertex-3))
		   (point (segment-midpoint vertex-1 vertex-3))
		   (old-triangle-1 
		     (vector-ref triangulation old-ptriangle-1))
		   (old-triangle-2 
		     (vector-ref triangulation old-ptriangle-2))
		   (pvertex-2 (opposite-vertex-to-edge 
				pvertex-1 pvertex-3 old-triangle-1))
		   (pvertex-4 (opposite-vertex-to-edge
				pvertex-1 pvertex-3 old-triangle-2))
		   (pointer-1 (pointer-associated-to-edge
				pvertex-1 pvertex-2 old-triangle-1))
		   (pointer-2 (pointer-associated-to-edge
				pvertex-2 pvertex-3 old-triangle-1))
		   (pointer-3 (pointer-associated-to-edge
				pvertex-3 pvertex-4 old-triangle-2))
		   (pointer-4 (pointer-associated-to-edge
				pvertex-4 pvertex-1 old-triangle-2))
		   (new-vertices (vector-append vertices point))
		   (new-segments 
		     (vector-append 
		       segments (make-segment 
				  vertices-size pvertex-3 #t)))
		   (new-triangle-1
		     (make-triangle pvertex-1 pvertex-2 vertices-size
				    pointer-1 triangulation-size 
				    old-ptriangle-2 
				    (triangle-status old-triangle-1)))
		   (new-triangle-2
		     (make-triangle pvertex-2 pvertex-3 vertices-size
				    pointer-2 (+ 1 triangulation-size)
				    old-ptriangle-1 
				    (triangle-status old-triangle-1)))
		   (new-triangle-3
		     (make-triangle pvertex-3 pvertex-4 vertices-size
				    pointer-3 old-ptriangle-2 
				    triangulation-size
				    (triangle-status old-triangle-2)))
		   (new-triangle-4
		     (make-triangle pvertex-4 pvertex-1 vertices-size
				    pointer-4 old-ptriangle-1
				    (+ 1 triangulation-size) 
				    (triangle-status old-triangle-2)))
		   (new-triangulation
		     (vector-append
		       (vector-append triangulation new-triangle-2)
		       new-triangle-3)))
	      (change-pointers 
		pointer-2 old-ptriangle-1 
		triangulation-size new-triangulation)
	      (change-pointers
		pointer-3 old-ptriangle-2 
		(+ 1 triangulation-size) new-triangulation)
	      (vector-set! 
		new-triangulation old-ptriangle-1 new-triangle-1)
	      (vector-set! 
		new-triangulation old-ptriangle-2 new-triangle-4)
	      (segment-to-set! segment vertices-size)
	      (legalize-edge pvertex-1 pvertex-2 old-ptriangle-1
			     new-triangulation new-vertices)
	      (legalize-edge pvertex-2 pvertex-3 triangulation-size
			     new-triangulation new-vertices)
	      (legalize-edge pvertex-3 pvertex-4 
			     (+ 1 triangulation-size)
			     new-triangulation new-vertices)
	      (legalize-edge pvertex-4 pvertex-1 old-ptriangle-2
			     new-triangulation new-vertices)
	      (helper 0
		(make-polygon 
		  new-vertices 
		  new-segments 
		  holes 
		  new-triangulation))))))
      
      (let* ((vertices (polygon-vertices polygon))
	     (segments (polygon-segments polygon))
	     (holes (polygon-holes polygon))
	     (triangulation (polygon-triangulation polygon))
	     (point (triangle-circumcenter triangle vertices))
	     (support (point-support point vertices triangulation))
	     (bad-segment 
	       (first-segment-encroached-upon-by 
		 point segments vertices)))
	(if (= bad-segment -1) 
	  (if (triangle? (car support))
	    (if (not (interior-triangle? (car support)))
	      (helper position polygon)
	      (helper 0 (split-this-triangle 
			point (cadr support) polygon)))
	    (if (not (interior-triangle? 
		       (vector-ref triangulation (cadr support))))
	      (helper position polygon)
	      (helper 0 (split-this-segment point support polygon))))
	  (split-this-edge bad-segment polygon position))))
    
      (let ((vertices (polygon-vertices polygon))
	    (segments (polygon-segments polygon))
	    (holes (polygon-holes polygon))
	    (triangulation (polygon-triangulation polygon)))
	(if (= triangle-entry (vector-length triangulation))
	  polygon
	  (let ((triangle (vector-ref triangulation triangle-entry)))
	    (if (not (interior-triangle? triangle))
	      (helper (+ 1 triangle-entry) polygon)
	      (let* ((triangle 
		       (vector-ref triangulation triangle-entry))
		     (pvertex-1 (triangle-vertex-1 triangle))
		     (vertex-1 (vector-ref vertices pvertex-1))
		     (pvertex-2 (triangle-vertex-2 triangle))
		     (vertex-2 (vector-ref vertices pvertex-2))
		     (pvertex-3 (triangle-vertex-3 triangle))
		     (vertex-3 (vector-ref vertices pvertex-3))
		     (angle-1 
		       (wedge-angle vertex-3 vertex-1 vertex-2))
		     (angle-2 
		       (wedge-angle vertex-1 vertex-2 vertex-3))
		     (angle-3 
		       (wedge-angle vertex-2 vertex-3 vertex-1))
		     (minimum-angle (min angle-1 angle-2 angle-3))
		     (maximum-angle (max angle-1 angle-2 angle-3)))
		(if (>= minimum-angle threshold-angle)
		  (helper (+ 1 triangle-entry) polygon)
		  (if (border-triangle? triangle segments)
		    (let ((corner 
			    (corner-angle 
			      triangle vertices segments)))
		      (if (< corner (* 2 threshold-angle))
			(if (>= (- 180 threshold-angle)
				(+ corner maximum-angle))
			  (helper (+ 1 triangle-entry) polygon)
			  (encroaching-step
			    triangle polygon (+ triangle-entry 1)))
			(encroaching-step 
			  triangle polygon (+ triangle-entry 1))))
		    (encroaching-step 
		      triangle polygon (+ triangle-entry 1))))))))))
(helper 0 (get-rid-of-exterior-triangles polygon)))

(define (split-this-triangle point triangle-entry polygon)
  (let* ((vertices (polygon-vertices polygon))
	 (segments (polygon-segments polygon))
	 (holes (polygon-holes polygon))
	 (triangulation (polygon-triangulation polygon))
	 (vertices-size (vector-length vertices))
	 (triangulation-size (vector-length triangulation))
	 (old-triangle (vector-ref triangulation triangle-entry))
	 (pvertex-1 (triangle-vertex-1 old-triangle))
	 (pvertex-2 (triangle-vertex-2 old-triangle))
	 (pvertex-3 (triangle-vertex-3 old-triangle))
	 (pointer-1 (triangle-pointer-1-2 old-triangle))
	 (pointer-2 (triangle-pointer-2-3 old-triangle))
	 (pointer-3 (triangle-pointer-3-1 old-triangle))
	 (new-triangle-1
	   (make-triangle pvertex-1 pvertex-2 vertices-size
			  pointer-1 triangulation-size
			  (+ 1 triangulation-size) 0))
	 (new-triangle-2
	   (make-triangle pvertex-2 pvertex-3 vertices-size
			  pointer-2 (+ 1 triangulation-size)
			  triangle-entry 0))
	 (new-triangle-3
	   (make-triangle pvertex-3 pvertex-1 vertices-size
			  pointer-3 triangle-entry
			  triangulation-size 0))
	 (new-vertices (vector-append vertices point))
	 (new-triangulation
	   (vector-append
	     (vector-append triangulation new-triangle-2)
	     new-triangle-3)))
    (change-pointers 
      pointer-2 triangle-entry 
      triangulation-size new-triangulation)
    (change-pointers
      pointer-3 triangle-entry 
      (+ 1 triangulation-size) new-triangulation)
    (vector-set! new-triangulation triangle-entry new-triangle-1)
    (legalize-edge pvertex-1 pvertex-2 triangle-entry
		   new-triangulation new-vertices)
    (legalize-edge pvertex-2 pvertex-3 triangulation-size
		   new-triangulation new-vertices)
    (legalize-edge pvertex-3 pvertex-1 (+ 1 triangulation-size)
		   new-triangulation new-vertices)
    (make-polygon new-vertices segments holes new-triangulation)))


(define (split-this-segment point support polygon)
  (let* ((vertices (polygon-vertices polygon))
	 (segments (polygon-segments polygon))
	 (holes (polygon-holes polygon))
	 (triangulation (polygon-triangulation polygon))
	 (vertices-size (vector-length vertices))
	 (triangulation-size (vector-length triangulation))
	 (old-ptriangle-1 (cadr support))
	 (old-ptriangle-2 (caddr support))
	 (old-triangle-1 (vector-ref triangulation old-ptriangle-1))
	 (old-triangle-2 (vector-ref triangulation old-ptriangle-2))
	 (pvertex-1 (segment-from (car support)))
	 (pvertex-3 (segment-to (car support)))
	 (pvertex-2 (opposite-vertex-to-edge 
		      pvertex-1 pvertex-3 old-triangle-1))
	 (pvertex-4 (opposite-vertex-to-edge
			pvertex-1 pvertex-3 old-triangle-2))
	   (pointer-1 (pointer-associated-to-edge
			pvertex-1 pvertex-2 old-triangle-1))
	   (pointer-2 (pointer-associated-to-edge
			pvertex-2 pvertex-3 old-triangle-1))
	   (pointer-3 (pointer-associated-to-edge
			pvertex-3 pvertex-4 old-triangle-2))
	   (pointer-4 (pointer-associated-to-edge
			pvertex-4 pvertex-1 old-triangle-2))
	   (new-vertices (vector-append vertices point))
	   (new-triangle-1
	     (make-triangle pvertex-1 pvertex-2 vertices-size
			    pointer-1 triangulation-size 
			    old-ptriangle-2 
			    (triangle-status old-triangle-1)))
	   (new-triangle-2
	     (make-triangle pvertex-2 pvertex-3 vertices-size
			    pointer-2 (+ 1 triangulation-size)
			    old-ptriangle-1 
			    (triangle-status old-triangle-1)))
	   (new-triangle-3
	     (make-triangle pvertex-3 pvertex-4 vertices-size
			    pointer-3 old-ptriangle-2 
			    triangulation-size
			    (triangle-status old-triangle-2)))
	   (new-triangle-4
	     (make-triangle pvertex-4 pvertex-1 vertices-size
			    pointer-4 old-ptriangle-1
			    (+ 1 triangulation-size) 
			    (triangle-status old-triangle-2)))
	   (new-triangulation
	     (vector-append
	       (vector-append triangulation new-triangle-2)
	       new-triangle-3)))
      (change-pointers 
	pointer-2 old-ptriangle-1 
	triangulation-size new-triangulation)
      (change-pointers
	pointer-3 old-ptriangle-2 
	(+ 1 triangulation-size) new-triangulation)
      (vector-set! new-triangulation old-ptriangle-1 new-triangle-1)
      (vector-set! new-triangulation old-ptriangle-2 new-triangle-4)
      (legalize-edge pvertex-1 pvertex-2 old-ptriangle-1
		     new-triangulation new-vertices)
      (legalize-edge pvertex-2 pvertex-3 triangulation-size
		     new-triangulation new-vertices)
      (legalize-edge pvertex-3 pvertex-4 (+ 1 triangulation-size)
		     new-triangulation new-vertices)
      (legalize-edge pvertex-4 pvertex-1 old-ptriangle-2
		     new-triangulation new-vertices)
      (make-polygon new-vertices segments holes new-triangulation)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Last Step: All together now
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (delaunay-triangulation threshold-angle poly-file gnu-file)
  (polygon->file 
    (time 
      (get-rid-of-angles-less-than-threshold
        threshold-angle
        (file->polygon poly-file)))
    gnu-file))
