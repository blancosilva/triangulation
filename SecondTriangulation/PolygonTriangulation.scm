;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                ;;
;;    Fisrt stage: Delaunay triangulation of input vertices.      ;;
;;                                                                ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; POLYGON-TRIANGULATION takes a polygon, creates a big box around it,
;; uses the first vertex of the polygon to create a naive triangulation
;; of that box, and then introduces the remaining points one at a time, 
;; splitting either triangle or edge (where it lies), and legalizing
;; all edges so that no vertex is inside the circumcircle of another
;; triangle.
(define (polygon-triangulation polygon)
  (let loop ((position 1) 
	     (triangulation (naive-triangulation polygon)))
    (if (= position (vector-length polygon))
      triangulation
      (loop (+ position 1) 
	    (introduce-point-in-triangulation
	      (vector-ref polygon position)
	      triangulation)))))

;; NAIVE-TRIANGULATION creates the initial triangulation.  It uses the
;; vertices of boundary-box and the first vertex of the polygon.
(define (naive-triangulation polygon)
  (let ((NE (rectangle-ne (boundary-box polygon)))
	(SE (rectangle-se (boundary-box polygon)))
	(SW (rectangle-sw (boundary-box polygon)))
	(NW (rectangle-nw (boundary-box polygon))))
    (list->vector
      (list
	(make-triangle NE SE (vector-ref polygon 0)
		       -1 1 3)
	(make-triangle SE SW (vector-ref polygon 0)
		       -1 2 0)
	(make-triangle SW NW (vector-ref polygon 0)
		       -1 3 1)
	(make-triangle NW NE (vector-ref polygon 0)
		       -1 0 2)))))

;; PTRIANGLE-CONTAINING-POINT finds the triangle in the triangulation in
;; which point lies in.  If no triangle, it returns a pointer to the 
;; exterior (-1); otherwise, a pointer to the entry in the vector of 
;; triangles (triangulation).
(define (ptriangle-containing-point point triangulation)
  (let loop ((position 0))
    (if (= position (vector-length triangulation))
      no-triangle
      (if (point-inside-triangle? point 
				  (vector-ref triangulation position))
	position
	(loop (+ position 1))))))

;; EDGE-CONTAINING-POINT finds the edge in the triangulation that
;; contains the given point.  In case it finds such edge, it returns a 
;; list with that edge, plus the pointers of the triangles that share   
;; that edge.
;; In case no edge is found, it returns the list with no-triangles and
;; no-segment.
(define (edge-containing-point point triangulation)
  (let loop ((position 0))
    (if (= position (vector-length triangulation))
      (list no-segment no-triangle no-triangle)
      (let* ((triangle (vector-ref triangulation position))
	     (segment-1 (make-segment (triangle-vertex-1 triangle)
				      (triangle-vertex-2 triangle)))
	     (segment-2 (make-segment (triangle-vertex-2 triangle)
				      (triangle-vertex-3 triangle)))
	     (segment-3 (make-segment (triangle-vertex-3 triangle)
				      (triangle-vertex-1 triangle))))
	(cond ((point-inside-segment? point segment-1)
	       (list segment-1
		     position 
		     (triangle-ptriangle-1-2 triangle)))
	      ((point-inside-segment? point segment-2)
	       (list segment-2 
		     position
		     (triangle-ptriangle-2-3 triangle)))
	      ((point-inside-segment? point segment-3)
	       (list segment-3 
		     position 
		     (triangle-ptriangle-3-1 triangle)))
	      (else (loop (+ position 1))))))))

; CLOSEST-SEGMENT is similar to the previous.  As it will be explained
; later on, sometimes we prefer looking for the closest segment to a
; given point, rather than checking for collinearity.
(define (closest-segment point triangulation)
  (let loop ((position 0) 
	     (result (list no-segment no-triangle no-triangle))
	     (threshold 0.1))
    (if (= position (vector-length triangulation))
      result
      (let* ((triangle (vector-ref triangulation position))
	     (segment-1
	       (make-segment (triangle-vertex-1 triangle)
			     (triangle-vertex-2 triangle)))
	     (distance-1
	       (distance point
			 (segment-midpoint segment-1)))
	     (segment-2
	       (make-segment (triangle-vertex-2 triangle)
			     (triangle-vertex-3 triangle)))
	     (distance-2
	       (distance point
			 (segment-midpoint segment-2)))
	     (segment-3
	       (make-segment (triangle-vertex-3 triangle)
			     (triangle-vertex-1 triangle)))
	     (distance-3
	       (distance point
			 (segment-midpoint segment-3)))
	     (min-distance (min distance-1 distance-2 distance-3)))
	(if (> min-distance threshold)
	  (loop (+ position 1) result threshold) 
	  (cond ((= min-distance distance-1)
		 (loop (+ position 1)
		       (list segment-1
			     position
			     (triangle-ptriangle-1-2 triangle))
		       min-distance))
		((= min-distance distance-2)
		 (loop (+ position 1)
		       (list segment-2
			     position
			     (triangle-ptriangle-2-3 triangle))
		       min-distance))
		((= min-distance distance-3)
		 (loop (+ position 1)
		       (list segment-3
			     position
			     (triangle-ptriangle-1-2 triangle))))))))))
			 
;; INTRODUCE-POINT-IN-TRIANGULATION is the main function.  It looks
;; for a triangle or an edge where the input point lies, and splits
;; the found object, updating ilegal edges if necessary.  This function
;; uses split-edge-by-point and split-triangle-by-point.
;; A little cheating: sometimes the point is clearly on a segment --no
;; problem.  Sometimes the point is clearly on a triangle --no problem
;; either.  But so it happens that in some cases there are points that
;; seem to be in both a triangle and a segment (in that case we assume
;; it is on a segment), or worse, the point will seem to be in neither
;; edge nor triangle.  In that case we assume again that the point is
;; on an edge, and perform a search on all edges of the triangulation.
;; We collect the ones encroached upon by that point, and will choose
;; the one that is closest.
(define (introduce-point-in-triangulation point triangulation)
  (let ((target-ptriangle 
	  (ptriangle-containing-point point triangulation))
	(target-edge
	  (car (edge-containing-point point triangulation))))
    (cond ((and (= target-ptriangle no-triangle)
	        (same-segments? target-edge no-segment)) 
	   (split-edge-by-point point
				(closest-segment point triangulation)
				triangulation))
	  ((same-segments? target-edge no-segment)
	   (split-triangle-by-point point 
				    target-ptriangle 
				    triangulation))
	  (else
	    (split-edge-by-point point 
				 (edge-containing-point 
				   point triangulation)
				 triangulation)))))

;; SPLIT-TRIANGLE replaces triangle in position ptriangle with one of 
;; the three new triangles, and appends the other two at the end of the
;; triangulation vector.  It legalizes edges if necessary, and returns
;; the updated triangulation.
;;                          
;;               _ _ _ _ _ _v2 _ _ _ _ _ _
;;                          /\
;;              |          /..\          |
;;                        / .. \
;;              |  pT12  /  ..  \   pT23 |
;;                      /   ..   \
;;              |      / pT ..  l \      |
;;                    /     .@.    \
;;              |    /    .. p ..   \    |
;;                  /   ..       ..  \ 
;;              |  /  ..           .. \  |
;;                / ..      l+1      ..\ 
;;              |/.                    .\|
;;           v1 *------------------------*v3
;;               +                      +
;;                  +      pT31      + 
;;                     +          +
;;                        +    + 
;;                            
;;
(define (split-triangle-by-point point ptriangle triangulation)
  (let* ((length (vector-length triangulation))
	 (old-triangle (vector-ref triangulation ptriangle))
	 (vertex-1 (triangle-vertex-1 old-triangle))
	 (vertex-2 (triangle-vertex-2 old-triangle))
	 (vertex-3 (triangle-vertex-3 old-triangle))
	 (pointer-1-2 (triangle-ptriangle-1-2 old-triangle))
	 (pointer-2-3 (triangle-ptriangle-2-3 old-triangle))
	 (pointer-3-1 (triangle-ptriangle-3-1 old-triangle))) 
    (let ((new-triangle-1 (make-triangle
			    vertex-1 vertex-2 point
			    pointer-1-2 length (+ length 1)))
	  (new-triangle-2 (make-triangle
			    vertex-2 vertex-3 point
			    pointer-2-3 (+ length 1) ptriangle))
	  (new-triangle-3 (make-triangle
			    vertex-3 vertex-1 point
			    pointer-3-1 ptriangle length)))
      (change-pointers-in-ptriangle 
	pointer-2-3 ptriangle length triangulation)
      (change-pointers-in-ptriangle 
	pointer-3-1 ptriangle (+ length 1) triangulation) 
      (vector-set! triangulation ptriangle new-triangle-1)
      (let ((new-triangulation 
	      (list->vector 
		(append 
		  (vector->list triangulation) 
		  (list new-triangle-2 new-triangle-3)))))
	(legalize-edge-with-respect-to-ptriangle 
	  (make-segment vertex-1 vertex-2) 
	  ptriangle 
	  new-triangulation) 
	(legalize-edge-with-respect-to-ptriangle 
	  (make-segment vertex-2 vertex-3) 
	  length 
	  new-triangulation) 
	(legalize-edge-with-respect-to-ptriangle 
	  (make-segment vertex-3 vertex-1) 
	  (+ length 1) 
	  new-triangulation)
	new-triangulation))))
	
;; SPLIT-EDGE-BY-POINT replaces the two triangles sharing the edge that
;; must be split with the new four triangles.  It also legalizes edges. 
;; It returns the updated triangulation.
;;
;;                ______v2______         ______v2______
;;               |      /\      |       |      /\      |
;;               | pT1 /..\ pT2 |       | pT1 /  \ pT2 |
;;               |    / .. \    |       |    /    \    |
;;               |   /  ..  \   |       |   /      \   |
;;               |  /   ..   \  |       |  /        \  |
;;               | /    .. l  \ |       | /   oT1    \ |
;;               |/ nT1 ..     \|       |/            \|
;;            v1 +======@@======+ v3 v1 +====segment===+ v3
;;               |\ nT2 p.     /|       |\            /|
;;               | \    .. l+1/ |       | \   oT2    / |
;;               |  \   ..   /  |       |  \        /  |
;;               |   \  ..  /   |       |   \      /   |
;;               |    \ .. /    |       |    \    /    |
;;               | pT4 \../ pT3 |       | pT4 \  / pT3 |
;;               |______\/______|       |______\/______|
;;                      v4                     v4
;;
(define (split-edge-by-point point edge-info triangulation)
  (let* ((length (vector-length triangulation))
	 (segment (car edge-info))
	 (vertex-1 (segment-point-1 segment))
	 (old-ptriangle-1 (cadr edge-info))
	 (old-ptriangle-2 (caddr edge-info))
	 (old-triangle-1 (vector-ref triangulation old-ptriangle-1))
	 (old-triangle-2 (vector-ref triangulation old-ptriangle-2))
	 (vertex-2 (opposite-vertex-to-segment
		     segment old-triangle-1))
	 (vertex-3 (segment-point-2 segment))
	 (vertex-4 (opposite-vertex-to-segment
		     segment old-triangle-2))
	 (pointer-1 (pointer-associated-to-segment
		      (make-segment vertex-1 vertex-2)
		      old-triangle-1))
	 (pointer-2 (pointer-associated-to-segment
		      (make-segment vertex-2 vertex-3)
		      old-triangle-1))
	 (pointer-3 (pointer-associated-to-segment
		      (make-segment vertex-3 vertex-4)
		      old-triangle-2))
	 (pointer-4 (pointer-associated-to-segment
		      (make-segment vertex-4 vertex-1)
		      old-triangle-2))
	 (new-triangle-1
	   (make-triangle vertex-1 vertex-2 point 
			  pointer-1 length old-ptriangle-2))
	 (new-triangle-2
	   (make-triangle vertex-1 vertex-4 point 
			  pointer-4 (+ length 1) old-ptriangle-1))
	 (new-triangle-3
	   (make-triangle vertex-4 vertex-3 point
			  pointer-3 length old-ptriangle-2))
	 (new-triangle-4
	   (make-triangle vertex-3 vertex-2 point
			  pointer-2 old-ptriangle-1 (+ length 1))))
    (vector-set! triangulation old-ptriangle-1 new-triangle-1)
    (vector-set! triangulation old-ptriangle-2 new-triangle-2)
    (change-pointers-in-ptriangle 
      pointer-2 old-ptriangle-1 length triangulation)
    (change-pointers-in-ptriangle
      pointer-3 old-ptriangle-2 (+ length 1) triangulation)
    (let ((new-triangulation 
	    (list->vector 
	      (append 
		(vector->list triangulation)
		(list new-triangle-4 new-triangle-3)))))
      (legalize-edge-with-respect-to-ptriangle
	(make-segment vertex-1 vertex-2)
	old-ptriangle-1
	new-triangulation)
      (legalize-edge-with-respect-to-ptriangle
	(make-segment vertex-2 vertex-3)
	length
	new-triangulation)
      (legalize-edge-with-respect-to-ptriangle
	(make-segment vertex-3 vertex-4)
	(+ length 1)
	new-triangulation)
      (legalize-edge-with-respect-to-ptriangle
	(make-segment vertex-4 vertex-1)
	old-ptriangle-2
	new-triangulation)
      new-triangulation)))
			  
;; LEGALIZE-EDGE is recursive: it checks if the input edge is legal; if
;; it is not, then it changes stuff, and checks if new edges are legal.
;; But fear not, it has been proved that this recursion must end.
;; Notice there are two cases:
;; Case 1: edge originated from splitting of an interior triangle or
;;         the edge of an interior triangle.  We change edge by segment
;;         [v2-v4], update all triangles, and check legality in the new
;;         triangles.
;;                    
;;              v5____v2____v6
;;               |    /\    |
;;               |pT1/  \pT2|
;;               |  /    \  |
;;               | /   <------------- ptriangle
;;               |/        \|
;;            v1 +===edge===+ v3
;;               |\        /|
;;               | \      / |
;;               |  \  <------------- opposite-ptriangle
;;               |pT4\  /pT3|
;;               |____\/____|
;;             v8     v4     v7
;;         
;;         If edge is not legal, we change to this:
;;
;;                ____v2____
;;               |    /\    |
;;               |pT1/||\pT2|
;;               |  / || \  |
;;               | /  ||  \ |                             
;;               |/ <---------------- ptriangle
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
(define (legalize-edge-with-respect-to-ptriangle 
	  segment ptriangle triangulation)
  (let ((opposite-ptriangle 
	  (adjacent-ptriangle-associated-to-segment
	    segment ptriangle triangulation)))
    (if (not (= opposite-ptriangle no-triangle))
      (let* ((triangle (vector-ref triangulation ptriangle)) 
	     (vertex-1 (segment-point-1 segment)) 
	     (vertex-2 (opposite-vertex-to-segment segment triangle)) 
	     (vertex-3 (segment-point-2 segment)) 
	     (opposite-triangle 
	       (vector-ref triangulation opposite-ptriangle)) 
	     (vertex-4 
	       (opposite-vertex-to-segment segment opposite-triangle)) 
	     (pointer-1 
	       (pointer-associated-to-segment 
		 (make-segment vertex-1 vertex-2) triangle)) 
	     (pointer-2 
	       (pointer-associated-to-segment 
		 (make-segment vertex-2 vertex-3) triangle)) 
	     (pointer-3 
	       (pointer-associated-to-segment 
		 (make-segment vertex-3 vertex-4) 
		 opposite-triangle)) 
	     (pointer-4 
	       (pointer-associated-to-segment 
		 (make-segment vertex-4 vertex-1) 
		 opposite-triangle)))
	(if (point-inside-circumcircle? vertex-4 triangle) 
	  (begin
	    (triangle-vertex-1-set! triangle vertex-1) 
	    (triangle-vertex-2-set! triangle vertex-2) 
	    (triangle-vertex-3-set! triangle vertex-4) 
	    (triangle-ptriangle-1-2-set! triangle pointer-1) 
	    (triangle-ptriangle-2-3-set! triangle opposite-ptriangle) 
	    (triangle-ptriangle-3-1-set! triangle pointer-4) 
	    (if (not (= pointer-2 no-triangle)) 
	      (change-pointers-in-ptriangle 
	        pointer-2 ptriangle opposite-ptriangle triangulation)) 
	    (if (not (= pointer-4 no-triangle)) 
	      (change-pointers-in-ptriangle 
	        pointer-4 opposite-ptriangle ptriangle triangulation)) 
	    (triangle-vertex-1-set! opposite-triangle vertex-2) 
	    (triangle-vertex-2-set! opposite-triangle vertex-3) 
	    (triangle-vertex-3-set! opposite-triangle vertex-4) 
	    (triangle-ptriangle-1-2-set! opposite-triangle pointer-2) 
	    (triangle-ptriangle-2-3-set! opposite-triangle pointer-3) 
	    (triangle-ptriangle-3-1-set! opposite-triangle ptriangle) 
	    (if (not (= pointer-3 no-triangle)) 
	      (legalize-edge-with-respect-to-ptriangle
	        (make-segment vertex-3 vertex-4)
	        pointer-3
	        triangulation))
	    (if (not (= pointer-4 no-triangle))
	      (legalize-edge-with-respect-to-ptriangle
	        (make-segment vertex-1 vertex-4)
	        pointer-4
	        triangulation))))))))
	      
;; ADJACENT-PTRIANGLE:  Given a segment and the pointer to a triangle
;; containing this segment in the triangulation, finds the pointer to
;; the adjacent triangle. If none, returns no-triangle.
(define (adjacent-ptriangle-associated-to-segment 
	  segment ptriangle triangulation)
  (let loop ((position 0) (pointer no-triangle))
    (cond ((= position (vector-length triangulation))
	   pointer)
	  ((= position ptriangle)
	   (loop (+ position 1) pointer))
	  (else
	    (let ((triangle (vector-ref triangulation position)))
	      (if (is-segment-a-triangle-edge?
		    segment triangle)
		(loop (+ position 1) position)
		(loop (+ position 1) pointer)))))))

;; CHANGE-POINTERS-IN-PTRIANGLE... do I need to explain this one?
(define (change-pointers-in-ptriangle 
	  ptriangle old-pointer new-pointer triangulation)
  (if (not (= ptriangle no-triangle))
    (let* ((triangle (vector-ref triangulation ptriangle))) 
      (cond ((= (triangle-ptriangle-1-2 triangle) old-pointer) 
	     (triangle-ptriangle-1-2-set! triangle new-pointer)) 
	    ((= (triangle-ptriangle-2-3 triangle) old-pointer) 
	     (triangle-ptriangle-2-3-set! triangle new-pointer)) 
	    ((= (triangle-ptriangle-3-1 triangle) old-pointer) 
	     (triangle-ptriangle-3-1-set! triangle new-pointer))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                ;;
;; Second stage: Inclusion of boundary edges in the triangulation ;;
;;                                                                ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; BRANCH takes a vector and creates a vector with one more entry
;; by changing one of the entries in the former by two new ones.
(define (branch items-vector 
		equality-function 
		element-to-branch
		new-element-1
		new-element-2)
  (let loop ((position 0) (result '()))
       (if (= position (vector-length items-vector))
	 (list->vector result)
	 (if (equality-function 
	       (vector-ref items-vector position)
	       element-to-branch)
	   (loop (+ position 1) 
		 (append result
			 (list new-element-1 new-element-2)))
	   (loop (+ position 1)
		 (append result 
			 (list (vector-ref items-vector position))))))))

(define (include-segments-in-triangulation super-polygon)
  (let loop ((position 0)
	     (vertices (super-polygon-vertices super-polygon))
	     (segments (super-polygon-segments super-polygon))
	     (triangulation 
	       (super-polygon-triangulation super-polygon)))
    (if (= position (vector-length segments))
      (make-super-polygon vertices segments triangulation)
      (let ((segment (vector-ref segments position))) 
	(if (is-segment-in-triangulation?  
	      (vector-ref segments position) 
	      triangulation) 
	  (loop (+ position 1) vertices segments triangulation) 
	  (loop 0
		(branch vertices
			same-points?
			(segment-point-1 segment)
			(segment-point-1 segment)
			(segment-midpoint segment))
	        (branch segments 
		        same-segments?
		        segment 
		        (make-segment (segment-point-1 segment)
				      (segment-midpoint segment))
		        (make-segment (segment-midpoint segment)
				      (segment-point-2 segment)))
	        (introduce-point-in-triangulation
	          (segment-midpoint segment)
	          triangulation)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                ;;
;; Third stage: Inclusion of angles in the triangulation          ;;
;;                                                                ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; DELAUNAY-TRIANGULATION is the only function the user should call.
;; Given a threshold angle and a polygon, it first gets a basic
;; Delaunay triangulation of the polygon without inserting any new
;; vertices.  It improves that triangulation by including the boundary
;; edges (inserting new vertices in the boundary edges if necessary),
;; and finally gets rid of all triangles with minimum angle less than
;; the threshold. 
(define (delaunay-triangulation threshold-angle polygon segments)
  (get-rid-of-angles-less-than-threshold
    threshold-angle
    (include-segments-in-triangulation
      (make-super-polygon polygon
			  segments
			  (polygon-triangulation polygon)))))
  
;; GET-RID assigns to each triangle a "social status".  The lowest class
;; triangles are either exterior to the polygon or with minimum angle 
;; greater than the threshold.  These do not have any right to order 
;; other triangles to be modified.   A middle class of triangles is 
;; constituted by corner interior triangles with their corner angle
;; lower than the threshold.  They will only order other triangles to
;; be modified if their minimum angle is less than the corner angle.
;; Finally, the first class of triangles is constituted by all other
;; interior angles such that their minimum angle is less than the
;; threshold.  These have the right (and the obligation) to order any
;; other triangle to split.
;; The splitting operations are as follows:
;; (1) A first or middle class triangle is found.
;; (2) Compute its circumcenter.
;; (3) Check if that circumcenter encroaches any boundary-edge. If that
;;     is the case, split the segment in half, update the triangulation
;;     and start all over with the new triangulation.  If no boundary-
;;     edge is encroached upon by this circumcenter, then introduce the
;;     point in the triangulation and start all over again with the new
;;     triangulation.
(define (get-rid-of-angles-less-than-threshold 
	  threshold-angle super-polygon)
  (define (helper position super-polygon)
    ;;;;;
    (define (encroaching-step triangle)
      (let* ((vertices (super-polygon-vertices super-polygon))
	     (segments (super-polygon-segments super-polygon))
	     (triangulation (super-polygon-triangulation
			      super-polygon))
	     (bad-segment
	       (first-segment-encroached-upon-by-triangle-circumcenter
		 triangle segments)))
	(if (same-segments? bad-segment no-segment)
	  (let* ((circumcenter (triangle-circumcenter triangle))
		 (new-vertices (branch vertices
				       same-points?
				       (vector-ref vertices 0)
				       (vector-ref vertices 0)
				       circumcenter))
		 (new-triangulation 
		   (introduce-point-in-triangulation
		     circumcenter triangulation))
		 (new-super-polygon
		   (make-super-polygon
		     new-vertices segments new-triangulation))
		 (newest-super-polygon
		   (include-segments-in-triangulation
		     new-super-polygon)))
	    (helper 0 newest-super-polygon))
	  (let* ((new-vertices (branch vertices
				       same-points?
				       (segment-point-1
					 bad-segment)
				       (segment-point-1
					 bad-segment)
				       (segment-midpoint
					 bad-segment)))
		 (new-segments (branch segments
				       same-segments?
				       bad-segment
				       (make-segment
					 (segment-point-1
					   bad-segment)
					 (segment-midpoint
					   bad-segment))
				       (make-segment
					 (segment-midpoint
					   bad-segment)
					 (segment-point-2
					   bad-segment))))
		 (new-triangulation 
		   (introduce-point-in-triangulation
		     (segment-midpoint bad-segment)
		     triangulation)))
	    (let* ((new-super-polygon 
		     (make-super-polygon 
		       new-vertices new-segments new-triangulation)) 
		   (newest-super-polygon 
		     (include-segments-in-triangulation 
		       new-super-polygon)))
	      (helper 0 newest-super-polygon))))))
    ;;;;;
    (if (= position 
	   (vector-length 
	     (super-polygon-triangulation super-polygon)))
      super-polygon
      (let* ((vertices (super-polygon-vertices super-polygon))
	     (segments (super-polygon-segments super-polygon))
	     (polygon  (segments->polygon segments))
	     (triangulation (super-polygon-triangulation 
			      super-polygon))
	     (triangle (vector-ref triangulation position))
	     (circumcenter (triangle-circumcenter triangle)))
	(cond ((or (>= (triangle-minimum-angle triangle)
		       threshold-angle)
		   (not (is-inner-triangle? triangle polygon))
		   (not (point? circumcenter)))
	       (helper (+ position 1) super-polygon))
	      ((is-corner-triangle? triangle segments)
	       (if (= (triangle-minimum-angle triangle)
		      (triangle-corner-angle-in-polygon
			triangle polygon))
		 (helper (+ position 1) super-polygon)
		 (encroaching-step triangle))) 
	      (else (encroaching-step triangle))))))
  (helper 0 super-polygon))
		 
;; FIRST-SEGMENT looks for a segment so that the circumcenter of the
;; input triangle is contained in the circle with that segment for
;; diameter.
(define (first-segment-encroached-upon-by-triangle-circumcenter
	  triangle segments)
  (let loop ((position 0))
    (if (= position (vector-length segments))
      no-segment
      (if (point-encroaches-segment?
	    (triangle-circumcenter triangle)
	    (vector-ref segments position))
	(vector-ref segments position)
	(loop (+ position 1))))))
