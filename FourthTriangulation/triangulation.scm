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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  .poly file parser
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (file->polygon file)

  (define (add-bounding-box vertices vertices-size segments holes)
    (declare (flonum))
    (let loop ((index 0) (xmax 0) (xmin 0) (ymax 0) (ymin 0))
      (if (= index vertices-size)
	(begin
	  (Dvector-append vertices
			 (make-point (- (* 2 xmax) xmin)
				     (- (* 2 ymax) ymin)))
	  (Dvector-append vertices
			 (make-point (- (* 2 xmax) xmin)
				     (- (* 2 ymin) ymax)))
	  (Dvector-append vertices
			 (make-point (- (* 2 xmin) xmax)
				     (- (* 2 ymin) ymax)))
	  (Dvector-append vertices
			 (make-point (- (* 2 xmin) xmax)
				     (- (* 2 ymax) ymin)))
	  (make-polygon (list xmin xmax ymin ymax)
			vertices
			segments holes
			(delaunay-triangulation vertices vertices-size)))
	(let ((x (point-x (Dvector-ref vertices index)))
	      (y (point-y (Dvector-ref vertices index))))
	  (loop (+ index` 1) 
		(max xmax x) (min xmin x)
		(max ymax y) (min ymin y))))))

  (let* ((port (open-input-file file))
	 (vertices-size (read port))
	 (vertices (build-Dvector vertices-size)))
    (read port)(read port)(read port)
    (do ((vertex-index 0 (+ 1 vertex-index)))
      ((= vertex-index vertices-size)
       
       (let* ((segments-size (read port))
	      (segments (build-Dvector segments-size)))
	 (read port)
	 (do ((segment-index 0 (+ 1 segment-index)))
	   ((= segment-index segments-size)
	    
	    (let* ((holes-size (read port))
		   (holes (make-vector holes-size)))
	      (do ((hole-index 0 (+ 1 hole-index)))
		((= hole-index holes-size)
		 (add-bounding-box vertices vertices-size segments holes))
		(vector-set! holes (- (read port) 1)
			     (make-point (read port) (read port))))))
	   
	   (Dvector-set! segments (- (read port) 1)
			 (make-segment (- (read port) 1)
				       (- (read port) 1))))))
      
      (Dvector-set! vertices (- (read port) 1)
		    (make-point (read port) (read port))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;;  First Step: Constrained Delaunay Triangulation of convex hull of input
;;  polygon.
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; DELAUNAY-TRIANGULATION computes the delaunay triangulation of the given 
;; Dvector of vertices.  It starts with a naive triangulation of the bounding
;; box, and inserts every point from the Dvector at a time, updating triangles
;; and legalizing edges where necessary.
;; 
(define (delaunay-triangulation vertices vertices-size)

  (define (naive-triangulation vertices-size)
    (make-Dvector
      (list->vector
	(list
	  (make-triangle 
	    vertices-size (+ 1 vertices-size) (+ 2 vertices-size)
	    -1 -1 1)
	  (make-triangle
	    (+ vertices-size 2) (+ vertices-size 3) vertices-size
	    -1 -1 0)))
      2))

  (let ((triangulation (naive-triangulation vertices-size)))
    (let loop ((index 0))
      (if (= index vertices-size)
	triangulation
	(begin
	  (insert-point-in-triangulation index vertices triangulation)
	  (loop (+ index 1)))))))

;; INSERT-POINT-IN-TRIANGULATION locates the triangle or edge where the input
;; point lies and calls either SPLIT-THIS-TRIANGLE or SPLIT-THIS-EDGE.
;; 
(define (insert-point-in-triangulation pvertex vertices triangulation)
  (let ((size (Dvector-length triangulation))
	(vertex (Dvector-ref vertices pvertex)))
    (let loop ((index 0))
      (if (= index size)
	(error "input point not in range")
	(let* ((triangle  (Dvector-ref triangulation index))
	       (pvertex-1 (triangle-pvertex-1 triangle))
	       (pvertex-2 (triangle-pvertex-2 triangle))
	       (pvertex-3 (triangle-pvertex-3 triangle))
	       (vertex-1  (Dvector-ref vertices pvertex-1))
	       (vertex-2  (Dvector-ref vertices pvertex-2))
	       (vertex-3  (Dvector-ref vertices pvertex-3)))
	  (cond ((point-inside-triangle? vertex vertex-1 vertex-2 vertex-3)
		 (split-this-triangle pvertex index vertices triangulation))
		((point-inside-segment? vertex vertex-1 vertex-2)
		 (split-this-edge pvertex pvertex-1 pvertex-2
				  index (triangle-pointer-1-2 triangle) 
				  vertices triangulation))
		((point-inside-segment? vertex vertex-2 vertex-3)
		 (split-this-edge pvertex pvertex-2 pvertex-2
				  index (triangle-pointer-2-3 triangle) 
				  vertices triangulation))
		((point-inside-segment? vertex vertex-3 vertex-1)
		 (split-this-edge pvertex pvertex-3 pvertex-1
				  index (triangle-pointer-3-1 triangle) 
				  vertices triangulation))
		(else (loop (+ 1 index)))))))))

;; SPLIT-THIS-TRIANGLE replaces triangle in position triangle-entry with one of
;; the three new triangles, and appends the other two at the end of the
;; triangulation vector.  It legalizes edges if necessary.
;;                          
;;               _ _ _ _ _ _v2 _ _ _ _ _ _
;;                          /\
;;              |          /..\          |
;;                        / .. \
;;              |   pt1  /  ..  \   pt2  |
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
;;                   +     pt3      + 
;;                      +        +
;;                         +  + 
;;                            
(define (split-this-triangle pvertex ptriangle vertices triangulation)
  (let* ((size (Dvector-length triangulation))
	 (triangle (Dvector-ref triangulation ptriangle))
	 (pvertex-1 (triangle-pvertex-1 triangle))
	 (pvertex-2 (triangle-pvertex-2 triangle))
	 (pvertex-3 (triangle-pvertex-3 triangle))
	 (pointer-1 (triangle-pointer-1-2 triangle))
	 (pointer-2 (triangle-pointer-2-3 triangle))
	 (pointer-3 (triangle-pointer-3-1 triangle))
	 (new-triangle-1
	   (make-triangle pvertex-1 pvertex-2 pvertex
			  pointer-1 size (+ size 1)))
	 (new-triangle-2
	   (make-triangle pvertex-2 pvertex-3 pvertex
			  pointer-2 (+ size 1) ptriangle))
	 (new-triangle-3
	   (make-triangle pvertex-3 pvertex-1 pvertex
			  pointer-3 ptriangle size)))
    (Dvector-set! triangulation ptriangle new-triangle-1)
    (Dvector-append triangulation new-triangle-2)
    (Dvector-append triangulation new-triangle-3)
    (change-pointers pointer-2 ptriangle size triangulation)
    (change-pointers pointer-3 ptriangle (+ size 1) triangulation)
    (legalize-edge pvertex-1 pvertex-2 ptriangle vertices triangulation)
    (legalize-edge pvertex-2 pvertex-3 size vertices triangulation)
    (legalize-edge pvertex-3 pvertex-1 (+ size 1) vertices triangulation)))


;; SPLIT-THIS-EDGE replaces the two adjacent triangles sharing the input edge
;; by two new triangles, and appends the other two at the end of the
;; triangulation vector.  As before, it legalizes edges if nesessary.
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
(define (split-this-edge pvertex pvertex-1 pvertex-3 
			 old-ptriangle-1 old-ptriangle-2
			 vertices triangulation)
  (let* ((size (Dvector-length triangulation))
	 (old-triangle-1 (Dvector-ref triangulation old-ptriangle-1))
	 (old-triangle-2 (Dvector-ref triangulation old-ptriangle-2))
	 (pvertex-2 
	   (opposite-vertex-to-segment pvertex-1 pvertex-3 old-triangle-1))
	 (pvertex-4
	   (opposite-vertex-to-segment pvertex-1 pvertex-3 old-triangle-2))
	 (pointer-1 
	   (pointer-associated-to-segment pvertex-1 pvertex-2 old-triangle-1))
	 (pointer-2 
	   (pointer-associated-to-segment pvertex-2 pvertex-3 old-triangle-1))
	 (pointer-3 
	   (pointer-associated-to-segment pvertex-3 pvertex-4 old-triangle-2))
	 (pointer-4 
	   (pointer-associated-to-segment pvertex-4 pvertex-1 old-triangle-2))
	 (new-triangle-1
	   (make-triangle pvertex-1 pvertex-2 pvertex
			  pointer-1 size old-ptriangle-2))
	 (new-triangle-2
	   (make-triangle pvertex-2 pvertex-3 pvertex
			  pointer-2 (+ size 1) old-ptriangle-1))
	 (new-triangle-3
	   (make-triangle pvertex-3 pvertex-4 pvertex
			  pointer-3 old-ptriangle-2 size))
	 (new-triangle-4
	   (make-triangle pvertex-4 pvertex-1 pvertex
			  pointer-4 old-ptriangle-1 (+ size 1))))
    (Dvector-set! triangulation old-ptriangle-1 new-triangle-1)
    (Dvector-set! triangulation old-ptriangle-2 new-triangle-4)
    (Dvector-append triangulation new-triangle-2)
    (Dvector-append triangulation new-triangle-3)
    (change-pointers pointer-2 old-ptriangle-1 size triangulation)
    (change-pointers pointer-3 old-ptriangle-2 (+ size 1) triangulation)
    (legalize-edge pvertex-1 pvertex-2 old-ptriangle-1 vertices triangulation)
    (legalize-edge pvertex-2 pvertex-3 size vertices triangulation)
    (legalize-edge pvertex-3 pvertex-4 (+ size 1) vertices triangulation)
    (legalize-edge pvertex-4 pvertex-1 old-ptriangle-2 vertices triangulation)))

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
;;                no-ptriangle
;;
(define (legalize-edge pvertex-1 pvertex-3 ptriangle vertices triangulation)
  (if (not (= ptriangle no-ptriangle))
    (let* ((triangle (Dvector-ref triangulation ptriangle))
	   (opposite-ptriangle 
	     (pointer-associated-to-segment pvertex-1 pvertex-3 triangle)))
      (if (not (= opposite-ptriangle no-ptriangle))
	(let* ((opposite-triangle 
		 (Dvector-ref triangulation opposite-ptriangle))
	       (pvertex-2 
		 (opposite-vertex-to-segment 
		   pvertex-1 pvertex-3 triangle))
	       (pvertex-4
		 (opposite-vertex-to-segment 
		   pvertex-1 pvertex-3 opposite-triangle))
	       (vertex-1 (Dvector-ref vertices pvertex-1))
	       (vertex-2 (Dvector-ref vertices pvertex-2))
	       (vertex-3 (Dvector-ref vertices pvertex-3))
	       (vertex-4 (Dvector-ref vertices pvertex-4)))
	  (if (point-inside-circumcircle? vertex-4 vertex-1 vertex-2 vertex-3)
	    (let* ((pointer-1 
		     (pointer-associated-to-segment 
		       pvertex-1 pvertex-2 triangle))
		   (pointer-2 
		     (pointer-associated-to-segment 
		       pvertex-2 pvertex-3 triangle))
		   (pointer-3 
		     (pointer-associated-to-segment 
		       pvertex-3 pvertex-4 opposite-triangle))
		   (pointer-4 
		     (pointer-associated-to-segment 
		       pvertex-4 pvertex-1 opposite-triangle))
		   (new-triangle-1
		     (make-triangle pvertex-1 pvertex-2 pvertex-4
				    pointer-1 opposite-ptriangle pointer-4))
		   (new-triangle-2
		     (make-triangle pvertex-2 pvertex-3 pvertex-4
				    pointer-2 pointer-3 ptriangle)))
	      (Dvector-set! triangulation ptriangle new-triangle-1)
	      (Dvector-set! triangulation opposite-ptriangle new-triangle-2)
	      (change-pointers 
		pointer-2 ptriangle opposite-ptriangle triangulation)
	      (change-pointers
		pointer-4 opposite-ptriangle ptriangle triangulation)
	      (legalize-edge 
		pvertex-1 pvertex-4 pointer-4 vertices triangulation)
	      (legalize-edge
		pvertex-3 pvertex-4 pointer-3 vertices triangulation))))))))
