(declare (standard-bindings) (extended-bindings) (block) (not safe))

(define-macro (FLOAT . rest) 
	      `(let ()
		 (declare (flonum))
		 ,@rest))

(define-macro (FIX . rest)
	      `(let ()
		 (declare (fixnum))
		 ,@rest))

(define-macro (GENERIC . rest)
	      `(let ()
		 (declare (generic))
		 ,@rest))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(define (read-file file-name)

  (define (read-points port points-size)
    (let ((result (make-vector points-size)))
      (do ((counter 0 (+ counter 1)))
	((= counter points-size) result)
	(let ((index        (- (read port) 1))
	      (x-coordinate (exact->inexact (read port)))
	      (y-coordinate (exact->inexact (read port))))
	  (vector-set! result (make-point index x-coordinate y-coordinate))))))

  (define (read-segments port segments-size points)
    (let ((result (make-vector segments-size)))
      (do ((counter 0 (+ counter 1)))
	((= counter segments-size) result)
	(let ((index      (- (read port) 1))
	      (from-index (- (read port) 1))
	      (to-index   (- (read port) 1)))
	  (vector-set! result index 
		       (make-segment (vector-ref points from-index)
				     (vector-ref points to-index)))))))

  (let* ((port (open-input-file file-name))
	 (vertices-size (read port))
	 (ignore (read port))
	 (ignore (read port))
	 (ignore (read port))
	 (vertices (read-points port vertices-size))
	 (segments-size (read port))
	 (ignore (read port))
	 (segments (read-segments port segments-size vertices))
	 (holes-size (read port))
	 (holes (read-points port holes-size)))
    (list vertices segments holes)))

(define (vertices-bounding-limits vertices)
  (let ((size (vector-length vertices)))
    (let loop ((index (- (vector-length vertices) 1)) 
	       (min-x +inf.) (max-x -inf.)
	       (min-y +inf.) (max-y -inf.))
      (if (< index 0)
	(list min-x max-x min-y max-y)
	(let* ((point (vector-ref vertices index))
	       (x (point-x point))
	       (y (point-y point)))
	  (loop (- index 1)
		(min min-x x) (max max-x x)
		(min min-y y) (max max-y y)))))))

(define (augment-vertices vertices)
  (let* ((limits (vertices-bounding-limits vertices))
	 (min-x (car limits))
	 (max-x (cadr limits))
	 (min-y (caddr limits))
	 (max-y (cadddr limits))
	 (NE (make-point (- (* 2.0 max-x) min-x)
			 (- (* 2.0 max-y) min-y)))
	 (SE (make-point (- (* 2.0 max-x) min-x)
			 (- (* 2.0 min-y) max-y)))
	 (SW (make-point (- (* 2.0 min-x) max-x)
			 (- (* 2.0 min-y) max-y)))
	 (NW (make-point (- (* 2.0 min-x) max-x)
			 (- (* 2.0 max-y) min-y)))
	 (new-vertices (make-Dvector vertices (vector-length vertices))))
    (Dvector-append new-vertices NE)
    (Dvector-append new-vertices SE)
    (Dvector-append new-vertices SW)
    (Dvector-append new-vertices NW)))

(define (naive-triangulation vertices)
  (
    
    
	      

(define (triangulate file triangulation-method removal-method)

  (define (add-bounding-box vertices holes method)
    (let ((vertices-size (Dvector-length vertices)))
      (let loop ((index 0) (xmin 0) (xmax 0) (ymin 0) (ymax 0))
	(if (= index vertices-size)
	  (begin
	    (Dvector-append vertices 
			    (make-point (- (* 2 xmax) xmin)
					(- (* 2 ymax) ymin) 
					(build-Dvector 0 no-point)
					(build-Dvector 0 no-triangle)))
	    (Dvector-append vertices 
			    (make-point (- (* 2 xmax) xmin)
					(- (* 2 ymin) ymax) 
					(build-Dvector 0 no-point)
					(build-Dvector 0 no-triangle)))
	    (Dvector-append vertices 
			    (make-point (- (* 2 xmin) xmax)
					(- (* 2 ymin) ymax) 
					(build-Dvector 0 no-point)
					(build-Dvector 0 no-triangle)))
	    (Dvector-append vertices 
			    (make-point (- (* 2 xmin) xmax)
					(- (* 2 ymax) ymin) 
					(build-Dvector 0 no-point)
					(build-Dvector 0 no-triangle)))
	    (make-polygon (list xmin xmax ymin ymax)
			  vertices
			  holes 
			  (removal-method 
			    (triangulation-method vertices) holes)))
	  (let ((x (point-x (Dvector-ref vertices index)))
		(y (point-y (Dvector-ref vertices index))))
	    (loop (+ index 1) 
		  (min xmin x) (max xmax x)
		  (min ymin y) (max ymax y)))))))

  (let* ((port (open-input-file file))
	 (vertices-size (read port))
	 (vertices (build-Dvector vertices-size)))
    (read port)(read port)(read port)
    (do ((vertex 0 (+ vertex 1)))
      ((= vertex vertices-size)

       (let ((segments-size (read port)))
	 (read port)
	 (do ((segment 0 (+ segment 1)))
	   ((= segment segments-size)
	    
	    (let* ((holes-size (read port))
		   (holes (make-vector holes-size)))
	      (do ((hole 0 (+ hole 1)))
		((= hole holes-size)
		 (add-bounding-box 
		   vertices holes triangulation-method))
		(vector-set! holes (- (read port) 1)
			     (make-point (read port)
					 (read port)
					 '()
					 (build-Dvector 0 no-triangle))))))

	   (read port)
	   (let* ((point-from (Dvector-ref vertices (- (read port) 1)))
		  (point-from-pointers (point-pointers point-from))
		  (point-to (Dvector-ref vertices (- (read port) 1))))
	     (Dvector-append (point-pointers point-from) point-to)))))

      (Dvector-set! vertices (- (read port) 1)
		    (make-point (read port) (read port) 
				(build-Dvector 0 no-point)
				(build-Dvector 0 no-triangle))))))

(define (delaunay-triangulation vertices)
  (let* ((vertices-size (Dvector-length vertices))
	 (triangle-1 
	   (make-triangle (Dvector-ref vertices (- vertices-size 4))
			  (Dvector-ref vertices (- vertices-size 3))
			  (Dvector-ref vertices (- vertices-size 2))
			  no-triangle no-triangle no-triangle))
	 (triangle-2
	   (make-triangle (Dvector-ref vertices (- vertices-size 2))
			  (Dvector-ref vertices (- vertices-size 1))
			  (Dvector-ref vertices (- vertices-size 4))
			  no-triangle no-triangle triangle-1))
	 (triangulation (build-Dvector 2)))
    (triangle-triangle-3-1-set! triangle-1 triangle-2)
    (Dvector-append 
      (point-triangles (Dvector-ref vertices (- vertices-size 4)))
      triangle-1)
    (Dvector-append 
      (point-triangles (Dvector-ref vertices (- vertices-size 4)))
      triangle-2)
    (Dvector-append 
      (point-triangles (Dvector-ref vertices (- vertices-size 3)))
      triangle-1)
    (Dvector-append 
      (point-triangles (Dvector-ref vertices (- vertices-size 2)))
      triangle-1)
    (Dvector-append 
      (point-triangles (Dvector-ref vertices (- vertices-size 2)))
      triangle-2)
    (Dvector-append 
      (point-triangles (Dvector-ref vertices (- vertices-size 1)))
      triangle-2)
    (Dvector-set! triangulation 0 triangle-1)
    (Dvector-set! triangulation 1 triangle-2)
    (do ((index 0 (+ index 1)))
      ((= index (- vertices-size 4))
       triangulation)
      (insert-point-in-triangulation 
	(Dvector-ref vertices index) triangulation))))

(define (insert-point-in-triangulation point triangulation)
  (let ((size (Dvector-length triangulation)))
    (let loop ((index 0))
      (if (= index size)
	(error "input point not in domain")
	(let* ((triangle (Dvector-ref triangulation index))
	       (vertex-1 (triangle-vertex-1 triangle))
	       (vertex-2 (triangle-vertex-2 triangle))
	       (vertex-3 (triangle-vertex-3 triangle)))
	  (cond ((point-inside-triangle? point vertex-1 vertex-2 vertex-3)
		 (split-this-triangle point triangle triangulation))
		((point-inside-segment? point vertex-1 vertex-2)
		 (split-this-segment point vertex-1 vertex-2
				     triangle 
				     (triangle-associated-to-segment
				       vertex-1 vertex-2 triangle)
				     triangulation))
		((point-inside-segment? point vertex-2 vertex-3)
		 (split-this-segment point vertex-2 vertex-3
				     triangle
				     (triangle-associated-to-segment
				       vertex-2 vertex-3 triangle)
				     triangulation))
		((point-inside-segment? point vertex-3 vertex-1)
		 (split-this-segment point vertex-3 vertex-1
				     triangle
				     (triangle-associated-to-segment
				       vertex-3 vertex-1 triangle)
				     triangulation))
		(else (loop (+ index 1)))))))))

(define (split-this-triangle point triangle triangulation)
  (let* ((vertex-1 (triangle-vertex-1 triangle))
         (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (triangle-2-3 (triangle-triangle-2-3 triangle))
	 (triangle-3-1 (triangle-triangle-3-1 triangle))
	 (new-triangle-1
	   (make-triangle vertex-2 vertex-3 point
			  triangle-2-3 no-triangle triangle))
	 (new-triangle-2
	   (make-triangle vertex-3 vertex-1 point
			  triangle-3-1 triangle new-triangle-1)))
    (Dvector-append (point-triangles point) triangle)
    (Dvector-append (point-triangles point) new-triangle-1)
    (Dvector-append (point-triangles point) new-triangle-2)
    (Dvector-append (point-triangles vertex-1) new-triangle-2)
    (Dvector-append (point-triangles vertex-2) new-triangle-1)
    (Dvector-append (point-triangles vertex-3) new-triangle-1)
    (change-pointers vertex-3 triangle new-triangle-2)
    (change-pointers triangle-2-3 triangle new-triangle-1)
    (change-pointers triangle-3-1 triangle new-triangle-2)
    (triangle-vertex-3-set! triangle point)
    (triangle-triangle-2-3-set! triangle new-triangle-1)
    (triangle-triangle-3-1-set! triangle new-triangle-2)
    (triangle-triangle-2-3-set! new-triangle-1 new-triangle-2)
    (Dvector-append triangulation new-triangle-1)
    (Dvector-append triangulation new-triangle-2)
    (legalize-edge vertex-1 vertex-2 triangle)
    (legalize-edge vertex-2 vertex-3 new-triangle-1)
    (legalize-edge vertex-3 vertex-1 new-triangle-2)))

(define (split-this-segment point vertex-1 vertex-3 
			    old-triangle-1 old-triangle-2
			    triangulation)
  (let* ((vertex-2 (opposite-vertex-to-segment 
		     vertex-1 vertex-3 old-triangle-1))
	 (vertex-4 (opposite-vertex-to-segment
		     vertex-1 vertex-3 old-triangle-2))
	 (triangle-1 (triangle-associated-to-segment
		       vertex-1 vertex-2 old-triangle-1))
	 (triangle-2 (triangle-associated-to-segment
		       vertex-2 vertex-3 old-triangle-1))
	 (triangle-3 (triangle-associated-to-segment
		       vertex-3 vertex-4 old-triangle-2))
	 (triangle-4 (triangle-associated-to-segment
		       vertex-4 vertex-1 old-triangle-2))
	 (new-triangle-1
	   (make-triangle vertex-2 vertex-3 point
			  triangle-2 no-triangle old-triangle-1))
	 (new-triangle-2
	   (make-triangle vertex-3 vertex-4 point
			  triangle-3 old-triangle-2 new-triangle-1)))
    (Dvector-append (point-triangles point) old-triangle-1)
    (Dvector-append (point-triangles point) old-triangle-2)
    (Dvector-append (point-triangles point) new-triangle-1)
    (Dvector-append (point-triangles point) new-triangle-2)
    (Dvector-append (point-triangles vertex-2) new-triangle-1)
    (change-pointers vertex-3 old-triangle-1 new-triangle-1)
    (change-pointers vertex-3 old-triangle-2 new-triangle-2)
    (Dvector-append (point-triangles vertex-4) new-triangle-2)
    (change-pointers triangle-2 old-triangle-1 new-triangle-1)
    (change-pointers triangle-3 old-triangle-2 new-triangle-2)
    (triangle-vertex-1-set! old-triangle-1 vertex-1)
    (triangle-vertex-2-set! old-triangle-1 vertex-2)
    (triangle-vertex-3-set! old-triangle-1 point)
    (triangle-triangle-1-2-set! old-triangle-1 triangle-1) 
    (triangle-triangle-2-3-set! old-triangle-1 new-triangle-1)
    (triangle-triangle-3-1-set! old-triangle-1 old-triangle-2)
    (triangle-vertex-1-set! old-triangle-2 vertex-4)
    (triangle-vertex-2-set! old-triangle-2 vertex-1)
    (triangle-vertex-3-set! old-triangle-2 point)
    (triangle-triangle-1-2-set! old-triangle-2 triangle-4)
    (triangle-triangle-2-3-set! old-triangle-2 old-triangle-1)
    (triangle-triangle-3-1-set! old-triangle-2 new-triangle-2)
    (triangle-triangle-2-3-set! new-triangle-1 new-triangle-2)
    (Dvector-append triangulation new-triangle-1)
    (Dvector-append triangulation new-triangle-2)
    (legalize-edge vertex-1 vertex-2 old-triangle-1)
    (legalize-edge vertex-2 vertex-3 new-triangle-1)
    (legalize-edge vertex-3 vertex-4 new-triangle-2)
    (legalize-edge vertex-4 vertex-1 old-triangle-2)))

(define (legalize-edge vertex-1 vertex-3 triangle)
  (let ((opposite-triangle 
	  (triangle-associated-to-segment vertex-1 vertex-3 triangle)))
    (if (triangle? opposite-triangle)
      (let ((vertex-2 (opposite-vertex-to-segment
		        vertex-1 vertex-3 triangle))
	    (vertex-4 (opposite-vertex-to-segment
			vertex-1 vertex-3 opposite-triangle)))
	(if (point-inside-circumcircle? vertex-4 vertex-1 vertex-2 vertex-3)
	  (let ((triangle-1 (triangle-associated-to-segment
			      vertex-1 vertex-2 triangle))
		(triangle-2 (triangle-associated-to-segment
			      vertex-2 vertex-3 triangle))
		(triangle-3 (triangle-associated-to-segment
			      vertex-3 vertex-4 opposite-triangle))
		(triangle-4 (triangle-associated-to-segment
			      vertex-4 vertex-1 opposite-triangle)))
	    (change-pointers vertex-1 opposite-triangle no-triangle)
	    (Dvector-append (point-triangles vertex-2) opposite-triangle)
	    (change-pointers vertex-3 triangle no-triangle)
	    (Dvector-append (point-triangles vertex-4) triangle)
	    (change-pointers triangle-2 triangle opposite-triangle)
	    (change-pointers triangle-4 opposite-triangle triangle)
	    (triangle-vertex-1-set! triangle vertex-1)
	    (triangle-vertex-2-set! triangle vertex-2)
	    (triangle-vertex-3-set! triangle vertex-4)
	    (triangle-triangle-1-2-set! triangle triangle-1)
	    (triangle-triangle-2-3-set! triangle opposite-triangle)
	    (triangle-triangle-3-1-set! triangle triangle-4)
	    (triangle-vertex-1-set! opposite-triangle vertex-2)
	    (triangle-vertex-2-set! opposite-triangle vertex-3)
	    (triangle-vertex-3-set! opposite-triangle vertex-4)
	    (triangle-triangle-1-2-set! opposite-triangle triangle-2)
	    (triangle-triangle-2-3-set! opposite-triangle triangle-3)
	    (triangle-triangle-3-1-set! opposite-triangle triangle)
	    (legalize-edge vertex-1 vertex-4 triangle)
	    (legalize-edge vertex-3 vertex-4 opposite-triangle)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Constrained Steiner-Delaunay-Triangulation
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (constrained-steiner-delaunay-triangulation vertices)
  (let ((triangulation (time (delaunay-triangulation vertices))))
    (let outer-loop ((index 0))
      (if (= index (Dvector-length vertices))
	triangulation
	(let* ((vertex (Dvector-ref vertices index))
	       (segments (point-pointers vertex))
	       (segments-size (Dvector-length segments))
	       (triangles (point-triangles vertex))
	       (triangles-size (Dvector-length triangles)))
	  (let segments-loop ((opposite-index 0))
	    (if (= opposite-index segments-size)
	      (outer-loop (+ index 1))
	      (let ((opposite-vertex (Dvector-ref segments opposite-index)))
		(let triangles-loop ((triangle-index 0))
		  (if (= triangle-index triangles-size)
		    (let ((midpoint
			    (midpoint vertex opposite-vertex)))
		      (change-pointers vertex opposite-vertex midpoint)
		      (Dvector-append vertices midpoint)
		      (insert-point-in-triangulation midpoint triangulation)
		      (outer-loop 0))
		    (let ((triangle (Dvector-ref triangles triangle-index)))
		      (if (is-triangle-vertex? opposite-vertex triangle)
			(segments-loop (+ opposite-index 1))
			(triangles-loop (+ triangle-index 1))))))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Removal of exterior triangles
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (remove-exterior-triangles triangulation holes)

  (define (erase-triangle-from-triangulation triangle triangulation)
    (if (triangle? triangle)
      (let ((size (Dvector-length triangulation)))
	(let loop ((index 0))
	  (if (= index size)
	    (error "triangle not present in triangulation")
	    (if (equal-triangles? triangle (Dvector-ref triangulation index))
	      (Dvector-set! triangulation index no-triangle)
	      (loop (+ index 1))))))))

  (define (eat-adjacent-exterior-triangles triangle triangulation)
    (if (triangle? triangle)
      (let ((size (Dvector-length triangulation))
	    (vertex-1 (triangle-vertex-1 triangle))
	    (vertex-2 (triangle-vertex-2 triangle))
	    (vertex-3 (triangle-vertex-3 triangle))
	    (triangle-1 (triangle-triangle-1-2 triangle))
	    (triangle-2 (triangle-triangle-2-3 triangle))
	    (triangle-3 (triangle-triangle-3-1 triangle)))
	(if (or (triangle? triangle-1)
		(triangle? triangle-2)
		(triangle? triangle-3))
	  (begin
	    (change-pointers triangle-1 triangle no-triangle)
	    (change-pointers triangle-2 triangle no-triangle)
	    (change-pointers triangle-3 triangle no-triangle)
	    (erase-triangle-from-triangulation triangle triangulation)
	    (if (not (is-a-segment? vertex-1 vertex-2))
	      (eat-adjacent-exterior-triangles triangle-1 triangulation))
	    (if (not (is-a-segment? vertex-2 vertex-3))
	      (eat-adjacent-exterior-triangles triangle-2 triangulation))
	    (if (not (is-a-segment? vertex-3 vertex-1))
	      (eat-adjacent-exterior-triangles triangle-3 triangulation)))))))
	
  (trace eat-adjacent-exterior-triangles)
  (let loop ((index 0))
    (if (= index (vector-length holes))
      triangulation
      (let* ((hole (vector-ref holes index))
	     (triangle (triangle-support hole triangulation))
	     (triangle-1 (triangle-triangle-1-2 triangle))
	     (triangle-2 (triangle-triangle-2-3 triangle))
	     (triangle-3 (triangle-triangle-3-1 triangle)))
	(if (and (not (triangle? triangle-1))
		 (not (triangle? triangle-2))
		 (not (triangle? triangle-3)))
	  (loop (+ index 1))
	  (begin
	    (eat-adjacent-exterior-triangles triangle triangulation)
	    (loop (+ index 1))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Abbreviations
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define DT delaunay-triangulation)
(define CSDT constrained-steiner-delaunay-triangulation)
(define R remove-exterior-triangles)
(define (NR triangulation holes) triangulation)
