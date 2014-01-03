#|                                       XX        XX           
                                          X         X           
                                          X         X           
 XXXXXX  XXX XX   XXXXX   XXXX   XXX X    XXXXX     X     XXXXX 
  X    X   XX  X X     X      X   X X X   X    X    X    X     X
  X    X   X     XXXXXXX  XXXXX   X X X   X    X    X    XXXXXXX
  X    X   X     X       X    X   X X X   X    X    X    X      
  X    X   X     X     X X    X   X X X   X    X    X    X     X
  XXXXX  XXXXX    XXXXX   XXXX X XX X XX XXXXXX   XXXXX   XXXXX 
  X                                                             
 XXX                                                            |#

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

#| These are the basic objects: point, segment, edge and triangle.
A point is represented by its coordinates in the x-y plane.

A segment is given by the two points defining it.  An edge is a segment with
extra structure, since "edge" is always part of at least one and at most two
triangles.

In the triangle structure we will gather its three defining vertices, segments
(no edges) and adjacent triangles. |#

(define-structure point x y)

(define-structure segment from to)

(define-structure edge segment triangle-1 triangle-2)

(define-structure triangle 
		  vertex-1 vertex-2 vertex-3 
		  ;; segment-k and triangle-k are opposite to vertex-k
		  segment-1 segment-2 segment-3    
		  triangle-1 triangle-2 triangle-3)

(define (build-triangle vertex-1 vertex-2 vertex-3)
  (let ((segment-1 (make-segment vertex-2 vertex-3))
	(segment-2 (make-segment vertex-3 vertex-1))
	(segment-3 (make-segment vertex-1 vertex-2)))
    (make-triangle vertex-1 vertex-2 vertex-3
		   segment-1 segment-2 segment-3
		   #f #f #f)))

#| We also need to handle "collections".  We can change the definition of
collection later to something different, without having to change the rest of
the code.  We will have a collection-builder, collection-accessor and
collection-information.  Maybe something else later on. |#

(define (collection . elements)
  (let* ((size (length elements))
	 (set (build-Dvector size #f)))
    (do ((index 0 (+ index 1)))
      ((= index size) set)
      (Dvector-set! set index (list-ref elements index)))))

(define (add-elements-to-collection set . elements)
  (let ((size (length elements)))
    (do ((index 0 (+ index 1)))
      ((= index (- size 1)) (Dvector-append set (list-ref elements (- size 1))))
      (Dvector-append set (list-ref elements index)))))

(define (remove-element-from-collection set element eq-function)
  ;; this is a bit strange.  This function does not know what elements there
  ;; are in the collection, or if there are relations among them.  It will only
  ;; locate the input element and eliminate it (like Terminator).
  (if element
    (let loop ((index (- (Dvector-length set) 1)))
      (if (< index 0) 
	(error "input element not in collection")
	(if (eq-function element (Dvector-ref set index))
	  (Dvector-set! set index #f)
	  (loop (- index 1)))))))

(define collection-length Dvector-length)

(define (collection-ref set index) (Dvector-ref set index))
	 
#|         XX                               X            XX             
            X                                      X      X             
            X                                      X      X             
  XXXX      X     XXXXXX  XXXXX  XXX XX   XXX     XXXX    X XX   XXX X  
      X     X    X    X  X     X   XX  X    X      X      XX  X   X X X 
  XXXXX     X    X    X  X     X   X        X      X      X   X   X X X 
 X    X     X    X    X  X     X   X        X      X      X   X   X X X 
 X    X     X     XXXXX  X     X   X        X      X  X   X   X   X X X 
  XXXX X  XXXXX       X   XXXXX  XXXXX    XXXXX     XX   XXX XXX XX X XX
                      X                                                 
                  XXXX                                                  |#

#| First step: Read a poly file into a vector of points, a vector of segments
and a vector of "holes", the latter understood as a point for every exterior
region to the given PSLG. |# 

(define (read-file file-name)

  (define (read-points port points-size)
    (let ((result (make-vector points-size)))
      (do ((counter 0 (+ counter 1)))
	((= counter points-size) result)
	(let ((index        (- (read port) 1))
	      (x-coordinate (read port))
	      (y-coordinate (read port)))
	  (vector-set! result index (make-point x-coordinate y-coordinate))))))

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

#| Second step: we need to create a region big enough to contain every vertex
of the input PSLG and the vertices we might eventually create in the future.
We have chosen the following method:
  (1) Find the smallest rectangular region with sides parallel to the
      coordinate axis in which the PSLG is contained (vertices-bounding-limits)
  (2) Consider a rectangle centered in the previous, three times bigger.

With the vertices of that box, which we label NE, SE, SW and NW, we create a
triangulation.  In our case we went for {{NE,SE,SW}, {SW,NW,NE}}. |#

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

(define (initial-triangulation vertices)
  (let* ((limits (vertices-bounding-limits vertices))
	 (min-x (car limits))
	 (max-x (cadr limits))
	 (min-y (caddr limits))
	 (max-y (cadddr limits)))
    (let ((NE (make-point (- (* 2.0 max-x) min-x)
			  (- (* 2.0 max-y) min-y)))
	  (SE (make-point (- (* 2.0 max-x) min-x)
			  (- (* 2.0 min-y) max-y)))
	  (SW (make-point (- (* 2.0 min-x) max-x)
			  (- (* 2.0 min-y) max-y)))
	  (NW (make-point (- (* 2.0 min-x) max-x)
			  (- (* 2.0 max-y) min-y))))
      (let ((triangle-1 (build-triangle NE SE SW))
	    (triangle-2 (build-triangle SW NW NE)))
	(triangle-adjacent-triangles-set! triangle-1 #f triangle-2 #f)
	(triangle-adjacent-triangles-set! triangle-2 #f triangle-1 #f) 
	(collection triangle-1 triangle-2)))))

#| Third step: Delaunay Triangulation of the input vertices, without constrains
caused by the segments.  This gives us a triangulation of the convex hull of
the set of vertices.  In order to create such triangulation we handle the
vertices one at a time: Each vertex is either in a triangle or an edge; if it
is in a triangle, we "split" the triangle into three new triangles --or four,
if the input point lies in an edge of the triangle--, get rid of the previous,
and "legalize edges".  The functions used by delaunay-triangulation are in the
file delaunay.scm.  |#

(define (delaunay-triangulation poly-file)

  (define (insert-point-in-triangulation point triangulation)
    ;; INSERT-POINT-IN-TRIANGULATION:
    ;; (1) Given the input point, find first a triangle in triangulation such 
    ;;     that the point lies either in the interior or one of its edges.  
    ;; (2) If the point lies in the interior of the triangle, split that 
    ;;     triangle into three new, get rid of the older and include the new 
    ;;     ones in the triangulation.  
    ;; (3) If the point lies in one of the edges, we have to split two triangles
    ;;     into four.  Get rid of the old ones and include the new ones in the
    ;;     triangulation.
    ;; (4) Either way, we have to make all edges legal afterwards.
    (let ((object (point-support point triangulation)))
      (if (triangle? object)
	(split-triangle-by-point point object triangulation)
	(split-edge-by-point point object triangulation))))

  (let* ((parsed-file (read-file poly-file))
	 (vertices (car parsed-file))
	 (triangulation (initial-triangulation vertices)))
    (do ((index (- (vector-length vertices) 1) (- index 1)))
      ((< index 0) triangulation)
      (insert-point-in-triangulation 
	(vector-ref vertices index) triangulation))))
