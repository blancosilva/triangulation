#| SPLIT-TRIANGLE-BY-POINT: Given input point and triangle, split the latter in
three new triangles.  One of them will substitute the old, and the other two
will be included in the triangulation.   This is low level; I guess this method
is the better way to deal with this procedure.  It avoids looking in the
collection of triangles for the old triangle, and focuses in performing
manipulations on it and in its adjacents --which we gather from the triangle
itself-- and then include stuff in the triangulation, no matter where or how. 

After this operation, we must check for ilegal edges.  The candidates are the
three edges of the old triangle. 

The following picture summarizes the operations:

  old triangle (ot)       new triangles/segments


         v2                        v2
          /\                        /\
   t3<-s3/  \ s1->t1               /|<------ns2
        /  . \         nt1------->/ .. \<-------nt2
       /  p   \            ns1--->.'p '.<---ns3
     v1--------v3              v1--------v3
        s2->t2                      ^ 
                                    |
                                   nt3                                  |#

(define (split-triangle-by-point point triangle triangulation)
  (display "Splitting triangle\n")
  (let ((vertex-1 (triangle-vertex-1 triangle))
	(vertex-2 (triangle-vertex-2 triangle))
	(vertex-3 (triangle-vertex-3 triangle))
	(triangle-1 (triangle-triangle-1 triangle))
	(triangle-2 (triangle-triangle-2 triangle))
	(triangle-3 (triangle-triangle-3 triangle)))
    (let ((new-triangle-1 (build-triangle vertex-1 vertex-2 point))
	  (new-triangle-2 (build-triangle vertex-2 vertex-3 point))
	  (new-triangle-3 (build-triangle vertex-3 vertex-1 point)))
      (triangle-adjacent-triangles-set! 
	new-triangle-1
	new-triangle-2 new-triangle-3 triangle-3)
      (triangle-adjacent-triangles-set!
	new-triangle-2
	new-triangle-3 new-triangle-1 triangle-1)
      (triangle-adjacent-triangles-set!
	new-triangle-3
	new-triangle-1 new-triangle-2 triangle-2)
      (change-triangle-pointers triangle-1 triangle new-triangle-2)
      (change-triangle-pointers triangle-2 triangle new-triangle-3)
      (change-triangle-pointers triangle-3 triangle new-triangle-1)
      (remove-old-triangle-and-add-new-ones 
	triangulation triangle
	new-triangle-1 new-triangle-2 new-triangle-3)
      (let ((edge-1
	      (make-edge
		(triangle-segment-1 new-triangle-1)
		new-triangle-1 
		(triangle-triangle-1 new-triangle-1)))
	    (edge-2
	      (make-edge
		(triangle-segment-2 new-triangle-1)
		new-triangle-1
		(triangle-triangle-2 new-triangle-1)))
	    (edge-3
	      (make-edge
		(triangle-segment-3 new-triangle-1)
		new-triangle-1
		(triangle-triangle-3 new-triangle-1)))
	    (edge-4
	      (make-edge
		(triangle-segment-1 new-triangle-2)
		new-triangle-2 
		(triangle-triangle-1 new-triangle-2)))
	    (edge-5
	      (make-edge
		(triangle-segment-2 new-triangle-2)
		new-triangle-2
		(triangle-triangle-2 new-triangle-2)))
	    (edge-6
	      (make-edge
		(triangle-segment-3 new-triangle-2)
		new-triangle-2
		(triangle-triangle-3 new-triangle-2)))
	    (edge-7
	      (make-edge
		(triangle-segment-1 new-triangle-3)
		new-triangle-3 
		(triangle-triangle-1 new-triangle-3)))
	    (edge-8
	      (make-edge
		(triangle-segment-2 new-triangle-3)
		new-triangle-3
		(triangle-triangle-2 new-triangle-3)))
	    (edge-9
	      (make-edge
		(triangle-segment-3 new-triangle-3)
		new-triangle-3
		(triangle-triangle-3 new-triangle-3))))
	(legalize-edge edge-1 triangulation)
	(legalize-edge edge-2 triangulation)
	(legalize-edge edge-3 triangulation)
	(legalize-edge edge-4 triangulation)
	(legalize-edge edge-5 triangulation)
	(legalize-edge edge-6 triangulation)
	(legalize-edge edge-7 triangulation)
	(legalize-edge edge-8 triangulation)
	(legalize-edge edge-9 triangulation)))))

#| SPLIT-EDGE-BY-POINT: Similar to the previous function, but in this case four
triangles have to be handled.  The picture explains the operations:

   old triangles/segment      new triangles/segments
  
           v2                          v2
            /\                          /\
     t1<-s1/  \s2->t2                  /||<---------ns2
          / ot1\              nt1---->/ || \<----nt3
         /      \                    /  ||  \
       v1===os===v3                v1===..===v3
         \      /                    \^ || ^/
          \ ot2/              nt2---->| || |<----nt4
     t4<-s4\  /s3->t3                 |\||<---------ns4
            \/                        | \/ |
            v4                        | v4 |
                                     ns1  ns2                             |#

(define (split-edge-by-point point edge triangulation)
  (display "Splitting edge\n")
  (let ((old-triangle-1 (edge-triangle-1 edge))
	(old-triangle-2 (edge-triangle-2 edge))
	(segment (edge-segment edge)))
    (let ((vertex-1 (segment-from segment))
	  (vertex-3 (segment-to segment))
	  (vertex-2 (opposite-vertex-to-segment segment old-triangle-1))
	  (vertex-4 (opposite-vertex-to-segment segment old-triangle-2)))
      (let ((triangle-1 (opposite-triangle-to-vertex vertex-3 old-triangle-1))
	    (triangle-2 (opposite-triangle-to-vertex vertex-1 old-triangle-1))
	    (triangle-3 (opposite-triangle-to-vertex vertex-1 old-triangle-2))
	    (triangle-4 (opposite-triangle-to-vertex vertex-3 old-triangle-2))
	    (new-triangle-1 (build-triangle vertex-1 vertex-2 point))
	    (new-triangle-2 (build-triangle vertex-1 vertex-4 point))
	    (new-triangle-3 (build-triangle vertex-2 vertex-3 point))
	    (new-triangle-4 (build-triangle vertex-3 vertex-4 point)))
	(triangle-adjacent-triangles-set!
	  new-triangle-1
	  new-triangle-2 new-triangle-3 triangle-1)
	(triangle-adjacent-triangles-set!
	  new-triangle-2
	  new-triangle-1 new-triangle-4 triangle-4)
	(triangle-adjacent-triangles-set!
	  new-triangle-3
	  new-triangle-1 new-triangle-4 triangle-2)
	(triangle-adjacent-triangles-set!
	  new-triangle-4
	  new-triangle-2 new-triangle-3 triangle-3)
	(change-triangle-pointers 
	  triangle-1 old-triangle-1 new-triangle-1)
	(change-triangle-pointers 
	  triangle-2 old-triangle-1 new-triangle-3)
	(change-triangle-pointers 
	  triangle-3 old-triangle-2 new-triangle-4)
	(change-triangle-pointers 
	  triangle-4 old-triangle-2 new-triangle-2)
	(remove-old-triangle-and-add-new-ones
	  triangulation old-triangle-1
	  new-triangle-3 new-triangle-1)
	(remove-old-triangle-and-add-new-ones
	  triangulation old-triangle-2
	  new-triangle-4 new-triangle-2)
	(let ((edge-1 
		(make-edge 
		  (triangle-segment-1 new-triangle-1) 
		  new-triangle-1 
		  (triangle-triangle-1 new-triangle-1)))
	      (edge-2 
		(make-edge 
		  (triangle-segment-2 new-triangle-1) 
		  new-triangle-1 
		  (triangle-triangle-2 new-triangle-1)))
	      (edge-3 
		(make-edge 
		  (triangle-segment-3 new-triangle-1) 
		  new-triangle-1 
		  (triangle-triangle-3 new-triangle-1)))
	      (edge-4 
		(make-edge 
		  (triangle-segment-1 new-triangle-2) 
		  new-triangle-2 
		  (triangle-triangle-1 new-triangle-2)))
	      (edge-5 
		(make-edge 
		  (triangle-segment-2 new-triangle-2) 
		  new-triangle-2 
		  (triangle-triangle-2 new-triangle-2)))
	      (edge-6 
		(make-edge 
		  (triangle-segment-3 new-triangle-2) 
		  new-triangle-2 
		  (triangle-triangle-3 new-triangle-2)))
	      (edge-7 
		(make-edge 
		  (triangle-segment-1 new-triangle-3) 
		  new-triangle-3 
		  (triangle-triangle-1 new-triangle-3)))
	      (edge-8 
		(make-edge 
		  (triangle-segment-2 new-triangle-3) 
		  new-triangle-3 
		  (triangle-triangle-2 new-triangle-3)))
	      (edge-9 
		(make-edge 
		  (triangle-segment-3 new-triangle-3) 
		  new-triangle-3 
		  (triangle-triangle-3 new-triangle-3)))
	      (edge-10
		(make-edge 
		  (triangle-segment-1 new-triangle-4) 
		  new-triangle-4 
		  (triangle-triangle-1 new-triangle-4)))
	      (edge-11
		(make-edge 
		  (triangle-segment-2 new-triangle-4) 
		  new-triangle-4 
		  (triangle-triangle-2 new-triangle-4)))
	      (edge-12
		(make-edge 
		  (triangle-segment-3 new-triangle-4) 
		  new-triangle-4 
		  (triangle-triangle-3 new-triangle-4))))
	  (legalize-edge edge-1 triangulation)
	  (legalize-edge edge-2 triangulation)
	  (legalize-edge edge-3 triangulation) 
	  (legalize-edge edge-4 triangulation)
	  (legalize-edge edge-5 triangulation)
	  (legalize-edge edge-6 triangulation)
	  (legalize-edge edge-7 triangulation)
	  (legalize-edge edge-8 triangulation) 
	  (legalize-edge edge-9 triangulation)
	  (legalize-edge edge-10 triangulation)
	  (legalize-edge edge-11 triangulation) 
	  (legalize-edge edge-12 triangulation))))))

#| LEGALIZE-EDGE: This is time-consuming; it is a recursive procedure.  For
each edge we check legality --in Delaunay sense: no vertex lies in the
circumcircle of any triangle.  If an edge is not legal we flip it, and check
legality of the new edges in the old adjacent triangle.

       ilegal edge              after flip
     (v1-v3 not legal)
   
           v2                       v2
           /\                       /\
    t1<-s1/  \s2->t2               /||\
         / tu \                   / || \
      v1+==os==+v3             v1+  ns  +v3
         \ td /                   \ || /
    t4<-s4\  /s3->t3        tl---->\||/<----tr
           \/                       \/
           v4                       v4                                   |#

(define (legalize-edge edge triangulation)

  (define (flip-segment! segment triangle-up triangle-down vertex-4)
    (display "Flipping segment\n")
    (let ((vertex-1 (segment-from segment))
	  (vertex-2 (opposite-vertex-to-segment segment triangle-up))
	  (vertex-3 (segment-to segment)))
      (let ((triangle-1 
	      (opposite-triangle-to-vertex vertex-3 triangle-up))
	    (triangle-2 
	      (opposite-triangle-to-vertex vertex-1 triangle-up))
	    (triangle-3 
	      (opposite-triangle-to-vertex vertex-1 triangle-down))
	    (triangle-4 
	      (opposite-triangle-to-vertex vertex-3 triangle-down))
	    (triangle-left  (build-triangle vertex-1 vertex-2 vertex-4))
	    (triangle-right (build-triangle vertex-2 vertex-3 vertex-4)))
	(display "********************************************\n")
	(pp vertex-1)(pp vertex-2)(pp vertex-3)(pp vertex-4)
	(display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	(display-sub-triangle triangle-1)
	(display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	(display-sub-triangle triangle-2)
	(display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	(display-sub-triangle triangle-3)
	(display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	(display-sub-triangle triangle-4)
	(display "********************************************\n")
	(change-triangle-pointers triangle-1 triangle-up triangle-left)
	(change-triangle-pointers triangle-2 triangle-up triangle-right)
	(change-triangle-pointers triangle-3 triangle-down triangle-right)
	(change-triangle-pointers triangle-4 triangle-down triangle-left)
	(triangle-adjacent-triangles-set!
	  triangle-left
	  triangle-1 triangle-4 triangle-right)
	(triangle-adjacent-triangles-set!
	  triangle-right
	  triangle-2 triangle-3 triangle-left)
	(copy-triangle triangle-right triangle-up)
	(copy-triangle triangle-left triangle-down)
	(let ((edge-1 
		(make-edge 
		  (triangle-segment-1 triangle-left)
		  triangle-left
		  (triangle-triangle-1 triangle-left)))
	      (edge-2
		(make-edge
		  (triangle-segment-2 triangle-left)
		  triangle-left
		  (triangle-triangle-2 triangle-left)))
	      (edge-3
		(make-edge
		  (triangle-segment-3 triangle-left)
		  triangle-left
		  (triangle-triangle-3 triangle-left)))
	      (edge-4 
		(make-edge 
		  (triangle-segment-1 triangle-right)
		  triangle-right
		  (triangle-triangle-1 triangle-right)))
	      (edge-5 
		(make-edge 
		  (triangle-segment-2 triangle-right)
		  triangle-right
		  (triangle-triangle-2 triangle-right)))
	      (edge-6 
		(make-edge 
		  (triangle-segment-3 triangle-right)
		  triangle-right
		  (triangle-triangle-3 triangle-right))))
	  (legalize-edge edge-1 triangulation)
	  (legalize-edge edge-2 triangulation)
	  (legalize-edge edge-3 triangulation)
	  (legalize-edge edge-4 triangulation)
	  (legalize-edge edge-5 triangulation)
	  (legalize-edge edge-6 triangulation)))))

  (let ((triangle-up (edge-triangle-1 edge))
	(triangle-down (edge-triangle-2 edge))
	(segment (edge-segment edge)))
    ;; notice that, by construction, triangle-up is always the target triangle
    ;; we found in the search performed by point-support; therefore, it is
    ;; always a proper triangle.  But triangle-down might be a #f.
    (if (triangle? triangle-down)
      (let ((vertex-4 (opposite-vertex-to-segment segment triangle-down)))
	(if (point-inside-circumcircle? vertex-4 triangle-up)
	  (flip-segment! segment triangle-up triangle-down vertex-4))))))

#| The previous functions perform many calls to the following procedures. All
of these below are low-level and based in the current structure of points,
segments, edges and triangles. |#

(define (point-support point triangulation)
  (let ((x (point-x point))
	(y (point-y point)))
    (let loop ((index (- (collection-length triangulation) 1)))
      (if (< index 0)
	(error "input point not in domain")
	(let ((triangle (collection-ref triangulation index)))
	  (if (triangle? triangle) ;; it could be a #f
	    (if (point-in-triangle-closure? point triangle)
	      (let ((segment-1 (triangle-segment-1 triangle))
		    (segment-2 (triangle-segment-2 triangle))
		    (segment-3 (triangle-segment-3 triangle)))
		(cond ((zero? ((segment-functional segment-1) x y))
		       (make-edge 
			 segment-1 triangle (triangle-triangle-1 triangle)))
		      ((zero? ((segment-functional segment-2) x y))
		       (make-edge 
			 segment-2 triangle (triangle-triangle-2 triangle)))
		      ((zero? ((segment-functional segment-3) x y))
		       (make-edge 
			 segment-3 triangle (triangle-triangle-3 triangle)))
		      (else triangle)))
	      (loop (- index 1)))
	    (loop (- index 1))))))))

(define (change-triangle-pointers triangle old-triangle new-triangle)
  (if (triangle? triangle)
    (cond ((same-triangles? old-triangle (triangle-triangle-1 triangle))
	   (triangle-triangle-1-set! triangle new-triangle))
	  ((same-triangles? old-triangle (triangle-triangle-2 triangle))
	   (triangle-triangle-2-set! triangle new-triangle))
	  ((same-triangles? old-triangle (triangle-triangle-3 triangle))
	   (triangle-triangle-3-set! triangle new-triangle)))))

(define (opposite-vertex-to-segment segment triangle)
  (cond ((same-segments? segment (triangle-segment-1 triangle))
	 (triangle-vertex-1 triangle))
	((same-segments? segment (triangle-segment-2 triangle))
	 (triangle-vertex-2 triangle))
	((same-segments? segment (triangle-segment-3 triangle))
	 (triangle-vertex-3 triangle))
	(else (error "input segment not present in triangle"))))

(define (opposite-triangle-to-vertex vertex triangle)
  (cond ((equal? vertex (triangle-vertex-1 triangle))
	 (triangle-triangle-1 triangle))
	((equal? vertex (triangle-vertex-2 triangle))
	 (triangle-triangle-2 triangle))
	((equal? vertex (triangle-vertex-3 triangle))
	 (triangle-triangle-3 triangle))
	(else (error "input vertex not present in triangle"))))

(define (triangle-adjacent-triangles-set! triangle . adjacent-triangles)

  (define (triangle-adjacent-triangle-set! triangle adjacent-triangle)
    (cond ((not (vertex-in-triangle? 
		  (triangle-vertex-1 triangle) adjacent-triangle))
	   (triangle-triangle-1-set! triangle adjacent-triangle))
	  ((not (vertex-in-triangle? 
		  (triangle-vertex-2 triangle) adjacent-triangle))
	   (triangle-triangle-2-set! triangle adjacent-triangle))
	  (else 
	    (triangle-triangle-3-set! triangle adjacent-triangle))))

  (do ((index (- (length adjacent-triangles) 1) (- index 1)))
    ((< index 0) void)
    (let ((adjacent-triangle (list-ref adjacent-triangles index)))
      (if (triangle? adjacent-triangle)
	(triangle-adjacent-triangle-set!
	  triangle adjacent-triangle)))))

(define (vertex-in-triangle? vertex triangle)
  (if (triangle? triangle)
    (or (equal? vertex (triangle-vertex-1 triangle))
	(equal? vertex (triangle-vertex-2 triangle))
	(equal? vertex (triangle-vertex-3 triangle)))
    (error "no triangle")))

(define (remove-old-triangle-and-add-new-ones 
	  triangulation old-triangle . new-triangles)
  (do ((index (- (length new-triangles) 1) (- index 1)))
    ((= index 0) (copy-triangle (list-ref new-triangles 0) old-triangle))
    (add-elements-to-collection triangulation (list-ref new-triangles index))))

(define (copy-triangle source-triangle target-triangle)
  (triangle-vertex-1-set!   target-triangle 
			    (triangle-vertex-1 source-triangle))
  (triangle-vertex-2-set!   target-triangle 
			    (triangle-vertex-2 source-triangle))
  (triangle-vertex-3-set!   target-triangle 
			    (triangle-vertex-3 source-triangle))
  (triangle-segment-1-set!  target-triangle 
			    (triangle-segment-1 source-triangle))
  (triangle-segment-2-set!  target-triangle 
			    (triangle-segment-2 source-triangle))
  (triangle-segment-3-set!  target-triangle 
			    (triangle-segment-3 source-triangle))
  (triangle-triangle-1-set! target-triangle 
			    (triangle-triangle-1 source-triangle))
  (triangle-triangle-2-set! target-triangle 
			    (triangle-triangle-2 source-triangle))
  (triangle-triangle-3-set! target-triangle 
			    (triangle-triangle-3 source-triangle)))
