(define (minimum-angle structure)
  (cond ((triangle? structure)
	 (let ((vertex-1 (triangle-vertex-1 structure))
	       (vertex-2 (triangle-vertex-2 structure))
	       (vertex-3 (triangle-vertex-3 structure)))
	   (min (abs (angle vertex-1 vertex-2 vertex-3))
		(abs (angle vertex-2 vertex-3 vertex-1))
		(abs (angle vertex-3 vertex-1 vertex-2)))))
	((list? structure)
	 (let* ((poly-vector (list->vector structure))
		(last-vertex   (vector-ref poly-vector (- (vector-length poly-vector) 1)))
		(first-vertex  (vector-ref poly-vector 0))
		(second-vertex (vector-ref poly-vector 1)))
	   (let loop ((position 1) (result (angle last-vertex first-vertex second-vertex)))
	     (if (= position (- (length structure) 1))
		 (let ((previous-vertex (vector-ref poly-vector (- position 1))))
		   (min result (angle previous-vertex last-vertex first-vertex)))
		 (let ((previous-vertex (vector-ref poly-vector (- position 1)))
		       (current-vertex  (vector-ref poly-vector position))
		       (next-vertex     (vector-ref poly-vector (+ position 1))))
		   (loop (+ position 1) (min result (angle previous-vertex current-vertex next-vertex))))))))
	((super-polygon? structure)
	 (let ((triangulation (super-polygon-triangulation structure))
	       (segments      (super-polygon-segments structure)))
	   (let loop ((position 0) (result 360) (entry 0))
	     (if (= position (vector-length triangulation))
		 (list result entry)
		 (let* ((triangle (vector-ref triangulation position))
			(test (minimum-angle triangle)))
		   (if (and (< test result)
			    (is-inner-triangle? triangle segments))
		       (loop (+ position 1) test position)
		       (loop (+ position 1) result entry)))))))
	(else (error "unable to extract minimum angle of structure"))))

(define (is-boundary-triangle? triangle segments)
  (let ((segment-1 (list (triangle-vertex-1 triangle) (triangle-vertex-2 triangle)))
	(segment-2 (list (triangle-vertex-2 triangle) (triangle-vertex-3 triangle)))
	(segment-3 (list (triangle-vertex-3 triangle) (triangle-vertex-1 triangle))))
    (let loop ((position 0))
      (if (= position (vector-length segments))
	  (= 0 1)
	  (let ((segment (vector-ref segments position)))
	    (if (and (not (same-segments? segment-1 segment))
		     (not (same-segments? segment-2 segment))
		     (not (same-segments? segment-3 segment)))
		(loop (+ position 1))
		(= 0 0)))))))

(define (is-corner-triangle? triangle segments)
  (let ((segment-1 (list (triangle-vertex-1 triangle) (triangle-vertex-2 triangle)))
	(segment-2 (list (triangle-vertex-2 triangle) (triangle-vertex-3 triangle)))
	(segment-3 (list (triangle-vertex-3 triangle) (triangle-vertex-1 triangle))))
    (let loop ((position 1))
      (if (= position (- (vector-length segments) 1))
	  (let ((segment-before (vector-ref segments (- position 1)))
		(segment-after  (vector-ref segments position)))
	    (or (and (equal? segment-1 segment-before)
		     (equal? segment-2 segment-after))
		(and (equal? segment-1 segment-after)
		     (equal? segment-2 segment-before))
		(and (equal? segment-2 segment-before)
		     (equal? segment-3 segment-after))
		(and (equal? segment-2 segment-after)
		     (equal? segment-3 segment-before))
		(and (equal? segment-3 segment-before)
		     (equal? segment-1 segment-after))
		(and (equal? segment-3 segment-after)
		     (equal? segment-1 segment-before))))
	  (let ((segment-before (vector-ref segments (- position 1)))
		(segment-after  (vector-ref segments position)))
	    (if (or (and (equal? segment-1 segment-before)
			 (equal? segment-2 segment-after))
		    (and (equal? segment-1 segment-after)
			 (equal? segment-2 segment-before))
		    (and (equal? segment-2 segment-before)
			 (equal? segment-3 segment-after))
		    (and (equal? segment-2 segment-after)
			 (equal? segment-3 segment-before))
		    (and (equal? segment-3 segment-before)
			 (equal? segment-1 segment-after))
		    (and (equal? segment-3 segment-after)
			 (equal? segment-1 segment-before)))
		(= 0 0)
		(loop (+ position 1))))))))

(define (corner-angle corner-triangle segments)
  (let ((segment-1 (list (triangle-vertex-1 corner-triangle) (triangle-vertex-2 corner-triangle)))
	(segment-2 (list (triangle-vertex-2 corner-triangle) (triangle-vertex-3 corner-triangle)))
	(segment-3 (list (triangle-vertex-3 corner-triangle) (triangle-vertex-1 corner-triangle))))
    (let loop ((position 1))
      (if (= position (vector-length segments))
	  (error "triangle is not a corner triangle")
	  (let ((segment-before (vector-ref segments (- position 1)))
		(segment-after  (vector-ref segments position)))
	    (cond ((and (equal? segment-1 segment-before)
			(equal? segment-2 segment-after))
		   (angle (car segment-1) (car segment-2) (cadr segment-2)))
		  ((and (equal? segment-1 segment-after)
			(equal? segment-2 segment-before))
		   (angle (car segment-2) (car segment-1) (cadr segment-1)))
		  ((and (equal? segment-2 segment-before)
			(equal? segment-3 segment-after))
		   (angle (car segment-2) (car segment-3) (cadr segment-3)))
		  ((and (equal? segment-2 segment-after)
			(equal? segment-3 segment-before))
		   (angle (car segment-3) (car segment-2) (cadr segment-2)))
		  ((and (equal? segment-3 segment-before)
			(equal? segment-1 segment-after))
		   (angle (car segment-3) (car segment-1) (cadr segment-1)))
		  ((and (equal? segment-3 segment-after)
			(equal? segment-1 segment-before))
		   (angle (car segment-1) (car segment-3) (cadr segment-3)))
		  (else (loop (+ position 1)))))))))

(define (easy-interior-point triangle)
  (let* ((vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (x1 (point-x vertex-1))
	 (y1 (point-y vertex-1))
	 (x2 (point-x vertex-2))
	 (y2 (point-y vertex-2))
	 (x3 (point-x vertex-3))
	 (y3 (point-y vertex-3)))
    (make-point (/ (+ x1 x2 x3) 3.0) (/ (+ y1 y2 y3) 3.0) 'I)))

(define (point-inside-polygon? point segments)
  (let loop ((position 0) (result 0))
    (if (= position (vector-length segments))
	(> (abs result) 10)
	(let ((segment (vector-ref segments position)))
	  (loop (+ position 1) (+ result (angle (car segment) point (cadr segment))))))))

(define (is-inner-triangle? triangle segments) (point-inside-polygon? (easy-interior-point triangle) segments))
		
    