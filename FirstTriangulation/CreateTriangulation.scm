(declare (flonum))

(define-structure point x y name)
(define-structure triangle vertex-1 vertex-2 vertex-3 ptriangle-1-2 ptriangle-2-3 ptriangle-3-1)
;; a triangulation is just a vector of triangles,
;; so we can access each triangle very fast.

(define no-triangle -1)
(define origen (make-point 0.0 0.0 'O))
(define no-segment (list origen origen))
(define pi 3.14159265358969)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (determinant-4 a b c d e f g h i j k l m n o p)
  (- (+ (* a f k p) (* a j h o) (* a n g l) (* e b l o)
	(* e j c p) (* e n d k) (* i b g p) (* i f d o)
	(* i n c h) (* m b h k) (* m f c l) (* m j d g)) 
     (+ (* a f l o) (* a j g p) (* a n h k) (* e b k p)
	(* e j d o) (* e n c l) (* i b h o) (* i f c p)
	(* i n d g) (* m b g l) (* m f d k) (* m j c h))))

(define (discriminant a0 b0 a1 b1 a2 b2);; the discriminant is just
  (- (* (- b0 b1)                       ;; plugging the point (a0,b0)
	(- a2 a1))                      ;; in the line defined by the
     (* (- b2 b1)                       ;; points (a1,b1) and (a2,b2). 
	(- a0 a1))))                    ;; If it's there, you get zero
                                        ;; if it's not, you get + or -

(define (sign a)
  (cond ((> a 0)  1)
	((< a 0) -1)
	(else     0)))

(define (not-legal? edge-vertex-1 edge-vertex-2 point opposite-point)
  ;; it checks if opposite-point is inside the circumcircle
  ;; defined by edge-vertices and point.
  (let* ((x0 (point-x opposite-point))
	 (y0 (point-y opposite-point))
	 (x1 (point-x edge-vertex-1))
	 (y1 (point-y edge-vertex-1))
	 (x2 (point-x edge-vertex-2))
	 (y2 (point-y edge-vertex-2))
	 (x3 (point-x point))
	 (y3 (point-y point))
	 (circumcircle (determinant-4 1 x0 y0 (+ (* x0 x0) (* y0 y0))
				      1 x1 y1 (+ (* x1 x1) (* y1 y1))
				      1 x2 y2 (+ (* x2 x2) (* y2 y2))
				      1 x3 y3 (+ (* x3 x3) (* y3 y3))))
	 (orientation (discriminant x1 y1 x2 y2 x3 y3)))
    (> (* circumcircle orientation) 0)))

(define (polygon-boundary polygon)
  (let ((poly-vector (list->vector polygon)))
    (let loop ((position 0) (xmin 0) (xmax 0) (ymin 0) (ymax 0))
      (if (= (vector-length poly-vector) position)
	  (let* ((xlength (/ (- xmax xmin) 2.0))
		 (ylength (/ (- ymax ymin) 2.0))
		 (N (make-point (+ xmax xlength xlength) (+ ymax ylength ylength) 'N))
		 (E (make-point (+ xmax xlength xlength) (- ymin ylength ylength) 'E))
		 (S (make-point (- xmin xlength xlength) (- ymin ylength ylength) 'S))
		 (W (make-point (- xmin xlength xlength) (+ ymax ylength ylength) 'W)))
	    (list->vector (list N E S W)))
	  (let ((vertex (vector-ref poly-vector position)))
	    (loop (+ position 1) 
		  (min xmin (point-x vertex)) (max xmax (point-x vertex))
		  (min ymin (point-y vertex)) (max ymax (point-y vertex))))))))

(define (naive-triangulation polygon)
  (let ((initial-point (car polygon))
	(N (vector-ref (polygon-boundary polygon) 0))
	(E (vector-ref (polygon-boundary polygon) 1))
	(S (vector-ref (polygon-boundary polygon) 2))
	(W (vector-ref (polygon-boundary polygon) 3)))
    (list->vector (list
		   (make-triangle N E initial-point -1 1 3) 
		   (make-triangle E S initial-point -1 2 0)
		   (make-triangle S W initial-point -1 3 1)
		   (make-triangle W N initial-point -1 0 2)))))

(define (the-other-vertex vertex-1 vertex-2 triangle)
  (cond ((and (not (equal? (triangle-vertex-1 triangle) vertex-1)) 
	      (not (equal? (triangle-vertex-1 triangle) vertex-2))) 
	 (triangle-vertex-1 triangle))
	((and (not (equal? (triangle-vertex-2 triangle) vertex-1)) 
	      (not (equal? (triangle-vertex-2 triangle) vertex-2))) 
	 (triangle-vertex-2 triangle))
	((and (not (equal? (triangle-vertex-3 triangle) vertex-1)) 
	      (not (equal? (triangle-vertex-3 triangle) vertex-2))) 
	 (triangle-vertex-3 triangle))))

(define (pointer-associated-to-edge edge-vertex-1 edge-vertex-2 triangle)
  (let ((triangle-vertex-1 (triangle-vertex-1 triangle))
	(triangle-vertex-2 (triangle-vertex-2 triangle))
	(triangle-vertex-3 (triangle-vertex-3 triangle)))
    (cond ((or (and (equal? edge-vertex-1 triangle-vertex-1) (equal? edge-vertex-2 triangle-vertex-2)) 
	       (and (equal? edge-vertex-1 triangle-vertex-2) (equal? edge-vertex-2 triangle-vertex-1)))
	   (triangle-ptriangle-1-2 triangle))
	  ((or (and (equal? edge-vertex-1 triangle-vertex-2) (equal? edge-vertex-2 triangle-vertex-3)) 
	       (and (equal? edge-vertex-1 triangle-vertex-3) (equal? edge-vertex-2 triangle-vertex-2)))
	   (triangle-ptriangle-2-3 triangle))
	  ((or (and (equal? edge-vertex-1 triangle-vertex-1) (equal? edge-vertex-2 triangle-vertex-3)) 
	       (and (equal? edge-vertex-1 triangle-vertex-3) (equal? edge-vertex-2 triangle-vertex-1)))
	   (triangle-ptriangle-3-1 triangle)))))

(define (change-pointers-in-ptriangle ptriangle old-pointer new-pointer triangulation)
  (if (not (= ptriangle no-triangle))
      (let ((triangle (vector-ref triangulation ptriangle)))
	(cond ((= (triangle-ptriangle-1-2 triangle) old-pointer)
	       (triangle-ptriangle-1-2-set! triangle new-pointer))
	      ((= (triangle-ptriangle-2-3 triangle) old-pointer)
	       (triangle-ptriangle-2-3-set! triangle new-pointer))
	      ((= (triangle-ptriangle-3-1 triangle) old-pointer)
	       (triangle-ptriangle-3-1-set! triangle new-pointer))))))

(define (update-edge-in-triangulation vertex-edge-1 vertex-edge-2 target-ptriangle triangulation)
  (let* ((target-triangle (vector-ref triangulation target-ptriangle))
	 (edge-ptriangle  (pointer-associated-to-edge vertex-edge-1 vertex-edge-2 target-triangle)))
    (if (not (= edge-ptriangle no-triangle))
	(let* ((edge-triangle   (vector-ref triangulation edge-ptriangle))
	       (target-vertex   (the-other-vertex vertex-edge-1 vertex-edge-2 target-triangle))
	       (opposite-vertex (the-other-vertex vertex-edge-1 vertex-edge-2 edge-triangle))
	       (alpha-ptriangle (pointer-associated-to-edge vertex-edge-1 target-vertex target-triangle))
	       (beta-ptriangle  (pointer-associated-to-edge vertex-edge-2 target-vertex target-triangle))
	       (gamma-ptriangle (pointer-associated-to-edge vertex-edge-1 opposite-vertex edge-triangle))
	       (delta-ptriangle (pointer-associated-to-edge vertex-edge-2 opposite-vertex edge-triangle)))
	  (begin
	    (triangle-vertex-1-set! target-triangle vertex-edge-1)
	    (triangle-vertex-2-set! target-triangle target-vertex)
	    (triangle-vertex-3-set! target-triangle opposite-vertex)
	    (triangle-ptriangle-1-2-set! target-triangle alpha-ptriangle)
	    (triangle-ptriangle-2-3-set! target-triangle edge-ptriangle)
	    (triangle-ptriangle-3-1-set! target-triangle gamma-ptriangle)
	    (if (not (= beta-ptriangle no-triangle))
		(change-pointers-in-ptriangle beta-ptriangle target-ptriangle edge-ptriangle triangulation))
	    (if (not (= gamma-ptriangle no-triangle))
		(change-pointers-in-ptriangle gamma-ptriangle edge-ptriangle target-ptriangle triangulation))
	    (triangle-vertex-1-set! edge-triangle target-vertex)
	    (triangle-vertex-2-set! edge-triangle vertex-edge-2)
	    (triangle-vertex-3-set! edge-triangle opposite-vertex)
	    (triangle-ptriangle-1-2-set! edge-triangle beta-ptriangle)
	    (triangle-ptriangle-2-3-set! edge-triangle delta-ptriangle)
	    (triangle-ptriangle-3-1-set! edge-triangle target-ptriangle)
	    (if (not (= gamma-ptriangle no-triangle))
		(if (not-legal? vertex-edge-1 
				opposite-vertex 
				target-vertex 
				(the-other-vertex vertex-edge-1 opposite-vertex 
						  (vector-ref triangulation gamma-ptriangle)))
		    (begin
		      ;(display "Updating edge [") (display (point-name vertex-edge-1))
		      ;(display "-") (display (point-name opposite-vertex)) (display "].\n")
		      (update-edge-in-triangulation vertex-edge-1 opposite-vertex 
						    target-ptriangle triangulation))))
	    (if (not (= delta-ptriangle no-triangle))
		(if (not-legal? vertex-edge-2
				opposite-vertex
				target-vertex
				(the-other-vertex vertex-edge-2 opposite-vertex
						  (vector-ref triangulation delta-ptriangle)))
		    (begin
		      ;(display "Updating edge [") (display (point-name vertex-edge-2))
		      ;(display "-") (display (point-name opposite-vertex)) (display "].\n")
		      (update-edge-in-triangulation vertex-edge-2 opposite-vertex 
						    edge-ptriangle triangulation)))))))))

(define (polygon-triangulation polygon)
  (let loop ((next-vertices (cdr polygon)) (triangulation (naive-triangulation polygon)))
    (if (= (length next-vertices) 1)
	(split-something (car next-vertices) triangulation)
	(loop (cdr next-vertices) (split-something (car next-vertices) triangulation)))))

(define (look-for-edge-in-triangulation point triangulation)
  (define (point-inside-edge? point vertex-1 vertex-2)
    (let* ((x1 (point-x vertex-1))
	   (y1 (point-y vertex-1))
	   (x2 (point-x vertex-2))
	   (y2 (point-y vertex-2))
	   (x0 (point-x point))
	   (y0 (point-y point))
	   (center-x (/ (+ x1 x2) 2))
	   (center-y (/ (+ y1 y2) 2))
	   (radius (+ (* (- x1 center-x) (- x1 center-x)) 
		      (* (- y1 center-y) (- y1 center-y))))
	   (distance (+ (* (- x0 center-x) (- x0 center-x)) 
			(* (- y0 center-y) (- y0 center-y)))))
      (and (= (discriminant x0 y0 x1 y1 x2 y2) 0)
	   (< distance radius))))
  (let loop ((position 0))
    (if (= position (vector-length triangulation))
	(list->vector (list no-triangle origen origen))
	(let* ((triangle (vector-ref triangulation position))
	       (vertex-1 (triangle-vertex-1 triangle))
	       (vertex-2 (triangle-vertex-2 triangle))
	       (vertex-3 (triangle-vertex-3 triangle)))
	  (cond ((point-inside-edge? point vertex-1 vertex-2) (list->vector (list position vertex-1 vertex-2)))
		((point-inside-edge? point vertex-2 vertex-3) (list->vector (list position vertex-2 vertex-3)))
		((point-inside-edge? point vertex-3 vertex-1) (list->vector (list position vertex-3 vertex-1)))
		(else (loop (+ position 1))))))))

(define (look-for-ptriangle-in-triangulation point triangulation)
  (define (point-inside-triangle? point triangle) 
    (let ((x1 (point-x (triangle-vertex-1 triangle))) 
	  (y1 (point-y (triangle-vertex-1 triangle)))
	  (x2 (point-x (triangle-vertex-2 triangle))) 
	  (y2 (point-y (triangle-vertex-2 triangle)))
	  (x3 (point-x (triangle-vertex-3 triangle))) 
	  (y3 (point-y (triangle-vertex-3 triangle)))
	  (x0 (point-x point))
	  (y0 (point-y point)))
      (let ((discriminant-1 (discriminant x0 y0 x1 y1 x2 y2))
	    (discriminant-2 (discriminant x0 y0 x2 y2 x3 y3))
	    (discriminant-3 (discriminant x0 y0 x3 y3 x1 y1)))
	(and (= (sign discriminant-1) (sign discriminant-2))
	     (= (sign discriminant-1) (sign discriminant-3))))))
  (let loop ((position 0))
    (if (= position (vector-length triangulation))
	(error "The point seems to be in neither edge nor triangle")
	(if (point-inside-triangle? point (vector-ref triangulation position))
	    position
	    (loop (+ position 1))))))

(define (split-something point triangulation)
  (if (= (vector-ref (look-for-edge-in-triangulation point triangulation) 0) no-triangle)
      (split-triangle-in-triangulation point 
				       (look-for-ptriangle-in-triangulation point triangulation)
				       triangulation)
      (split-edge-in-triangulation point (look-for-edge-in-triangulation point triangulation) triangulation)))

(define (split-triangle-in-triangulation point ptriangle triangulation)
  (let* ((length (vector-length triangulation))
	 (triangle (vector-ref triangulation ptriangle))
	 (vertex-1 (triangle-vertex-1 triangle))
	 (vertex-2 (triangle-vertex-2 triangle))
	 (vertex-3 (triangle-vertex-3 triangle))
	 (ptriangle-1-2 (triangle-ptriangle-1-2 triangle))
	 (ptriangle-2-3 (triangle-ptriangle-2-3 triangle))
	 (ptriangle-3-1 (triangle-ptriangle-3-1 triangle))
	 (new-triangle-1 (make-triangle vertex-2 vertex-3 point
					ptriangle-2-3 (+ length 1) ptriangle))
	 (new-triangle-2 (make-triangle vertex-3 vertex-1 point
					ptriangle-3-1 ptriangle length))
	 (new-triangulation (list->vector (append (vector->list triangulation) 
						  (list new-triangle-1 new-triangle-2)))))
    (begin
      (triangle-vertex-1-set! triangle vertex-1)
      (triangle-vertex-2-set! triangle vertex-2)
      (triangle-vertex-3-set! triangle point)
      (triangle-ptriangle-2-3-set! triangle length)
      (triangle-ptriangle-3-1-set! triangle (+ length 1))
      (change-pointers-in-ptriangle ptriangle-2-3 ptriangle length new-triangulation)
      (change-pointers-in-ptriangle ptriangle-3-1 ptriangle (+ length 1) new-triangulation)
      (if (not (= ptriangle-1-2 no-triangle))
	  (if (not-legal? vertex-1 vertex-2 
			  point (the-other-vertex vertex-1 vertex-2 (vector-ref new-triangulation ptriangle-1-2)))
	      (update-edge-in-triangulation vertex-1 vertex-2 ptriangle new-triangulation)))
      (if (not (= ptriangle-2-3 no-triangle))
	  (if (not-legal? vertex-2 vertex-3 
			  point (the-other-vertex vertex-2 vertex-3 (vector-ref new-triangulation ptriangle-2-3)))
	      (update-edge-in-triangulation vertex-2 vertex-3 length new-triangulation)))
      (if (not (= ptriangle-3-1 no-triangle))
	  (if (not-legal? vertex-3 vertex-1
			  point (the-other-vertex vertex-3 vertex-1 (vector-ref new-triangulation ptriangle-3-1)))
	      (update-edge-in-triangulation vertex-3 vertex-1 (+ length 1) new-triangulation)))
      new-triangulation)))

(define (split-edge-in-triangulation point triple triangulation)
  (let* ((length (vector-length triangulation))
	 (ptriangle      (vector-ref triple 0))
	 (edge-vertex-1  (vector-ref triple 1))
	 (edge-vertex-2  (vector-ref triple 2))
	 (triangle       (vector-ref triangulation ptriangle))
	 (target-vertex  (the-other-vertex edge-vertex-1 edge-vertex-2 triangle))
	 (edge-ptriangle (pointer-associated-to-edge edge-vertex-1 edge-vertex-2 triangle))
	 (alpha-ptriangle (pointer-associated-to-edge edge-vertex-1 target-vertex triangle))
	 (beta-ptriangle  (pointer-associated-to-edge edge-vertex-2 target-vertex triangle)))
    (cond ((= edge-ptriangle no-triangle)
	   (let* ((new-triangle (make-triangle target-vertex edge-vertex-2 point
					       beta-ptriangle no-triangle ptriangle))
		  (new-triangulation (list->vector (append (vector->list triangulation) (list new-triangle)))))
	     (begin
	       (triangle-vertex-1-set! triangle edge-vertex-1)
	       (triangle-vertex-2-set! triangle target-vertex)
	       (triangle-vertex-3-set! triangle point)
	       (triangle-ptriangle-1-2-set! triangle alpha-ptriangle)
	       (triangle-ptriangle-2-3-set! triangle length)
	       (triangle-ptriangle-3-1-set! triangle no-triangle)
	       (change-pointers-in-ptriangle beta-ptriangle ptriangle length new-triangulation)
	       (if (not (= alpha-ptriangle no-triangle))
		   (if (not-legal? target-vertex edge-vertex-1
				   point (the-other-vertex target-vertex edge-vertex-1 
							   (vector-ref new-triangulation alpha-ptriangle)))
		       (update-edge-in-triangulation target-vertex edge-vertex-1 ptriangle new-triangulation)))
	       (if (not (= beta-ptriangle no-triangle))
		   (if (not-legal? target-vertex edge-vertex-2
				   point (the-other-vertex target-vertex edge-vertex-2
							   (vector-ref new-triangulation beta-ptriangle)))
		       (update-edge-in-triangulation target-vertex edge-vertex-2 
						     length new-triangulation)))
	       new-triangulation)))
	  (else (let* ((length (vector-length triangulation))
		       (edge-triangle   (vector-ref triangulation edge-ptriangle))
		       (opposite-vertex (the-other-vertex edge-vertex-1 edge-vertex-2 edge-triangle))
		       (gamma-ptriangle (pointer-associated-to-edge edge-vertex-1 opposite-vertex edge-triangle))
		       (delta-ptriangle (pointer-associated-to-edge edge-vertex-2 opposite-vertex edge-triangle))
		       (up-triangle     (make-triangle target-vertex edge-vertex-2 point
						       beta-ptriangle (+ length 1) ptriangle))
		       (down-triangle   (make-triangle opposite-vertex point edge-vertex-2
						       edge-ptriangle length delta-ptriangle))
		       (new-triangulation (list->vector (append (vector->list triangulation) 
								(list up-triangle down-triangle)))))
		  (begin
		    (triangle-vertex-1-set! triangle edge-vertex-1)
		    (triangle-vertex-2-set! triangle target-vertex)
		    (triangle-vertex-3-set! triangle point)
		    (triangle-ptriangle-1-2-set! triangle alpha-ptriangle)
		    (triangle-ptriangle-2-3-set! triangle length)
		    (triangle-ptriangle-3-1-set! triangle edge-ptriangle)
		    (triangle-vertex-1-set! edge-triangle edge-vertex-1)
		    (triangle-vertex-2-set! edge-triangle point)
		    (triangle-vertex-3-set! edge-triangle opposite-vertex)
		    (triangle-ptriangle-1-2-set! edge-triangle ptriangle)
		    (triangle-ptriangle-2-3-set! edge-triangle (+ length 1))
		    (triangle-ptriangle-3-1-set! edge-triangle gamma-ptriangle)
		    (change-pointers-in-ptriangle beta-ptriangle ptriangle length new-triangulation)
		    (change-pointers-in-ptriangle delta-ptriangle edge-ptriangle (+ length 1) new-triangulation)
		    (if (not-legal? edge-vertex-1 point target-vertex opposite-vertex)
			(update-edge-in-triangulation edge-vertex-1 point ptriangle new-triangulation))
		    (if (not-legal? edge-vertex-2 point target-vertex opposite-vertex)
			(update-edge-in-triangulation edge-vertex-2 point length new-triangulation))
		    (if (not (= alpha-ptriangle no-triangle))
			(if (not-legal? target-vertex edge-vertex-1
					point (the-other-vertex target-vertex edge-vertex-1 
								(vector-ref new-triangulation alpha-ptriangle)))
			    (update-edge-in-triangulation target-vertex edge-vertex-1
							  ptriangle new-triangulation)))
		    (if (not (= beta-ptriangle no-triangle))
			(if (not-legal? target-vertex edge-vertex-2
					point (the-other-vertex target-vertex edge-vertex-2
								(vector-ref new-triangulation beta-ptriangle)))
			    (update-edge-in-triangulation target-vertex edge-vertex-2 
							  length new-triangulation)))
		    (if (not (= gamma-ptriangle no-triangle))
			(if (not-legal? opposite-vertex edge-vertex-1
					point (the-other-vertex opposite-vertex edge-vertex-1 
								(vector-ref new-triangulation gamma-ptriangle)))
			    (update-edge-in-triangulation opposite-vertex edge-vertex-1 
							  edge-ptriangle new-triangulation)))
		    (if (not (= delta-ptriangle no-triangle))
			(if (not-legal? opposite-vertex edge-vertex-2
					point (the-other-vertex opposite-vertex edge-vertex-2
								(vector-ref new-triangulation delta-ptriangle)))
			    (update-edge-in-triangulation opposite-vertex edge-vertex-2
							  (+ length 1) new-triangulation)))
		    new-triangulation))))))
