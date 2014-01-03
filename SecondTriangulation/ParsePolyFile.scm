(define (PolyFile+angle->delaunay-triangulation 
	  file-name threshold-angle)
  (let* ((port  (open-input-file file-name))
	 (size  (read port))
	 (dim   (read port))
	 (dummy (read port))
	 (dummy (read port)))
    (let loop ((position 0) (result '()))
      (if (= position size)
	(let ((polygon (list->vector result))
	      (size  (read port))
	      (dummy (read port)))
	  (let loop ((position 0) (result '()))
	    (if (= position size)
	      (with-output-to-file
		"triangulation.gnu"
		(lambda () (write-triangulation 
			     (delaunay-triangulation threshold-angle 
						     polygon 
						     (list->vector result)))))
	      (let ((dummy   (read port))
		    (point-1 (- (read port) 1))
		    (point-2 (- (read port) 1)))
		(loop (+ position 1)
		      (append 
			result
			(list (make-segment
				(vector-ref polygon point-1)
				(vector-ref polygon point-2)))))))))
	(let ((name (read port))
	      (x    (read port))
	      (y    (read port)))
	  (loop (+ position 1)
		(append result
			(list (make-point x y name)))))))))
    
(define (PolyFile->gnu file-name)
  (let* ((port  (open-input-file file-name))
	 (size  (read port))
	 (dim   (read port))
	 (dummy (read port))
	 (dummy (read port)))
    (let loop ((position 0) (result '()))
      (if (= position size)
	(let ((polygon (list->vector result)) 
	      (size (read port))
	      (dummy (read port)))
	  (let loop ((position 0))
	    (if (= position size)
	      (newline)
	      (let ((dummy (read port))
		    (point-1 (- (read port) 1))
		    (point-2 (- (read port) 1)))
		(begin
		  (display (point-x (vector-ref polygon point-1)))
		  (display " ")
		  (display (point-y (vector-ref polygon point-1)))
		  (newline)
		  (display (point-x (vector-ref polygon point-2)))
		  (display " ")
		  (display (point-y (vector-ref polygon point-2)))
		  (newline)
		  (newline)
		  (loop (+ position 1)))))))
	(let ((name (read port))
	      (x    (read port))
	      (y    (read port)))
	  (loop (+ position 1)
		(append
		  result
		  (list (make-point x y name)))))))))
	  
(define (PolyFile->gnuscript poly-file script-file)
  (with-output-to-file "polygon.gnu"
		       (lambda () (PolyFile->gnu poly-file)))
  (with-output-to-file 
    script-file
    (lambda ()
      (display "set terminal postscript\n")
      (display "plot \"polygon.gnu\" w l\n")
      (display "quit"))))

(define (file->polygon-triangulation file)
  (let* ((port (open-input-file file))
	 (size (read port))
	 (dim  (read port))
	 (dummy (read port))
	 (dummy (read port)))
    (let loop ((position 0) (result '()))
      (if (= position size)
	(let ((polygon (list->vector result))
	      (size (read port))
	      (dummy (read port)))
	  (let loop ((position 0) (result '()))
	    (if (= position size)
	      (make-super-polygon 
		polygon 
		(list->vector result) 
		(polygon-triangulation polygon))
	      (let ((dummy (read port))
		    (point-1 (- (read port) 1))
		    (point-2 (- (read port) 1)))
		(loop (+ position 1)
		      (append
			result
			(list (make-segment (vector-ref polygon point-1)
					    (vector-ref polygon point-2)))))))))
	(let ((name (read port))
	      (x (read port))
	      (y (read port)))
	  (loop (+ position 1)
		(append
		  result
		  (list (make-point x y name)))))))))
	      
	 
