;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Definition of dynamic-vector (Dvector) and related basic functions
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-structure Dvector body valid-entries)

(define (build-Dvector n #!optional initial-value)
  (make-Dvector 
    (make-vector (one-half-more-plus-two n) initial-value)
    0))

(define (Dvector-ref Dvector index)
  (declare (fixnum))
  (if (< index (Dvector-valid-entries Dvector))
    (vector-ref (Dvector-body Dvector) index)
    (error "out of range")))

(define (Dvector-set! Dvector index value)
  (declare (fixnum))
  (let ((size (Dvector-valid-entries Dvector)))
    (cond ((< index size)
	   (vector-set! (Dvector-body Dvector) index value))
	  ((< index (vector-length (Dvector-body Dvector)))
	   (vector-set! (Dvector-body Dvector) index value)
	   (Dvector-valid-entries-set! Dvector (+  index 1)))
	  (else
	    (let ((new-body
		    (make-vector 
		      (one-half-more-plus-two (+ index 1))
		      #f))
		  (body (Dvector-body Dvector))
		  (size (Dvector-valid-entries Dvector)))
	      (do ((counter (- size 1) (- counter 1)))
		((< counter 0)
		 (Dvector-body-set! Dvector new-body)
		 (vector-set! new-body index value)
		 (Dvector-valid-entries-set! Dvector (+ index 1)))
		(vector-set! 
		  new-body counter (vector-ref body counter))))))))

(define Dvector-length Dvector-valid-entries)

(define (Dvector-append Dvector value)
  (Dvector-set! Dvector (Dvector-valid-entries Dvector) value))

(define (one-half-more-plus-two x)
  (declare (fixnum))
  (+ 2 (inexact->exact (round (* 1.5 x)))))
