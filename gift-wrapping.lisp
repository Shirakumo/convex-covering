(in-package #:org.shirakumo.fraf.convex-covering)

(defvar *data*)
(defvar *3d*)
(defvar *2d*)
(defvar *leftmost*)
(defvar *hull* )

(defvar *hull-v*)
(defvar *hull-f*)

(defun orthonormal-base2 (base-a base-b base-c normal)
  (let* ((origin base-a)
         (base-x (vunit (v- base-b origin)))
         (base-y (let ((y* (v- base-c origin)))
                   (vunit (v- y* (v* base-x (v. base-x y*)))))))
    (values base-x base-y)
    #+no (if (< (v. (vc base-x base-y) normal) 0)
        (values base-x base-y)
        (values base-x (v- base-y)))))

(defun gift-wrap (vertices base-a base-b base-c normal)
  (multiple-value-bind (base-x base-y)
      (orthonormal-base2 base-a base-b base-c normal)
    #+no (format *trace-output* "O ~A~%X ~A ~A~%Y ~A ~A~%. ~A~%" o x (vlength x) y (vlength y) (v. x y))
    (loop with count = (floor (length vertices) 3)
          with centroid = (dvec3)
          for i below count
          for v3 = (manifolds:v vertices i)
          for v2 = (vec (v. base-x v3) (v. base-y v3))
          do (nv+ centroid v3)
          collect v2 into v2d
          :finally (let ((boundary (jarvis (coerce v2d 'vector)))
                         (flip (< (v. (vc base-x base-y) normal) 0)))
                     (return (case (length boundary)
                               (3
                                (one-triangle vertices boundary flip))
                               (4 ; TODO(jmoringe): two-triangles
                                (triangle-fan vertices centroid boundary flip))
                               (t
                                (triangle-fan vertices centroid boundary flip))))))))

(defun one-triangle (vertices boundary flip)
  (let ((vertices* (make-array (* 3 3) :element-type 'double-float))
        (faces (make-array 3 :element-type 'manifolds:u32
                             :initial-contents (if flip
                                                   '(0 1 2)
                                                   '(0 2 1)))))
    (loop for i from 0 below 3
          for j = (aref boundary i)
          for vertex-fill-pointer from 0 by 3
          for v = (manifolds:v vertices j)
          do (setf (aref vertices* (+ vertex-fill-pointer 0)) (vx v)
                   (aref vertices* (+ vertex-fill-pointer 1)) (vy v)
                   (aref vertices* (+ vertex-fill-pointer 2)) (vz v)))
    (values vertices* faces)))

(defun triangle-fan (vertices centroid boundary flip)
  (let* ((count (floor (length vertices) 3))
         (centroid (v/ centroid count))
         (vertices* (make-array (* 3 (1+ (length boundary))) :element-type 'double-float))
         (faces (make-array (* 3 (length boundary)) :element-type 'manifolds:u32)))
    (setf (aref vertices* 0) (vx centroid)
          (aref vertices* 1) (vy centroid)
          (aref vertices* 2) (vz centroid))
    (loop for j from 0 below (length boundary)
          for k = (mod (1+ j) (length boundary))
          for l = (aref boundary j)
          for vertex-fill-pointer from 3 by 3
          for face-fill-pointer   from 0 by 3
          do (let ((v (manifolds:v vertices l)))
               (setf (aref vertices* (+ vertex-fill-pointer 0)) (vx v)
                     (aref vertices* (+ vertex-fill-pointer 1)) (vy v)
                     (aref vertices* (+ vertex-fill-pointer 2)) (vz v)))
             (setf (values (aref faces (+ face-fill-pointer 0))
                           (aref faces (+ face-fill-pointer 1))
                           (aref faces (+ face-fill-pointer 2)) )
                   (if flip
                       (values 0 (1+ j) (1+ k))
                       (values 0 (1+ k) (1+ j))))
          finally (progn
                    #+no (format *trace-output* "~A~%~A%" vertices* faces)
                    (setf *hull-v* vertices* *hull-f* faces)))
    (values vertices* faces)))

;; (manifolds:centroid )
;;; https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
(defun jarvis (vertices) ; 2d points
  (let ((leftmost (loop with index = 0
                        with min = (vx (aref vertices index))
                        for i from 1 below (length vertices)
                        for x = (vx (aref vertices i))
                        when (< x min)
                        do (setf min x index i)
                        finally (progn
                                  (setf *leftmost* index)
                                  (return
                                    index))))
        (hull-points (make-array 0 :adjustable T :fill-pointer 0)))
    (loop with point-on-hull = leftmost
          for count :from 0
          for vpoint-on-hull = (aref vertices point-on-hull)
          for end-point = (loop with best = 0
                                for j below (length vertices)
                                when (or (= best point-on-hull)
                                         (let* ((vbest (aref vertices best))
                                                (vcandidate (aref vertices j))
                                                (line (v- vbest vpoint-on-hull))
                                                (line* (v- vcandidate vpoint-on-hull)))
                                        ; (format *trace-output* "  ~A ~A => ~A~%" point-on-hull j (vz (vc (vec line 0.0) (vec line* 0.0))))
                                           (let ((cross (vc (vec line 0.0) (vec line* 0.0))))
                                             (and (plusp (vz cross))
                                                  (or (and (= j leftmost)
                                                           #+no (> (vlength (vc (vunit (vec line 0.0))
                                                                           (vunit (vec line* 0.0))))
                                                              .00001))
                                                      (and (/= j point-on-hull)
                                                           (not (find j hull-points))))))))
                                do (setf best j)
                                finally (progn
                                        ; (format *trace-output* "best ~A~%" best)
                                          (return best)))
          do (vector-push-extend point-on-hull hull-points)
             (setf point-on-hull end-point)
             (when (or (= point-on-hull leftmost ; (aref hull-points 0)
                          )
                       (= (length hull-points) (length vertices))) ; TODO: last resort
               #+no (setf *data* (list :v2d vertices
                                  :hull (loop for i across hull-points
                                              collect (aref vertices i))
                                  :point vpoint-on-hull
                                  ; :candidate (aref vertices j)
                                  ; :best (aref vertices best)
                                  ))
               #+no (break)
               (return hull-points)))))

(jarvis (coerce (setf *2d* (loop :repeat 100
                                 :collect (vrand (vec 0 0) (vec 5 5))))
                'vector))

(v. (vec 0.0 1.0 0.0) (vec 5.0 0.0 0.0))
