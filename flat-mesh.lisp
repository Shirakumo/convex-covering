(in-package #:org.shirakumo.fraf.convex-covering)

;;; From VERTICES and the plane described by BASE-{A,B,C} and NORMAL,
;;; compute a flat mesh. Return as two values the vector of vertex
;;; positions and the vector of face indices.
(defun flat-mesh (vertices base-a base-b base-c normal
                  &key (triangulation :triangle-fan))
  (check-type triangulation (member :triangle-fan :subdivide))
  (multiple-value-bind (base2x base2y) (orthonormal-base2 base-a base-b base-c)
    (loop with count = (floor (length vertices) 3)
          with points2 = (make-array count)
          for i below count
          for v3 = (manifolds:v vertices i)
          for v2 = (vec (v. base2x v3) (v. base2y v3))
          do (setf (aref points2 i) v2)
          finally (let ((boundary (jarvis points2))
                        (flip (< (v. (vc base2x base2y) normal) 0)))
                    (return
                      (ecase triangulation
                        (:subdivide
                         (subdivide vertices boundary flip))
                        (:triangle-fan
                         (case (length boundary)
                           (3
                            (one-triangle vertices boundary flip))
                           (T
                            (triangle-fan vertices boundary flip))))))))))

;;; Compute and return two 2d vectors that form an orthonormal base in
;;; the 2d plane which is embedded according to BASE-{A,B,C},
(defun orthonormal-base2 (base-a base-b base-c)
  (let* ((origin base-a)
         (base-x (vunit (v- base-b origin)))
         (base-y (let ((y* (v- base-c origin)))
                   (vunit (v- y* (v* base-x (v. base-x y*)))))))
    (values base-x base-y)))

;;; Given the set of 2d points, POINTS, return sequence of indices
;;; into POINTS that describe a convex boundary which all POINTS.
;;; See https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
(defun jarvis (points) ; 2d points
  (let ((leftmost (loop with index = 0
                        with min = (vx (aref points index))
                        for i from 1 below (length points)
                        for x = (vx (aref points i))
                        when (< x min)
                        do (setf min x index i)
                        finally (return index)))
        (hull-points (make-array 0 :adjustable T :fill-pointer 0)))
    (loop with point-on-hull = leftmost
          for count from 0
          for vpoint-on-hull = (aref points point-on-hull)
          for end-point
             = (loop with best = 0
                     for j below (length points)
                     when (or (= best point-on-hull)
                              (let* ((vbest (aref points best))
                                     (vcandidate (aref points j))
                                     (line (v- vbest vpoint-on-hull))
                                     (line* (v- vcandidate vpoint-on-hull)))
                                (let ((cross (vc (vec line 0.0) (vec line* 0.0))))
                                  (and (plusp (vz cross))
                                       (or (and (= j leftmost))
                                           (and (/= j point-on-hull)
                                                (not (find j hull-points))))))))
                       do (setf best j)
                     finally (return best))
          do (vector-push-extend point-on-hull hull-points)
             (setf point-on-hull end-point)
             (when (or (= point-on-hull leftmost)
                       (= (length hull-points) (length points))) ; last resort
               (return hull-points)))))

;;; Triangulation methods

(defun subdivide (all-vertices boundary-indices flip)
  (let* ((vertex-count (length boundary-indices))
         (face-count (- vertex-count 2))
         (vertices (make-array (* 3 vertex-count)
                               :element-type (array-element-type all-vertices)))
         (vertex-fill-pointer 0)
         (faces (make-array (* 3 face-count)))
         (face-fill-pointer 0))
    (labels ((add-vertex (i)
               (setf (manifolds:v vertices vertex-fill-pointer)
                     (manifolds:v all-vertices i))
               (1- (incf vertex-fill-pointer)))
             (add-triangle (i j k)
               (setf (aref faces (+ face-fill-pointer 0)) i
                     (aref faces (+ face-fill-pointer 1)) j
                     (aref faces (+ face-fill-pointer 2)) k)
               (incf face-fill-pointer 3))
             (rec (remaining)
               (let ((count (length remaining)))
                 (case count
                   (3 (destructuring-bind (i j k) remaining
                        (if flip
                            (add-triangle i j k)
                            (add-triangle i k j))))
                   (t (let ((split (floor count 2)))
                        (rec (subseq remaining 0 (1+ split)))
                        (rec (append (subseq remaining split)
                                     (list (elt remaining 0))))))))))
      (rec (map 'list #'add-vertex boundary-indices))
      (values vertices faces))))

(defun one-triangle (all-vertices boundary flip)
  (let ((vertices (make-array (* 3 3) :element-type 'double-float))
        (faces (make-array 3 :element-type 'manifolds:u32
                             :initial-contents (if flip
                                                   '(0 1 2)
                                                   '(0 2 1)))))
    (loop for i from 0 below 3
          for j = (aref boundary i)
          do (setf (manifolds:v vertices i)
                   (manifolds:v all-vertices j)))
    (values vertices faces)))

(defun triangle-fan (all-vertices boundary flip)
  (let* ((centroid (dvec3))
         (vertices (make-array (* 3 (1+ (length boundary))) :element-type 'double-float))
         (faces (make-array (* 3 (length boundary)) :element-type 'manifolds:u32)))
    (loop for j from 0 below (length boundary)
          for k = (mod (1+ j) (length boundary))
          for l = (aref boundary j)
          for vertex-fill-pointer from 1
          for face-fill-pointer from 0 by 3
          for v = (manifolds:v all-vertices l)
          do (nv+ centroid v)
             (setf (manifolds:v vertices vertex-fill-pointer) v)
             (setf (values (aref faces (+ face-fill-pointer 0))
                           (aref faces (+ face-fill-pointer 1))
                           (aref faces (+ face-fill-pointer 2)))
                   (if flip
                       (values 0 (1+ j) (1+ k))
                       (values 0 (1+ k) (1+ j)))))
    (setf (manifolds:v vertices 0) (v/ centroid (length boundary)))
    (values vertices faces)))
