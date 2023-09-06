(in-package #:org.shirakumo.fraf.convex-covering)

;;; Merge validity

(defun valid-patch-p (patch all-vertices all-faces)
  #+no (let ((hull (patch-hull patch)))
         (when (null (face-normals hull))
           (loop :with hull-vertices = (vertices hull)
                 :with hull-faces = (faces hull)
                 :for face-index :below (/ (length hull-faces) 3)
                 :for face-vertex-index = (aref hull-faces (* 3 face-index))
                 :for face-vertex = (manifolds:v hull-vertices face-vertex-index)
                 :for face-normal = (handler-case
                                        (vunit (manifolds:face-normal hull-vertices hull-faces face-index))
                                      (error ()
                                        (dvec 1 0 0)))
                 :do (push face-normal (face-normals hull))
                 :do (loop :with sample-count = 32
                           :with d = (v/ face-normal (float sample-count 1.0d0))
                           :repeat sample-count
                           :for v = face-vertex :then (v+ v d)
                           :do (push (cons v :normal) (problem hull))))))

  (let ((valid
          (let ((hull (patch-hull patch)))
            (and (not (find-vertex-in-hull hull all-vertices all-faces))
                 (not (find-edge-in-hull hull all-vertices all-faces))
                 (normals-matching-p patch hull all-vertices all-faces)
                 (not (find-boundary-constraint-bar hull))
                 (not (find-negative-side-touching-triangle hull))))))
    (d "; hull ~A is ~:[invalid~;valid~]"
       (global-faces (patch-hull patch)) valid)
    valid))

;;; Individual criteria

(defun compute-facet-normals (hull) ; TODO these are just face normals for now
  (loop :with hull-vertices = (vertices hull)
        :with hull-faces = (faces hull)
        :for face-index :below (/ (length hull-faces) 3)
        :for face-vertex-index = (aref hull-faces (* 3 face-index))
        :for face-vertex = (manifolds:v hull-vertices face-vertex-index)
        :for face-normal = (vunit (manifolds:face-normal hull-vertices hull-faces face-index))
        :for face-centroid = (manifolds:centroid hull-vertices (subseq hull-faces (* 3 face-index) (+ (* 3 face-index) 3)))
        :when *debug-visualizations*
          :do (push (cons face-centroid :facet-centroid) (problem hull))
              (debug-line face-centroid (v* face-normal .1) :facet-normal hull)
        :collect (cons face-centroid face-normal)))

(defun ensure-facet-normals (hull)
  (or (face-normals hull)
      (setf (face-normals hull) (compute-facet-normals hull))))


(defun vertex-in-hull-p (vertex hull)
  (handler-case
   (loop :for (facet-centroid . facet-normal) :in (ensure-facet-normals hull) ; TODO rename face-normals -> %face-normals; ensure-face-normals -> face-normals
         :always
                                        ; (minusp (v. facet-normal (v- vertex facet-centroid)))
                                        ; (not (> (v. facet-normal (v- vertex facet-centroid)) .0001))
            (< (v. facet-normal (v- vertex facet-centroid)) -.0001))
    (floating-point-invalid-operation () ; HACK
      t))


  #+no (loop :with hull-vertices = (vertices hull)
        :with hull-faces = (faces hull)
        :for face-index :below (/ (length hull-faces) 3)
        :for face-vertex-index = (aref hull-faces (* 3 face-index))
        :for face-vertex = (manifolds:v hull-vertices face-vertex-index)
        :for face-normal = (handler-case
                               (vunit (manifolds:face-normal hull-vertices hull-faces face-index))
                             (error ()
                               (dvec 1 0 0)))
           #+no :do #+no (format *trace-output* "~A ~A~&  => v ~A n ~A~&  => ~A ~A~%"
                                 vertex face-vertex
                                 (v- face-vertex vertex) face-normal
                                 (v. face-normal (v- vertex face-vertex))
                                 (plusp (v. face-normal (v- vertex face-vertex))))
        :always
           ; (minusp (v. face-normal (v- vertex face-vertex)))
           (not (plusp (v. face-normal (v- vertex face-vertex))))
           ; (< (v. face-normal (v- vertex face-vertex)) .000001)
        ))


(defun find-vertex-in-hull (hull all-vertices all-faces)
  (let ((result (loop :with hull-faces = (global-faces hull)
                      :for vertex-index :from 0 :below (/ (length all-vertices) 3)
                      :for vertex = (manifolds:v all-vertices vertex-index)
                      :when (and *debug-visualizations*
                                 (not (find vertex-index hull-faces))
                                 (vertex-in-hull-p vertex hull))
                        :do (push (cons vertex :bad) (annotations hull))
                      :thereis (and (not (find vertex-index hull-faces))
                                    (vertex-in-hull-p vertex hull)))))
    (when result
      (setf (problem hull) :vertex))
    result))

(defun find-edge-in-hull (hull all-vertices all-faces)
  (when (let (                          ; (hull-faces (faces hull))
              (hull-faces/global (global-faces hull)))
          (manifolds:do-faces (a b c all-faces nil) ; TODO tests some edges multiple times
            #+no (format *trace-output* "; [~A ~A ~A] vs hull ~A~%"
                         a b c
                         hull-faces/global
                                        ; (not (not (search (vector a b c) hull-faces/global)))
                         )
            (when (and ; (not (search (vector a b c) hull-faces/global))
                   (let (          ; (va (manifolds:v all-vertices a))
                                        ; (vb (manifolds:v all-vertices b))
                                        ; (vc (manifolds:v all-vertices c))
                         )
                     (flet ((edge-in-hull-p (vi1 vi2)
                              #+no (format *trace-output* "; [~A ~A] part of hull ~A => ~A~%"
                                           vi1 vi2
                                           hull-faces/global
                                           (loop :for i :below (length hull-faces/global) :by 3
                                                 :thereis (and (find vi1 hull-faces/global :start i :end (+ i 3))
                                                               (find vi2 hull-faces/global :start i :end (+ i 3)))))
                              (and (not (let ((result (loop :for i :below (length hull-faces/global) :by 3
                                                            :thereis (and (find vi1 hull-faces/global :start i :end (+ i 3))
                                                                          (find vi2 hull-faces/global :start i :end (+ i 3))))))
                                          (when (and result *debug-visualizations*)
                                            (debug-line* (manifolds:v all-vertices vi1)
                                                         (manifolds:v all-vertices vi2)
                                                         :boundary-edge
                                                         hull))
                                          result))
                                   (loop :with v1 = (manifolds:v all-vertices vi1)
                                         :with v2 = (manifolds:v all-vertices vi2)
                                         :with bad = nil
                                         :with sample-count = 8
                                         :repeat sample-count
                                         :with d = (v/ (v- v2 v1) sample-count)
                                         :for v = (v+ v1 (v* d (/ .5d0 sample-count))) :then (v+ v d)
                                         :for in? = (vertex-in-hull-p v hull)
                                         :when in?
                                           :do (d  "; [~A ~A]~&;~2@T~A - ~A~&;~2@Tat ~A~%"
                                                   vi1 vi2 v1 v2 v)
                                               (setf bad t)
                                         ; :do (push (cons v (if in? :bad :good)) (problem hull))
                                        ; :thereis (vertex-in-hull-p v hull all-vertices)
                                         :finally (return bad)))))
                       (or (edge-in-hull-p a b)
                           (edge-in-hull-p b c)
                           (edge-in-hull-p c a)))))
              (return t))))
    (setf (problem hull) :edge)
    t))

(defun face-key (vertices a b c)
  (labels ((v-lex-< (v1 v2)
             (cond ((= (vx v1) (vx v2))
                    (cond ((= (vy v1) (vy v2))
                           (cond ((= (vz v1) (vz v2))
                                  nil)
                                 ((< (vz v1) (vz v2))
                                  t)
                                 (t nil)))
                          ((< (vy v1) (vy v2))
                           t)
                          (t nil)))
                   ((< (vx v1) (vx v2))
                    t)
                   (t nil))))
    (let ((a (manifolds:v vertices a))
          (b (manifolds:v vertices b))
          (c (manifolds:v vertices c)))
      (sort (vector a b c) #'v-lex-<))))

(defun normals-matching-p-old (patch hull all-vertices all-faces)
  (loop :with patch-faces = (patch-faces patch)
        :for face-index :below (/ (length patch-faces) 3)
        :for face-normal = (vunit (manifolds:face-normal all-vertices patch-faces face-index))
        :for face-centroid = (manifolds:centroid all-vertices (subseq patch-faces (* 3 face-index) (+ (* 3 face-index) 3)))
        :do (debug-line face-centroid (v* face-normal .1) :face-normal hull))


  (or (<= (length (vertices hull)) (* 3 4)) ; HACK to work around degenerate 2D case
      (let* ((indices->normal (make-hash-table :test #'equalp)))
        (manifolds:do-faces (a b c all-faces)
          (let ((key (face-key all-vertices a b c)))
            (assert (not (gethash key indices->normal)))
            (setf (gethash key indices->normal)
                  (let ((a (manifolds:v all-vertices a))
                        (b (manifolds:v all-vertices b))
                        (c (manifolds:v all-vertices c)))
                    (vunit (vc (v- b a) (v- c a)))))))
        (let ((result (loop :with hull-vertices = (vertices hull)
                                        ; :with hull-faces = (faces hull)
                            :with hull-faces = (faces hull)
                            :for face-index :below (/ (length hull-faces) 3)
                            :for face-vertex-indices = (subseq hull-faces (* 3 face-index) (+ (* 3 face-index) 3))
                            :for face-centroid = (manifolds:centroid hull-vertices (subseq hull-faces (* 3 face-index) (+ (* 3 face-index) 3)))
                            :for key = (face-key hull-vertices
                                                 (aref face-vertex-indices 0)
                                                 (aref face-vertex-indices 1)
                                                 (aref face-vertex-indices 2))
                            :for (nil . facet-normal) :in (face-normals hull)
                            :for face-normal = (gethash key indices->normal)
                            :for diff = (when face-normal
                                          (abs (- 1 (v. face-normal facet-normal))))
                            :for good? = (or (null diff)
                                             (< diff .00001))
                            :do (d "facet ~A~&~2@T-> ~A~&~2@Tvs ~A~&~2@T-> ~A~%"
                                   face-vertex-indices facet-normal face-normal diff)
                            :when (and (not good?) *debug-visualizations*)
                              :do (loop :with sample-count = 32
                                        :with d = (v/ face-normal (float sample-count 1.0d0))
                                        :repeat sample-count
                                        :for v = face-centroid :then (v+ v d)
                                        :do (push (cons v :bad-normal) (annotations hull)))
                                  (loop :with sample-count = 32
                                        :with d = (v/ facet-normal (float sample-count 1.0d0))
                                        :repeat sample-count
                                        :for v = face-centroid :then (v+ v d)
                                        :do (push (cons v :bad-normal) (annotations hull)))
                            :always good?)))
          (unless result
            (setf (problem hull) :normals))
          result))))

(defun normals-matching-p (patch hull all-vertices all-faces)
  #+no (loop :with patch-faces = (patch-faces patch)
        :for face-index :below (/ (length patch-faces) 3)
        :for face-normal = (vunit (manifolds:face-normal all-vertices patch-faces face-index))
        :for face-centroid = (manifolds:centroid all-vertices (subseq patch-faces (* 3 face-index) (+ (* 3 face-index) 3)))
        :do (loop :with sample-count = 32
                  :with d = (v/ face-normal (float sample-count 1.0d0))
                  :repeat sample-count
                  :for v = face-centroid :then (v+ v d)
                  :do (push (cons v :face-normal) (problem hull))))
  (or (<= (length (vertices hull)) (* 3 4)) ; HACK to work around degenerate 2D case
      (manifolds:do-faces (ia ib ic all-faces t)
        (let ((a (manifolds:v all-vertices ia))
              (b (manifolds:v all-vertices ib))
              (c (manifolds:v all-vertices ic)))
          (let ((hull-vertices (vertices hull))
                (hull-faces (faces hull)))
            (flet ((point-in-facet-p (point facet-index)
                     (let ((closest (manifolds:closest-point-on-triangle
                                     hull-vertices hull-faces facet-index point)))
                       (< (v2norm (v- closest point)) .0001))))
              (when (not (loop :for i :below (/ (length hull-faces) 3)
                               :for (facet-centroid . facet-normal) :in (face-normals hull)
                               :always (if (or (point-in-facet-p a i)
                                               (point-in-facet-p b i)
                                               (point-in-facet-p c i))
                                           (progn
                                             (d "[~A ~A ~A] possibly intersects facet ~A~%"
                                                ia ib ic i)
                                             (let ((face-normal (vunit (vc (v- b a) (v- c a)))))
                                               (d "~2@Tface normal ~A~&~2@Tfacet normal ~A~&~2@T=> ~A~%"
                                                  face-normal
                                                  facet-normal
                                                  (v. face-normal facet-normal))
                                               (cond ((< (abs (- -1 (v. face-normal facet-normal))) .0001)
                                                      (when *debug-visualizations*
                                                        (debug-line facet-centroid face-normal :bad-normal hull)
                                                        (debug-line facet-centroid facet-normal :bad-normal hull))
                                                      nil)
                                                     (t
                                                      t))))
                                           t)))
                (setf (problem hull) :normals)
                (return nil))))))))

(defun find-boundary-constraint-bar (hull)
  nil)

(defun find-negative-side-touching-triangle (hull)
  nil)

;;; Merge score

(defun compute-compactness (all-vertices patch)
  (assert (< (abs (- (manifolds:surface-area all-vertices (patch-faces patch))
                     (patch-surface-area patch)))
             .0001))
  (let* ((faces (patch-faces patch))
         (boundary-length (manifolds:boundary-length all-vertices faces)))
    (if (zerop boundary-length)
        most-positive-double-float ; TODO
        (/ (+                      (sqrt (patch-surface-area patch))
            #+no (expt (manifolds:convex-volume all-vertices faces) 1/3)
            #+no (let ((v (make-array (* (length faces) 3))))
              (loop :for i :below (length faces)
                    :for vi = (aref faces i)
                    :do (setf (aref v (+ (* 3 i) 0)) (aref all-vertices (+ vi 0))
                              (aref v (+ (* 3 i) 1)) (aref all-vertices (+ vi 1))
                              (aref v (+ (* 3 i) 2)) (aref all-vertices (+ vi 2))))
              (let ((bb (nth-value 1 (manifolds:bounding-box v))))
                (print bb)
                (expt (* (vx bb) (vy bb) (vz bb)) 1/3))))
           boundary-length))))
