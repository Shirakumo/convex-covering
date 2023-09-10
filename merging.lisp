(in-package #:org.shirakumo.fraf.convex-covering)

;;; Merge validity

(defun valid-patch-p (patch all-vertices all-faces)
  (when *debug-visualizations*
    (handler-case
        (progn
          #+no (let ((hull (patch-hull patch)))
            (loop :for face-index :below (/ (length (patch-faces patch)) 3)
                  :do (debug-face* all-vertices (patch-faces patch) face-index :patch-edge hull :sample-count 11)))

          (let ((hull (patch-hull patch)))
            (loop :with hull-vertices = (vertices hull)
                  :with hull-faces = (faces hull)
                  :for face-index :below (/ (length hull-faces) 3)
                  :for face-vertex-index = (aref hull-faces (* 3 face-index))
                  :for face-vertex = (manifolds:v hull-vertices face-vertex-index)
                  :for face-normal = (vunit (manifolds:face-normal hull-vertices hull-faces face-index))
                  :for face-centroid = (manifolds:centroid hull-vertices (subseq hull-faces (* 3 face-index) (+ (* 3 face-index) 3)))
                  :do (debug-face* hull-vertices hull-faces face-index :facet-edge hull :sample-count 12)
                  :do (push (cons face-centroid :facet-centroid) (problem hull))
                      (debug-line face-centroid (v* face-normal .2) :facet-normal hull :sample-count 11)
                  :collect (cons face-centroid face-normal))))
      (floating-point-invalid-operation ()
        )))

  (let ((valid
          (handler-case
              (let ((hull (patch-hull patch)))
                (and (not (find-vertex-in-hull hull all-vertices all-faces))
                     (not (find-edge-in-hull hull all-vertices all-faces))
                     (normals-matching-p patch hull all-vertices all-faces)
                     (not (find-boundary-constraint-bar hull))
                     (not (find-negative-side-touching-triangle hull))))
            (floating-point-invalid-operation (error)
              (format *error-output* "~A~%" error)
              nil))))
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
        ; :when *debug-visualizations*
          ; :do (push (cons face-centroid :facet-centroid) (problem hull))
            ;   (debug-line face-centroid (v* face-normal .3) :facet-normal hull)
        :collect (cons face-centroid face-normal)))

(defun ensure-facet-normals (hull)
  (or (face-normals hull)
      (setf (face-normals hull) (compute-facet-normals hull))))


(defun vertex-in-hull-p (vertex hull &key (eps -.001))
  (handler-case
   (loop :for (facet-centroid . facet-normal) :in (ensure-facet-normals hull) ; TODO rename face-normals -> %face-normals; ensure-face-normals -> face-normals
         :always
                                        ; (minusp (v. facet-normal (v- vertex facet-centroid)))
                                        ; (not (> (v. facet-normal (v- vertex facet-centroid)) .0001))
            (< (v. facet-normal (v- vertex facet-centroid)) eps))
    (nil #+no floating-point-invalid-operation () ; HACK
      t)))

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
                                                         hull :sample-count 13))
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

(defun face-matches-facet-p (h1 h2 h3 v1 v2 v3
                             &key (eps .000001))
  (flet ((d (p1 p2)
           (v2norm (v- p1 p2))))
    (flet ((matching-facet-vertex (p)
             (loop :for h :in (list h1 h2 h3)
                   :when (< (d p h) eps)
                   :return h)))
      (let ((m1 (matching-facet-vertex v1))
            (m2 (matching-facet-vertex v2))
            (m3 (matching-facet-vertex v3)))
        (and m1 m2 m3)))))

(defun face-shares-vertex-or-edge-with-facet-p (h1 h2 h3 v1 v2 v3
                                                &key (eps .00001) hull)
  (flet ((d (p1 p2)
           (v2norm (v- p1 p2))))
    (flet ((matching-facet-vertex (p)
             (loop :for h :in (list h1 h2 h3)
                   :when (< (d p h) eps)
                   :return h))
           (bad-point-p (p)
             (cond ((not (point-on-triangle-p h1 h2 h3 p))
                    (when (and hull *debug-visualizations*)
                      (push (cons (v+ p (dvec (random .1) (random .1) (random .1))) :bad) (annotations hull)))
                    t)
                   (t
                    (when (and hull *debug-visualizations*)
                      (push (cons (v+ p (dvec (random .1) (random .1) (random .1))) :good) (annotations hull)))
                    nil))))
      (let ((m1 (matching-facet-vertex v1))
            (m2 (matching-facet-vertex v2))
            (m3 (matching-facet-vertex v3)))
        #+no (when (and hull *debug-visualizations*)
          (when m1
            (push (cons m1 :good) (annotations hull)))
          (when m2
            (push (cons m2 :good) (annotations hull)))
          (when m3
            (push (cons m3 :good) (annotations hull))))
        (and (<= 1 (count-if #'identity (list m1 m2 m3)) 2)
             (or (and (null m1) (bad-point-p v1))
                 (and (null m2) (bad-point-p v2))
                 (and (null m3) (bad-point-p v3))))))))

(defun face-overlaps-or-matches-facet-old-p (hull-vertices hull-faces facet-index v1 v2 v3
                                         &key (eps .000001))
  (flet ((d (p1 p2)
           (v2norm (v- p1 p2))))
    (let* ((h1 (manifolds:v hull-vertices (aref hull-faces (+ (* 3 facet-index) 0))))
           (h2 (manifolds:v hull-vertices (aref hull-faces (+ (* 3 facet-index) 1))))
           (h3 (manifolds:v hull-vertices (aref hull-faces (+ (* 3 facet-index) 2))))
           (c1 (manifolds:closest-point-on-triangle
                hull-vertices hull-faces facet-index v1))
           (i1 (< (d c1 v1) eps)) ; vertex inside triangle
           (c2 (manifolds:closest-point-on-triangle
                hull-vertices hull-faces facet-index v2))
           (i2 (< (d c2 v2) eps))
           (c3 (manifolds:closest-point-on-triangle
                hull-vertices hull-faces facet-index v3))
           (i3 (< (d c3 v3) eps)))
      (and (or i1 i2 i3)      ; at least one must be on or in the face
           (flet ((matching-facet-vertex (p)
                    (loop :for h :in (list h1 h2 h3)
                          :when (< (d p h) eps)
                          :return h)))
             (let ((m1 (matching-facet-vertex c1))
                   (m2 (matching-facet-vertex c2))
                   (m3 (matching-facet-vertex c3)))
               (or (and m1 m2 m3 (not (eq m1 m2)) (not (eq m2 m3)) (not (eq m3 m1))) ; match
                   (or (and i1 (null m1))
                       (and i2 (null m2))
                       (and i3 (null m3))) ; proper overlap (i.e. not just sharing an edge)
                   )))))))

(defun face-overlaps-or-matches-facet-p (hull-vertices hull-faces facet-index v1 v2 v3
                                         &key (threshold 0) hull)
  (let ((h1 (manifolds:v hull-vertices (aref hull-faces (+ (* 3 facet-index) 0))))
        (h2 (manifolds:v hull-vertices (aref hull-faces (+ (* 3 facet-index) 1))))
        (h3 (manifolds:v hull-vertices (aref hull-faces (+ (* 3 facet-index) 2)))))
    #+alternative (and (let ((facet-normal (vunit (vc (v- h2 h1) (v- h3 h1))))
                             (face-normal (vunit (vc (v- v2 v1) (v- v3 v1)))))
                         (< (abs (- 1 (abs (v. facet-normal face-normal)))) .00001))
                       (not (face-shares-vertex-or-edge-with-facet-p h1 h2 h3 v1 v2 v3 :hull hull))
                       (triangles-intersect-p h1 h2 h3 v1 v2 v3 :threshold threshold))
    #+old (or (face-matches-facet-p h1 h2 h3 v1 v2 v3)
              (and (not (face-shares-vertex-or-edge-with-facet-p h1 h2 h3 v1 v2 v3 :hull hull))
                   (triangles-intersect-p h1 h2 h3 v1 v2 v3 :threshold threshold)))
    (or (face-matches-facet-p h1 h2 h3 v1 v2 v3)
        (multiple-value-bind (intersectp constellation contact)
            (triangles-intersect-p h1 h2 h3 v1 v2 v3 :threshold .0001)
          #+no (when (and (eq constellation :coplanar) (eq contact :penetrating))
            (let ((*debug-output* t))
              (triangles-intersect-p h1 h2 h3 v1 v2 v3)))
          (and (eq constellation :coplanar) (eq contact :penetrating))))))

(defun triangles-coplanar-and-intersect-p (u1 u2 u3 v1 v2 v3 &key threshold)
  (let* ((normal (vc (v- u2 u1) (v- u3 u1)))
         (w1     (v. normal (v- u1 v1))))
    ))

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
            (flet (#+unused (point-in-facet-p (point facet-index)
                              (let ((closest (manifolds:closest-point-on-triangle
                                              hull-vertices hull-faces facet-index point)))
                                (< (v2norm (v- closest point)) .001)))) ; TODO magic number
              (when (not (loop :for i :below (/ (length hull-faces) 3)
                               :for (facet-centroid . facet-normal) :in (face-normals hull)
                               :always (if #+no (or (point-in-facet-p a i)
                                                    (point-in-facet-p b i)
                                                    (point-in-facet-p c i))
                                           (face-overlaps-or-matches-facet-p
                                            hull-vertices hull-faces i a b c :hull hull)
                                           (let ( ; (*debug-output* t)
                                                 )
                                             (d "[~A ~A ~A] possibly intersects facet ~A~%"
                                                ia ib ic i)
                                             (let ((face-normal (vunit (vc (v- b a) (v- c a)))))
                                               (when *debug-visualizations*
                                                 (debug-line* a b :face-edge hull)
                                                 (debug-line* b c :face-edge hull)
                                                 (debug-line* c a :face-edge hull)
                                                 (let ((face-centroid (v/ (v+ a b c) 3)))
                                                   (push (cons face-centroid :face-centroid) (annotations hull))
                                                   (debug-line face-centroid (v* face-normal .3) :face-normal hull :sample-count 13)))
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
