(in-package #:org.shirakumo.fraf.convex-covering)

;;; Merge validity

(defun valid-patch-p (patch context)
  (when *debug-visualizations*
    (handler-case
        (progn
          #+no (let ((hull (patch-hull patch)))
            (loop :for face-index :below (/ (length (patch-faces patch)) 3)
                  :do (debug-face* all-vertices (patch-faces patch) face-index :patch-edge hull :sample-count 11)))

          (let ((hull (patch-hull patch)))
            (loop with hull-vertices = (vertices hull)
                  with hull-faces = (facets hull)
                  for face-index below (/ (length hull-faces) 3)
                  for face-vertex-index = (aref hull-faces (* 3 face-index))
                  for face-vertex = (manifolds:v hull-vertices face-vertex-index)
                  for face-normal = (vunit (manifolds:face-normal hull-vertices hull-faces face-index))
                  for face-centroid = (manifolds:centroid hull-vertices (subseq hull-faces (* 3 face-index) (+ (* 3 face-index) 3)))
                  do (debug-face* hull-vertices hull-faces face-index :facet-edge hull :sample-count 12)
                  do (push (cons face-centroid :facet-centroid) (annotations hull))
                     (debug-line face-centroid (v* face-normal .2) :facet-normal hull :sample-count 11))))
      (floating-point-invalid-operation ()
        )))

  (let ((valid
          (handler-case
              (let ((hull (patch-hull patch)))
                (and ; (not (find-vertex-in-hull hull context))
                     (not (find-edge-in-hull hull context))
                     #+TODO (normals-matching-p/overlap-per-facet patch hull context)
                     (normals-matching-p patch hull context)
                     (not (find-boundary-constraint-bar hull context))
                     (not (find-negative-side-touching-triangle hull))))
            (floating-point-invalid-operation (error)
              (format *error-output* "~a~%" error)
              NIL))))
    (d "; hull ~a is ~:[invalid~;valid~]"
       (global-faces (patch-hull patch)) valid)
    valid))

;;; Individual criteria

(defun vertex-in-hull-p (vertex hull &key (eps -.01))
  (declare (type convex-hull hull)
           (type dvec3 vertex))
  (destructuring-bind (center . size/2) (bounding-box hull)
    (declare (type vec3 center size/2))
    (and (v<= (v- center size/2) (vec3 (vx vertex) (vy vertex) (vz vertex)) (v+ center size/2)) ; TODO avoid conversion
         (loop for (facet-centroid . facet-normal) in (facet-normals hull)
               always
                                        ; (minusp (v. facet-normal (v- vertex facet-centroid)))
                                        ; (not (> (v. facet-normal (v- vertex facet-centroid)) .0001))
                  (< (v. (the dvec3 facet-normal) (v- vertex (the dvec3 facet-centroid))) eps)))))

(defun find-vertex-in-hull (hull context)
  (declare (optimize speed (safety 1) (debug 1)))
  (let* ((all-vertices (context-vertices context))
         (result (loop with hull-faces = (global-faces hull)
                       for vertex-index from 0 below (/ (length all-vertices) 3)
                       for vertex = (manifolds:v all-vertices vertex-index)
                       when (and *debug-visualizations*
                                 (not (find vertex-index hull-faces))
                                 (vertex-in-hull-p vertex hull))
                         do (push (cons vertex :bad) (annotations hull))
                      thereis (and (vertex-in-hull-p vertex hull)
                                   (not (find vertex-index hull-faces))))))
    (when result
      (setf (problem hull) :vertex))
    result))

(defvar *overlap-counts* (make-array 0 :adjustable T :fill-pointer 0))

(defun find-edge-in-hull (hull context)
  (when (let ((all-vertices (context-vertices context))
              (all-faces (context-faces context))
                                        ; (hull-faces (facets hull))
                                        ; (hull-faces/global (global-faces hull))
              (overlapping (make-array 0 :adjustable T :fill-pointer 0)))
          (unwind-protect
               (space:do-overlapping (face-info (context-face-index context) hull NIL)
                 ;; manifolds:do-faces (a b c all-faces NIL) ; TODO tests some edges multiple times
                 (let* #+no ()
                   ((face-index (face-info-index face-info))
                    (a (aref all-faces (+ (* 3 face-index) 0)))
                    (b (aref all-faces (+ (* 3 face-index) 1)))
                    (c (aref all-faces (+ (* 3 face-index) 2))))
                   (when (flet ((edge-in-hull-p (v1 v2)
                                  (let ()
                                    #+TODO (loop :for (facet-centroid . facet-normal) :in (facet-normals hull)
                                                 :thereis
                                        ; (minusp (v. facet-normal (v- vertex facet-centroid)))
                                        ; (not (> (v. facet-normal (v- vertex facet-centroid)) .0001))
                                                    (< (v. (the dvec3 facet-normal) (v- vertex (the dvec3 facet-centroid))) eps))
                                    (or (vertex-in-hull-p v1 hull)
                                        (vertex-in-hull-p v2 hull)
                                        (vertex-in-hull-p (v* (v+ v1 v2) .5d0) hull)))

                                  #+old (and (not (let ((result (loop :for i :below (length hull-faces/global) :by 3
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
                                                   :with bad = NIL
                                                   :with sample-count = 8
                                                   :repeat sample-count
                                                   :with d = (v/ (v- v2 v1) sample-count)
                                                   :for v = (v+ v1 (v* d (/ .5d0 sample-count))) :then (v+ v d)
                                                   :for in? = (vertex-in-hull-p v hull)
                                                   :when in?
                                                   :do (d  "; [~A ~A]~&;~2@T~A - ~A~&;~2@Tat ~A~%"
                                                           vi1 vi2 v1 v2 v)
                                                       (setf bad T)
                                        ; :do (push (cons v (if in? :bad :good)) (problem hull))
                                        ; :thereis (vertex-in-hull-p v hull all-vertices)
                                                   :finally (return bad)))))
                           (let ((v1 (manifolds:v all-vertices a))
                                 (v2 (manifolds:v all-vertices b))
                                 (v3 (manifolds:v all-vertices c)))
                             (or (edge-in-hull-p v1 v2)
                                 (edge-in-hull-p v2 v3)
                                 (edge-in-hull-p v3 v1))))
                     (return T))))
            ; (vector-push-extend (cons hull overlapping) *overlap-counts*)
            ))
    (setf (problem hull) :edge)
    T))

#+unused (defun face-matches-facet-p (h1 h2 h3 v1 v2 v3
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
    (or                    ; (face-matches-facet-p h1 h2 h3 v1 v2 v3)
     (multiple-value-bind (intersectp constellation contact)
         (triangles-intersect-p h1 h2 h3 v1 v2 v3 :threshold threshold)
       (declare (ignore intersectp))
       (and (eq constellation :coplanar) (eq contact :penetrating))))))

(defun normals-matching-p (patch hull context)
  (or (<= (length (vertices hull)) (* 3 4)) ; HACK to work around degenerate 2D case
      (let ((all-vertices (context-vertices context))
            (all-faces (context-faces context))
            (hull-vertices (vertices hull))
            (hull-faces (facets hull)))
        (;; space:do-all (face-info (context-face-index context) T)
         space:do-overlapping (face-info (context-face-index context) hull T)
         ;; manifolds:do-faces (ia ib ic all-faces T)
          (let* #+no ()
            ((face-index (face-info-index face-info))
             (face-normal (face-info-normal face-info))
             (ia (aref all-faces (+ (* 3 face-index) 0)))
             (ib (aref all-faces (+ (* 3 face-index) 1)))
             (ic (aref all-faces (+ (* 3 face-index) 2)))
             (a (manifolds:v all-vertices ia))
             (b (manifolds:v all-vertices ib))
             (c (manifolds:v all-vertices ic))
             ; (face-normal (vunit (vc (v- b a) (v- c a))))
             )
            (when (not (loop for i below (/ (length hull-faces) 3)
                             for (facet-centroid . facet-normal) in (facet-normals hull)
                             always (if (and (< (abs (- -1 (v. face-normal (the dvec3 facet-normal)))) .0001)
                                             (face-overlaps-or-matches-facet-p
                                              hull-vertices hull-faces i a b c :threshold .0001 :hull hull))
                                        (let ( ; (*debug-output* T)
                                              )
                                          (d "[~a ~a ~a] possibly intersects facet ~a~%"
                                             ia ib ic i)
                                          (let ()
                                            (when *debug-visualizations*
                                              (debug-line* a b :face-edge hull)
                                              (debug-line* b c :face-edge hull)
                                              (debug-line* c a :face-edge hull)
                                              (let ((face-centroid (v/ (v+ a b c) 3)))
                                                (push (cons face-centroid :face-centroid) (annotations hull))
                                                (debug-line face-centroid (v* face-normal .3) :face-normal hull :sample-count 13)))
                                            (d "~2@tface normal ~a~&~2@tfacet normal ~a~&~2@t=> ~a~%"
                                               face-normal
                                               facet-normal
                                               (v. face-normal facet-normal))
                                            (cond (T ; (< (abs (- -1 (v. face-normal facet-normal))) .0001)
                                                   (when *debug-visualizations*
                                                     (debug-line facet-centroid face-normal :bad-normal hull)
                                                     (debug-line facet-centroid facet-normal :bad-normal hull))
                                                   NIL)
                                                  (T
                                                   T))))
                                         T)))
              (setf (problem hull) :normals)
              (return NIL)))))))

(defun normals-matching-p/overlap-per-facet (patch hull context)
  (declare (optimize speed))
  (or (<= (length (vertices hull)) (* 3 4)) ; HACK to work around degenerate 2D case
      (let ((face-container (context-face-index context))
            (all-vertices   (context-vertices context))
            (all-faces      (context-faces context))
            (hull-vertices  (vertices hull))
            (hull-faces     (facets hull)))
        (loop for i below (/ (length hull-faces) 3)
              for (facet-centroid . facet-normal) in (facet-normals hull)
              for facet-bounding-box = (face-bounding-box hull-vertices hull-faces i)
              always (space:do-all (face-info face-container T)
                       ;; space:do-overlapping (face-info face-container facet-bounding-box T)
                        (let ((face-normal (face-info-normal face-info)))
                          (declare (type dvec3 face-normal))
                          (when (< (abs (- -1 (v. face-normal (the dvec3 facet-normal)))) .0001)
                            (let* ((face-index (face-info-index face-info))
                                   (ia (aref all-faces (+ (* 3 face-index) 0)))
                                   (ib (aref all-faces (+ (* 3 face-index) 1)))
                                   (ic (aref all-faces (+ (* 3 face-index) 2)))
                                   (a (manifolds:v all-vertices ia))
                                   (b (manifolds:v all-vertices ib))
                                   (c (manifolds:v all-vertices ic))
                                   ;; (face-normal (vunit (vc (v- b a) (v- c a))))
                                   )
                              (when (face-overlaps-or-matches-facet-p
                                     hull-vertices hull-faces i a b c :hull hull)
                                (let (  ; (*debug-output* T)
                                      )
                                  (d "[~A ~A ~A] possibly intersects facet ~A~%"
                                     ia ib ic i)
                                  (let ()
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
                                    (when *debug-visualizations*
                                      (debug-line facet-centroid face-normal :bad-normal hull)
                                      (debug-line facet-centroid facet-normal :bad-normal hull))
                                    (return NIL))))))))))))

(defun find-boundary-constraint-bar (hull context)
  ;; We get a permuted version of the face indices but the order
  ;; (mod 3) and therefore the normal direction should be the same
  ;; as for the original face (and the normal direction does not
  ;; matter for boundary constraint bars anyway).
  ;;
  ;; The paper does not say how to scale constrain bars :(
  (space:do-overlapping (info (context-boundary-edge-index context) hull NIL)
    (let ((bar1 (boundary-bar-info-bar1 info))
          (bar2 (boundary-bar-info-bar2 info)))
      (when (or (vertex-in-hull-p bar1 hull)
                (vertex-in-hull-p bar2 hull)
                (vertex-in-hull-p (v* (v+ bar1 bar2) .5d0) hull)))))
  #+old (let ((all-vertices (context-vertices context))
        (all-faces (context-faces context))
        (hull-faces (facets hull)))
    (loop for edge across (context-boundary-edges context)
          for a = (manifolds:start edge)
          for b = (manifolds:end edge)
          for c = (manifolds:opposite edge)
          for v1 = (manifolds:v all-vertices a)
          for v2 = (manifolds:v all-vertices b)
          for v3 = (manifolds:v all-vertices c)
          for edge-center = (v/ (v+ v1 v2) 2)
          for face-normal = (vunit (vc (v- v2 v1) (v- v3 v1)))
          for offset = (v* face-normal .1)
          for bar1 = (v- edge-center offset)
          for bar2 = (v+ edge-center offset)
          do (when *debug-visualizations*
               #+ no (let ((c (v/ (v+ v1 v2 v3) 3))
                           (n (vunit (vc (v- v2 v1) (v- v3 v1)))))
                       (debug-line c n :face-normal hull))
               (debug-line* bar1 bar2 :face-normal hull))
          thereis (or (vertex-in-hull-p bar1 hull)
                      (vertex-in-hull-p bar2 hull)
                      (vertex-in-hull-p (v* (v+ bar1 bar2) .5d0) hull)))))

(defun find-negative-side-touching-triangle (hull)
  NIL)

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
