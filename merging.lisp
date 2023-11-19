(in-package #:org.shirakumo.fraf.convex-covering)

;;; Merge costs

(deftype merge-cost-function ()
  '(function (patch patch patch) (values (or null double-float) &optional nil)))

(defun coerce-to-cost-function (thing)
  (typecase thing
    (symbol (fdefinition thing))
    (function thing)
    ((cons t null) (coerce-to-cost-function (first thing)))
    (list (apply #'combine (mapcar #'coerce-to-cost-function thing)))))

(declaim (ftype merge-cost-function /compactness patch-size-symmetry))

(defun /compactness (merged-patch patch1 patch2)
  (declare (ignore patch1 patch2))
  (/ (patch-compactness merged-patch)))

(defun patch-size-symmetry (merged-patch patch1 patch2)
  (declare (optimize speed (safety 1))
           (ignore merged-patch))
  (flet ((patch-hull-facet-count (patch)
           (let ((hull (patch-hull patch)))
             (if hull
                 (length (hull-facets hull))
                 1))))
    (coerce (log (1+ (abs (- (patch-hull-facet-count patch1)
                             (patch-hull-facet-count patch2)))))
            'double-float)))

(defun make-patch-size-limit (limit)
  (declare (type manifolds:u32 limit))
  (lambda (merged-patch patch1 patch2)
    (declare (optimize speed (safety 1))
             (ignore patch1 patch2))
    (let ((hull (patch-hull merged-patch)))
      (cond ((null hull)
             0)
            ((< (length (hull-facets hull)) limit)
             0)
            (T
             NIL)))))

(let ((cache (make-hash-table :test #'equal)))
  (defun combine (&rest functions)
    (or (gethash functions cache)
        (setf (gethash functions cache)
              (let ((variables '()))
                (labels ((call (function)
                           (let ((variable (gensym "COMPONENT")))
                             (push variable variables)
                             `(,variable (funcall ,function merged-patch patch1 patch2)))))
                  (compile nil `(lambda (merged-patch patch1 patch2)
                                  (let (,@(mapcar #'call functions))
                                    (when (and ,@variables)
                                      (+ ,@variables)))))))))))

;;; Merge validity

(defun valid-patch-p (patch context)
  (let ((valid
          (handler-case
              (let ((hull (patch-hull patch)))
                (and ; (not (find-vertex-in-hull hull context))
                     ; (not (find-edge-in-hull hull context))
                     (not (find-edge-in-hull/edge-index hull context))
                     (normals-matching-p patch hull context)
                     (not (find-boundary-constraint-bar hull context))
                     (not (find-negative-side-touching-triangle hull))))
            (floating-point-invalid-operation (error)
              (format *error-output* "~a~%" error)
              NIL))))
    (d "; hull ~a is ~:[invalid~;valid~]"
       (hull-global-faces (patch-hull patch)) valid)
    valid))

;;; Individual criteria

(defun vertex-in-hull-p (vertex hull &key (eps -.01))
  (declare (type hull hull)
           (type dvec3 vertex))
  (destructuring-bind (center . size/2) (bounding-box hull)
    (declare (type vec3 center size/2))
    (and (v<= (v- center size/2) (vec3 (vx vertex) (vy vertex) (vz vertex)) (v+ center size/2)) ; TODO avoid conversion
         (loop for (facet-centroid . facet-normal) in (facet-normals hull)
               always
                                        ; (minusp (v. facet-normal (v- vertex facet-centroid)))
                                        ; (not (> (v. facet-normal (v- vertex facet-centroid)) .0001))
                  (< (v. (the dvec3 facet-normal) (v- vertex (the dvec3 facet-centroid))) eps)))))

;;;
(defun vertex-in-hull-p* (vertex hull &key (eps 1e-6))
  (declare (type hull hull)
           (type dvec3 vertex))
  (destructuring-bind (center . size/2) (bounding-box hull)
    (declare (type vec3 center size/2))
    (let* ((abs-eps (abs eps))
           (offset  (vec abs-eps abs-eps abs-eps)))
      (and (v<= (v- center size/2 offset)
                (vec3 (vx vertex) (vy vertex) (vz vertex))
                (v+ center size/2 offset)) ; TODO avoid conversion
           (loop for (facet-centroid . facet-normal) in (facet-normals hull)
                 for depth = (v. (the dvec3 facet-normal) (v- vertex (the dvec3 facet-centroid)))
                 do (when *debug-output*
                      (format *trace-output* "Depth ~A~%" depth))
                 always (<= depth eps))))))

;;; This version assumes that VERTEX is within the bounding box of
;;; HULL and thus only checks VERTEX against all hull facets.
(defun vertex-in-hull-p** (vertex hull &key (eps -.01))
  (declare (type hull hull)
           (type vec3 vertex))
  (let ((vertex (dvec vertex)))
    (loop for (facet-centroid . facet-normal) in (facet-normals hull)
          always (< (v. (the dvec3 facet-normal) (v- vertex (the dvec3 facet-centroid))) eps))))

(defun find-vertex-in-hull (hull context)
  (declare (type hull hull)
           (optimize speed (safety 1) (debug 1)))
  (let ((hull-faces (hull-global-faces hull)))
    (space:do-overlapping (vertex (context-vertex-index context) hull nil)
      (let ((index    (vertex-info-index vertex))
            (position (vertex-info-location vertex)))
        (when (and *debug-visualizations*
                   (vertex-in-hull-p** position hull))
          (push (cons vertex :bad) (hull-annotations hull)))
        (when (vertex-in-hull-p** position hull)
          ;; TODO(jmoringe): remove later; vertices that are part of
          ;; the hull can never be properly inside the hull
          (assert (not (find index hull-faces)))
          (return T)))))
  #+no (let* ((all-vertices (context-vertices context))
              (result (loop with hull-faces = (hull-global-faces hull)
                            for vertex-index from 0 below (/ (length all-vertices) 3)
                            for vertex = (manifolds:v all-vertices vertex-index)
                            when (and *debug-visualizations*
                                      (not (find vertex-index hull-faces))
                                      (vertex-in-hull-p vertex hull))
                            do (push (cons vertex :bad) (hull-annotations hull))
                            thereis (and (vertex-in-hull-p vertex hull)
                                         (not (find vertex-index hull-faces))))))
         (when result
           (setf (hull-problem hull) :vertex))
         result))

(defun find-edge-in-hull (hull context)
  (declare (type hull hull))
  (when (let ((all-vertices (context-vertices context))
              (all-faces (context-faces context))
              (hull-vertices (hull-vertices hull))
              (hull-facets (hull-facets hull)))
          (unwind-protect
               ;; TODO(jmoringe): could do overlapping queries for individual hull facets
               (space:do-overlapping (face-info (context-face-index context) hull NIL)
                 (let* ((face-index (face-info-index face-info))
                        (a (aref all-faces (+ (* 3 face-index) 0)))
                        (b (aref all-faces (+ (* 3 face-index) 1)))
                        (c (aref all-faces (+ (* 3 face-index) 2))))
                   (when (flet ((edge-in-hull-p (v1 v2)
                                  (manifolds:do-faces (ai bi ci hull-facets nil)
                                    (let ((a (manifolds:v hull-vertices ai))
                                          (b (manifolds:v hull-vertices bi))
                                          (c (manifolds:v hull-vertices ci)))
                                      (multiple-value-bind (result constellation contact)
                                          (line-intersects-triangle-p a b c v1 v2 :threshold 1d-4)
                                        (when (and result
                                                   (case contact
                                                     (:touching
                                                      ;; The edge is just touching the hull facet. Could be touching
                                                      ;; at a vertex of the hull facet or in the (2d) interior of the
                                                      ;; hull facet. For each end point of the edge, generate a sample
                                                      ;; point that is offset towards the center of the edge from the
                                                      ;; end point and test whether that point is properly inside the
                                                      ;; hull. Chose the offset large enough so handle an edge that is
                                                      ;; almost coplanar to the hull facet (assume 1° incidence angle
                                                      ;; as the worst case).
                                                      (let* ((threshold -1e-2)
                                                             (diff      (v- v2 v1))
                                                             (delta     (v* diff (min .5 (* 1.5
                                                                                            (/ (sin (/ pi 180)))
                                                                                            (abs threshold)
                                                                                            (/ (vlength diff)))))))
                                                        (or (vertex-in-hull-p* (v+ v1 delta) hull :eps threshold)
                                                            (vertex-in-hull-p* (v- v2 delta) hull :eps threshold))))
                                                     (:penetrating
                                                      ;; The edge is penetrating the hull face. If the edge and the hull
                                                      ;; face are coplanar, this is not a problem since the two may just
                                                      ;; be parts of different triangulations of a non-triangle mesh face.
                                                      ;; Other forms of penetrating contact make the hull invalid.
                                                      (not (eq constellation :coplanar))
                                                      #+no (or
                                                            (vertex-in-hull-p* v1 hull :eps -1e-4)
                                                            (vertex-in-hull-p* v2 hull :eps -1e-4)
                                                            (vertex-in-hull-p* (v* (v+ v1 v2) .5) hull :eps -1e-4)
                                                            #+no (and (not (eq constellation :coplanar))
                                                                      (or (not (vertex-in-hull-p* v1 hull :eps 1e-4))
                                                                          (not (vertex-in-hull-p* v2 hull :eps 1e-4))))
                                                            (not (vertex-in-hull-p* v1 hull :eps 1e-4))
                                                            (not (vertex-in-hull-p* v2 hull :eps 1e-4)))))
                                                   #+no (or
                                                         #+no (and (eq contact :touching)
                                                                   (or (and (multiple-value-bind (result contact)
                                                                                (point-on-triangle-p a b c v1)
                                                                              result ; (and result (not (eq contact :touching)))
                                                                              )
                                                                            (vertex-in-hull-p* v2 hull :eps -1e-2))
                                                                       (and (multiple-value-bind (result contact)
                                                                                (point-on-triangle-p a b c v2)
                                                                              result ; (and result (not (eq contact) :touching))
                                                                              )
                                                                            (vertex-in-hull-p* v1 hull :eps -1e-2))))))
                                          #+no (let ((*debug-output* t))
                                                 (line-intersects-triangle-p a b c v1 v2 :threshold 1d-4)
                                                 (vertex-in-hull-p* (v* (v+ v1 v2) .5) hull :eps -1e-2))
                                          (when (member :edges *debug-visualizations*)
                                            (push (list :triangle a b c :color '(1 .5 0))
                                                  (hull-annotations hull))
                                            (push (list :line v1 v2 :color '(1 0 .5))
                                                  (hull-annotations hull)))
                                          #+no (format *trace-output* "~A ~A ~A v1 ~A v2 ~A~%"
                                                       result constellation contact
                                                       (vertex-in-hull-p* v1 hull :eps -1e-2)
                                                       (vertex-in-hull-p* v2 hull :eps -1e-2))
                                          (return T)))))))
                           (let ((v1 (manifolds:v all-vertices a))
                                 (v2 (manifolds:v all-vertices b))
                                 (v3 (manifolds:v all-vertices c)))
                             #+no (when (member :edges *debug-visualizations*)
                                    (push (list :triangle v1 v2 v3 :color '(.5 .5 1))
                                          (hull-annotations hull)))
                             (or (edge-in-hull-p v1 v2)
                                 (edge-in-hull-p v2 v3)
                                 (edge-in-hull-p v3 v1))))
                     (return T))))))
    (setf (hull-problem hull) :edge)
    T))

(defun find-edge-in-hull/edge-index (hull context)
  (declare (type hull hull))
  (let ((hull-vertices (hull-vertices hull))
        (hull-facets (hull-facets hull)))
    ;; TODO(jmoringe): could do overlapping queries for
    ;; individual hull facets
    #+no (format *trace-output* "~&; Hull ~X, facets ~:D~%"
            (sxhash hull)
            (/ (length hull-facets) 3))
    (flet ((check-pair (a b c v1 v2)
             (multiple-value-bind (result constellation contact)
                 (line-intersects-triangle-p a b c v1 v2 :threshold 1d-4)
               (when (member :edges *debug-visualizations*)
                 (push (list :triangle a b c :color '(.4 .4 .2))
                       (hull-annotations hull)))
               (when (and result
                          (case contact
                            (:touching
                             (when (member :edges *debug-visualizations*)
                               (push (list :line v1 v2 :color '(.2 .2 .4))
                                     (hull-annotations hull)))
                             ;; The edge is just touching the hull facet. Could be touching
                             ;; at a vertex of the hull facet or in the (2d) interior of the
                             ;; hull facet. For each end point of the edge, generate a sample
                             ;; point that is offset towards the center of the edge from the
                             ;; end point and test whether that point is properly inside the
                             ;; hull. Chose the offset large enough so handle an edge that is
                             ;; almost coplanar to the hull facet (assume 1° incidence angle
                             ;; as the worst case).
                             (let* ((threshold -1e-2)
                                    (diff      (v- v2 v1))
                                    (delta     (v* diff (min .5 (* 1.5
                                                                   (/ (sin (/ pi 180)))
                                                                   (abs threshold)
                                                                   (/ (vlength diff)))))))
                               (or (vertex-in-hull-p* (v+ v1 delta) hull :eps threshold)
                                   (vertex-in-hull-p* (v- v2 delta) hull :eps threshold))))
                            (:penetrating
                             (when (member :edges *debug-visualizations*)
                               (push (list :line v1 v2 :color '(.2 .4 .4))
                                     (hull-annotations hull)))
                             ;; The edge is penetrating the hull face. If the edge and the hull
                             ;; face are coplanar, this is not a problem since the two may just
                             ;; be parts of different triangulations of a non-triangle mesh face.
                             ;; Other forms of penetrating contact make the hull invalid.
                             (not (eq constellation :coplanar)))))
                 (when (member :edges *debug-visualizations*)
                   (push (list :triangle a b c :color '(1 .5 0))
                         (hull-annotations hull))
                   (push (list :line v1 v2 :color '(1 0 .5))
                         (hull-annotations hull)))
                 (return-from find-edge-in-hull/edge-index T)))))
      (if (< (length hull-facets) (* 3 10))
          (space:do-overlapping (edge-info (context-edge-index context) hull NIL)
            (let ((v1 (edge-info-vertex1 edge-info))
                  (v2 (edge-info-vertex2 edge-info)))
              (manifolds:do-faces (ai bi ci hull-facets) ; TODO use bounding box when few facets
                (let* ((a   (manifolds:v hull-vertices ai))
                       (b   (manifolds:v hull-vertices bi))
                       (c   (manifolds:v hull-vertices ci)))
                  (check-pair a b c v1 v2)))))
          (manifolds:do-faces (ai bi ci hull-facets NIL) ; TODO use bounding box when few facets
            (let* ((a   (manifolds:v hull-vertices ai))
                   (b   (manifolds:v hull-vertices bi))
                   (c   (manifolds:v hull-vertices ci))
                   (min (vec (v- (vmin a b c) (dvec .001 .001 .001))))
                   (max (vec (v+ (vmax a b c) (dvec .001 .001 .001)))))
              (let ((region (space::%region (org.shirakumo.fraf.math.vectors::varr3 min) (v- max min))))
                (space:do-overlapping (edge-info (context-edge-index context) region)
                  (let ((v1 (edge-info-vertex1 edge-info))
                        (v2 (edge-info-vertex2 edge-info)))
                    (check-pair a b c v1 v2))))))))))

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
  (declare (type hull hull))
  (or (hull-flat-p hull)
      (let ((all-vertices (context-vertices context))
            (all-faces (context-faces context))
            (hull-vertices (hull-vertices hull))
            (hull-faces (hull-facets hull)))
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
           (when (member :normals *debug-visualizations*)
             (push (list :triangle a b c) (hull-annotations hull))
             (let ((face-centroid (v/ (v+ a b c) 3)))
               (push (list :line face-centroid (v+ face-centroid (v* face-normal .3)))
                     (hull-annotations hull))))

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
                                           (when (member :normals *debug-visualizations*)
                                        ; (debug-line* a b :face-edge hull)
                                        ; (debug-line* b c :face-edge hull)
                                        ; (debug-line* c a :face-edge hull)
                                             (push (list :triangle a b c :color '(.5 1 .5))
                                                   (hull-annotations hull))
                                             (let ((face-centroid (v/ (v+ a b c) 3)))
                                               (push (list :point face-centroid :color '(.5 1 .5)) (hull-annotations hull))
                                        ; (debug-line face-centroid (v* face-normal .3) :face-normal hull :sample-count 13)
                                               ))
                                           (d "~2@tface normal ~a~&~2@tfacet normal ~a~&~2@t=> ~a~%"
                                              face-normal
                                              facet-normal
                                              (v. face-normal facet-normal))
                                           (cond (T ; (< (abs (- -1 (v. face-normal facet-normal))) .0001)
                                                  (when (member :normals *debug-visualizations*)
                                                    (push (list :line facet-centroid face-normal :color '(.8 .8 0))
                                                          (hull-annotations hull))
                                                    (push (list :line facet-centroid facet-normal :color '(.8 .8 0))
                                                          (hull-annotations hull)))
                                                  NIL)
                                                 (T
                                                  T))))
                                       T)))
             (setf (hull-problem hull) :normals)
             (return NIL)))))))

(defun normals-matching-p/overlap-per-facet (patch hull context)
  (declare (type hull hull)
           (optimize speed))
  (or (<= (length (hull-vertices hull)) (* 3 4)) ; HACK to work around degenerate 2D case
      (let ((face-container (context-face-index context))
            (all-vertices   (context-vertices context))
            (all-faces      (context-faces context))
            (hull-vertices  (hull-vertices hull))
            (hull-faces     (hull-facets hull)))
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
                                        (push (cons face-centroid :face-centroid) (hull-annotations hull))
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
  (declare (type (manifolds:vertex-array manifolds:f64) all-vertices)
           (type patch patch))
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
