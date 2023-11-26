(in-package #:org.shirakumo.fraf.convex-covering)

;;; Vertex position index
;;;
;;; Lookup index in global vertex array given the vertex coordinates.

(defun make-vertex-index ()
  (make-hash-table :test #'equalp))

(defun vertex-position (vertex index)
  (check-type vertex dvec3)
  (the (or null manifolds:u32) (gethash vertex index)))

(defun (setf vertex-position) (new-value vertex index)
  (check-type vertex dvec3)
  (setf (gethash vertex index) new-value))

(defun index-vertex-positions (vertices)
  (check-type vertices manifolds:vertex-array)
  (loop with index = (make-vertex-index)
        for i below (/ (length vertices) 3)
        for vertex = (manifolds:v vertices i)
        do (setf (vertex-position vertex index) i)
        finally (return index)))

;;; Spatial indices

(defun make-spatial-index ()
  (ecase :kd-tree
    (:grid    (org.shirakumo.fraf.trial.space.grid3:make-grid .3 :bsize (vec3 2 2 2)))
    (:kd-tree (org.shirakumo.fraf.trial.space.kd-tree:make-kd-tree :dimensions 3))))

;;; Vertex index
;;;
;;; This allows efficiently finding vertices that are contained in a
;;; given axis-aligned box.

(defstruct (vertex-info
            (:constructor make-vertex-info (index location))
            (:predicate NIL)
            (:copier NIL))
  (index    (error "required") :type manifolds:u32 :read-only T)
  (location (error "required") :type vec3          :read-only T))

(defmethod space:location ((object vertex-info))
  (vertex-info-location object))

(defmethod space:bsize ((object vertex-info))
  (load-time-value (vec 0 0 0)))

(defun index-vertices (vertices)
  (check-type vertices manifolds:vertex-array)
  (loop with index = (make-spatial-index)
        for i below (/ (length vertices) 3)
        for vertex = (manifolds:v vertices i)
        for info = (make-vertex-info i (vec vertex))
        do (space:enter info index)
        finally (return index)))

;;; Edge index
;;;
;;; This allows efficiently finding edges that intersect a given
;;; axis-aligned box.

(defstruct (edge-info
            (:constructor make-edge-info (key vertex1 vertex2 center size/2))
            (:predicate NIL)
            (:copier NIL))
  ;; The key is used to identify the edge when testing whether it is
  ;; part of a hull. The key consists of the indices into the mesh
  ;; vertex array encoded as (logior (ash IS 32) IL) where IS is the
  ;; smaller of the two indices.
  (key     (error "required") :type (unsigned-byte 64) :read-only T)
  (vertex1 (error "required") :type dvec3 :read-only T)
  (vertex2 (error "required") :type dvec3 :read-only T)
  (center  (error "required") :type vec3  :read-only T)
  (size/2  (error "required") :type vec3  :read-only T))

(defmethod space:location ((object edge-info))
  (edge-info-center object))

(defmethod space:bsize ((object edge-info))
  (edge-info-size/2 object))

(defun index-edges (vertices faces)
  (check-type vertices manifolds:vertex-array)
  (check-type faces manifolds:face-array)
  (let ((index (make-spatial-index)))
    (flet ((add-edge (key vertex1 vertex2)
             (let ((min (vmin vertex1 vertex2))
                   (max (vmax vertex1 vertex2)))
               (multiple-value-bind (center size/2)
                   (center-and-size-from-min-and-max min max)
                 (let ((info (make-edge-info key vertex1 vertex2 center size/2)))
                   (space:enter info index))))))
      (loop for edge across (manifolds:edge-list faces)
            for i1 = (manifolds:start edge)
            for i2 = (manifolds:end edge)
            for key = (if (< i1 i2) ; TODO(jmoringe): make a function
                          (logior (ash i1 32) i2)
                          (logior (ash i2 32) i1))
            for v1 = (manifolds:v vertices i1)
            for v2 = (manifolds:v vertices i2)
            do (add-edge key v1 v2)))
    (space:reoptimize index) ; TODO(jmoringe): use `space:enter-all' instead
    index))

;;; Face index
;;;
;;; This allows efficiently finding faces that intersect a given
;;; axis-aligned box.

(defstruct (face-info
            (:constructor make-face-info (index center size/2 normal))
            (:predicate NIL)
            (:copier NIL))
  (index  (error "required") :type manifolds:u32 :read-only T)
  (center (error "required") :type vec3          :read-only T)
  (size/2 (error "required") :type vec3          :read-only T)
  (normal (error "required") :type dvec3         :read-only T))

(defmethod space:location ((object face-info))
  (face-info-center object))

(defmethod space:bsize ((object face-info))
  (face-info-size/2 object))

(defun index-faces (vertices faces)
  (check-type vertices manifolds:vertex-array)
  (check-type faces manifolds:face-array)
  (let ((index (make-spatial-index)))
    (loop for i below (/ (length faces) 3)
          for (center size/2) = (multiple-value-list ; SBCL optimizes this away
                                 (face-bounding-box vertices faces i))
          for normal = (vunit (manifolds:face-normal vertices faces i))
          for info = (make-face-info i center size/2 normal)
          do (space:enter info index))
    (space:reoptimize index) ; TODO(jmoringe): use `space:enter-all' instead
    index))

;;; Boundary edge index
;;;
;;; This allows efficiently finding boundary constraint bars that
;;; intersect a given axis-aligned box.

(defstruct (boundary-bar-info
            (:constructor make-boundary-bar-info (bar1 bar2)))
  (bar1 (error "required") :type dvec3 :read-only t)
  (bar2 (error "required") :type dvec3 :read-only t))

(defmethod space:location ((object boundary-bar-info))
  (let ((a (boundary-bar-info-bar1 object))
        (b (boundary-bar-info-bar2 object)))
    ;; TODO avoid conversion
    (vec (v/ (v+ a b) 2))))

(defmethod space:bsize ((object boundary-bar-info))
  (let ((a (boundary-bar-info-bar1 object))
        (b (boundary-bar-info-bar2 object)))
    ;; TODO avoid conversion
    (vec (v/ (vabs (v- a b)) 2))))

(defun compute-boundary-bar (all-vertices boundary-edge)
  (let* ((a (manifolds:start boundary-edge))
         (b (manifolds:end boundary-edge))
         (c (manifolds:opposite boundary-edge))
         (v1 (manifolds:v all-vertices a))
         (v2 (manifolds:v all-vertices b))
         (v3 (manifolds:v all-vertices c))
         (edge-center (v/ (v+ v1 v2) 2))
         (face-normal (vunit (vc (v- v2 v1) (v- v3 v1))))
         ;; TODO(jmoringe): the paper doesn't say how to scale this
         (offset (v* face-normal .1))
         (bar1 (v- edge-center offset))
         (bar2 (v+ edge-center offset)))
    (values bar1 bar2)))

(defun index-boundary-edges (all-vertices boundary-edges)
  (let ((index (org.shirakumo.fraf.trial.space.kd-tree:make-kd-tree :dimensions 3)))
    (loop for boundary-edge across boundary-edges
          for (bar1 bar2) = (multiple-value-list (compute-boundary-bar
                                                  all-vertices boundary-edge))
          for info = (make-boundary-bar-info bar1 bar2)
          do (space:enter info index))
    (space:reoptimize index) ; TODO enter-all
    index))

;;; Context
;;;
;;; A context instance stores various pieces of information that are
;;; needed for computing the convex decomposition such as vertices,
;;; faces and (boundary) edges of the mesh being decomposed as well as
;;; spatial index structures.

(defstruct (context
            (:constructor %make-context (vertices faces cost-function
                                         edge-tolerance normals-tolerance
                                         vertex-position-index
                                         vertex-index
                                         edge-index
                                         face-index
                                         boundary-edges
                                         boundary-edge-index))
            (:predicate nil)
            (:copier nil))
  ;; Mesh
  (vertices              (error "required") :type (manifolds:vertex-array manifolds:f64) :read-only T)
  (faces                 (error "required") :type manifolds:face-array :read-only T)
  (boundary-edges        (error "required") :read-only T) ; TODO types
  ;; Cost function and other parameters
  (cost-function         (error "required") :type function :read-only T)
  (edge-tolerance        (error "required") :type double-float :read-only T)
  (normals-tolerance     (error "required") :type double-float :read-only T)
  ;; Index structures
  (vertex-position-index (error "required") :read-only T)
  (vertex-index          (error "required") :read-only T)
  (edge-index            (error "required") :read-only T)
  (face-index            (error "required") :read-only T)
  (boundary-edge-index   (error "required") :read-only T))

(defun make-context (vertices faces cost-function
                     &key (edge-tolerance    1d-4)
                          (normals-tolerance 1d-4))
  (lparallel:plet ((vertex-position-index (index-vertex-positions vertices))
                   (vertex-index (index-vertices vertices))
                   (edge-index (index-edges vertices faces))
                   (face-index (index-faces vertices faces))
                   ((boundary-edges boundary-edge-index)
                    (let ((boundary-edges (manifolds:boundary-list faces)))
                      (values boundary-edges
                              (index-boundary-edges vertices boundary-edges)))))
    (%make-context vertices faces cost-function
                   edge-tolerance normals-tolerance
                   vertex-position-index
                   vertex-index
                   edge-index
                   face-index
                   boundary-edges
                   boundary-edge-index)))
